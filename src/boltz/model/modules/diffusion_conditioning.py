from __future__ import annotations

import torch
from torch import nn
from torch.nn import Module

from boltz.model.modules.encodersv2 import (
    AtomEncoder,
    PairwiseConditioning,
)


class DiffusionConditioning(Module):
    def __init__(
        self,
        token_s: int,
        token_z: int,
        atom_s: int,
        atom_z: int,
        atoms_per_window_queries: int = 32,
        atoms_per_window_keys: int = 128,
        atom_encoder_depth: int = 3,
        atom_encoder_heads: int = 4,
        token_transformer_depth: int = 24,
        token_transformer_heads: int = 8,
        atom_decoder_depth: int = 3,
        atom_decoder_heads: int = 4,
        atom_feature_dim: int = 128,
        conditioning_transition_layers: int = 2,
        use_no_atom_char: bool = False,
        use_atom_backbone_feat: bool = False,
        use_residue_feats_atoms: bool = False,
    ) -> None:
        super().__init__()

        self.pairwise_conditioner = PairwiseConditioning(
            token_z=token_z,
            dim_token_rel_pos_feats=token_z,
            num_transitions=conditioning_transition_layers,
        )

        self.atom_encoder = AtomEncoder(
            atom_s=atom_s,
            atom_z=atom_z,
            token_s=token_s,
            token_z=token_z,
            atoms_per_window_queries=atoms_per_window_queries,
            atoms_per_window_keys=atoms_per_window_keys,
            atom_feature_dim=atom_feature_dim,
            structure_prediction=True,
            use_no_atom_char=use_no_atom_char,
            use_atom_backbone_feat=use_atom_backbone_feat,
            use_residue_feats_atoms=use_residue_feats_atoms,
        )

        self.atom_enc_proj_z = nn.ModuleList()
        for _ in range(atom_encoder_depth):
            self.atom_enc_proj_z.append(
                nn.Sequential(
                    nn.LayerNorm(atom_z),
                    nn.Linear(atom_z, atom_encoder_heads, bias=False),
                )
            )

        self.atom_dec_proj_z = nn.ModuleList()
        for _ in range(atom_decoder_depth):
            self.atom_dec_proj_z.append(
                nn.Sequential(
                    nn.LayerNorm(atom_z),
                    nn.Linear(atom_z, atom_decoder_heads, bias=False),
                )
            )

        self.token_trans_proj_z = nn.ModuleList()
        for _ in range(token_transformer_depth):
            self.token_trans_proj_z.append(
                nn.Sequential(
                    nn.LayerNorm(token_z),
                    nn.Linear(token_z, token_transformer_heads, bias=False),
                )
            )

    def forward(
        self,
        s_trunk,  # Float['b n ts']
        z_trunk,  # Float['b n n tz']
        relative_position_encoding,  # Float['b n n tz']
        feats,
    ):
        z = self.pairwise_conditioner(
            z_trunk,
            relative_position_encoding,
        )

        # cdr3_token_trans_delta is returned alongside token_trans_bias for time-varying
        # beta support (L+.4). None when CDR3 features are absent.
        cdr3_token_trans_delta = None

        if "cdr3_beta_token" in feats:
            # --- CDR3 Beta-Scaling v2 (L+.2, L+.3, L+.4) ---
            cdr3_beta_token = feats["cdr3_beta_token"].to(z.device)
            if cdr3_beta_token.dim() == 2:
                cdr3_beta_token = cdr3_beta_token.squeeze(0)

            # L+.2: Asymmetric pair beta via geometric mean.
            # H3-H3 → beta_H3, L3-L3 → beta_L3, H3-L3 → sqrt(beta_H3 * beta_L3)
            beta_i = cdr3_beta_token.unsqueeze(-1)   # [n, 1]
            beta_j = cdr3_beta_token.unsqueeze(-2)   # [1, n]
            cdr3_pair_beta = torch.sqrt(beta_i * beta_j)  # [n, n], 0 for non-CDR3 pairs

            # L+.3: CDR3-antigen interface pair scaling (attract CDR3 toward antigen)
            antigen_pair_beta = torch.zeros_like(cdr3_pair_beta)
            if "cdr3_antigen_token_mask" in feats and "cdr3_antigen_beta" in feats:
                antigen_mask = feats["cdr3_antigen_token_mask"].to(z.device).to(torch.bool)
                antigen_beta_val = feats["cdr3_antigen_beta"].to(z.device).flatten()[0].item()
                if antigen_mask.dim() == 2:
                    antigen_mask = antigen_mask.squeeze(0)
                if abs(antigen_beta_val) > 1e-6 and antigen_mask.any():
                    cdr3_any_mask = cdr3_beta_token > 0  # any CDR3 token
                    # Both CDR3→antigen and antigen→CDR3 directions
                    cdr3_antigen_pair = (
                        (cdr3_any_mask.unsqueeze(-1) & antigen_mask.unsqueeze(-2)) |
                        (antigen_mask.unsqueeze(-1) & cdr3_any_mask.unsqueeze(-2))
                    )  # [n, n]
                    antigen_pair_beta = antigen_beta_val * cdr3_antigen_pair.float()

            # Combined full pair scaling factor [n, n]
            full_pair_scale_2d = 1.0 + cdr3_pair_beta + antigen_pair_beta
            scale_4d = full_pair_scale_2d.unsqueeze(0).unsqueeze(-1)  # [1, n, n, 1]

            # Use fully-scaled z for atom encoder (static at beta_max; approximation for L+.4)
            z_for_atom = z * scale_4d

            # Compute token_trans_bias from BASE z (no CDR3 scaling) so that L+.4 can
            # apply time-varying interpolation: effective_bias = base + t * delta
            token_trans_bias_parts = []
            for layer in self.token_trans_proj_z:
                token_trans_bias_parts.append(layer(z))
            token_trans_bias = torch.cat(token_trans_bias_parts, dim=-1)

            # L+.4: delta = (scale - 1) * base_bias
            # In DiffusionModule: effective = base + t_scale * delta
            # t_scale=1.0 (static) → full CDR3 scaling; t_scale=steering_t → time-varying
            cdr3_token_trans_delta = (scale_4d - 1.0) * token_trans_bias

            q, c, p, to_keys = self.atom_encoder(
                feats=feats,
                s_trunk=s_trunk,
                z=z_for_atom,
            )

        elif "cdr3_token_mask" in feats and "cdr3_beta_value" in feats:
            # Legacy fallback: original uniform CDR3 beta scaling
            cdr3_mask = feats["cdr3_token_mask"].to(z.device).to(torch.bool)
            cdr3_beta = feats["cdr3_beta_value"].to(z.device)
            beta_val = float(cdr3_beta.item())
            if cdr3_mask.dim() == 2:
                cdr3_mask = cdr3_mask.squeeze(0)
            if abs(beta_val) > 1e-6:
                cdr3_mask_i = cdr3_mask.unsqueeze(-1)
                cdr3_mask_j = cdr3_mask.unsqueeze(-2)
                cdr3_pair_mask = (cdr3_mask_i & cdr3_mask_j)
                scaling_tensor = cdr3_pair_mask.unsqueeze(0).unsqueeze(-1).float()
                z = z * (1.0 + beta_val * scaling_tensor)

            q, c, p, to_keys = self.atom_encoder(
                feats=feats,
                s_trunk=s_trunk,
                z=z,
            )
            token_trans_bias_parts = []
            for layer in self.token_trans_proj_z:
                token_trans_bias_parts.append(layer(z))
            token_trans_bias = torch.cat(token_trans_bias_parts, dim=-1)

        else:
            q, c, p, to_keys = self.atom_encoder(
                feats=feats,
                s_trunk=s_trunk,
                z=z,
            )
            token_trans_bias_parts = []
            for layer in self.token_trans_proj_z:
                token_trans_bias_parts.append(layer(z))
            token_trans_bias = torch.cat(token_trans_bias_parts, dim=-1)

        atom_enc_bias = []
        for layer in self.atom_enc_proj_z:
            atom_enc_bias.append(layer(p))
        atom_enc_bias = torch.cat(atom_enc_bias, dim=-1)

        atom_dec_bias = []
        for layer in self.atom_dec_proj_z:
            atom_dec_bias.append(layer(p))
        atom_dec_bias = torch.cat(atom_dec_bias, dim=-1)

        return q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias, cdr3_token_trans_delta

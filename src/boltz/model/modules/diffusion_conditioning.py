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

        # Apply CDR3 beta scaling if enabled
        # Skip static beta scaling when hierarchical steering is active (it handles beta per-step)
        hierarchical_active = (
            "hierarchical_cdr_token_mask" in feats
            and feats["hierarchical_cdr_token_mask"].any()
        )
        if not hierarchical_active and "cdr3_token_mask" in feats and "cdr3_beta_value" in feats:
            cdr3_mask = feats["cdr3_token_mask"].to(z.device).to(torch.bool)
            cdr3_beta = feats["cdr3_beta_value"].to(z.device)
            beta_val = float(cdr3_beta.item())

            # cdr3_mask might have batch dimension [batch, n_tokens] or just [n_tokens]
            if cdr3_mask.dim() == 2:
                # Remove batch dimension
                cdr3_mask = cdr3_mask.squeeze(0)

            if abs(beta_val) > 1e-6:  # Only apply if non-zero
                # Create pair mask: True for (i,j) where both i and j are in CDR3
                # cdr3_mask shape: [n_tokens], z shape: [batch, n_tokens, n_tokens, tz]
                cdr3_mask_i = cdr3_mask.unsqueeze(-1)  # [n_tokens, 1]
                cdr3_mask_j = cdr3_mask.unsqueeze(-2)  # [1, n_tokens]
                cdr3_pair_mask = (cdr3_mask_i & cdr3_mask_j)  # [n_tokens, n_tokens]

                # Apply scaling: z[CDR3] *= (1 + beta)
                # Add batch and feature dimensions for broadcasting: [1, n_tokens, n_tokens, 1]
                scaling_tensor = cdr3_pair_mask.unsqueeze(0).unsqueeze(-1).float()
                scaling_factor = 1.0 + beta_val * scaling_tensor

                z = z * scaling_factor

        q, c, p, to_keys = self.atom_encoder(
            feats=feats,
            s_trunk=s_trunk,  # Float['b n ts'],
            z=z,  # Float['b n n tz'],
        )

        atom_enc_bias = []
        for layer in self.atom_enc_proj_z:
            atom_enc_bias.append(layer(p))
        atom_enc_bias = torch.cat(atom_enc_bias, dim=-1)

        atom_dec_bias = []
        for layer in self.atom_dec_proj_z:
            atom_dec_bias.append(layer(p))
        atom_dec_bias = torch.cat(atom_dec_bias, dim=-1)

        token_trans_bias = []
        for layer in self.token_trans_proj_z:
            token_trans_bias.append(layer(z))
        token_trans_bias = torch.cat(token_trans_bias, dim=-1)

        return q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias

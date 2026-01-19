"""Run Boltz model submodules step by step.

This module provides functions to execute individual submodules of the Boltz model
using data from inspect_model_and_data. This allows detailed inspection of
intermediate outputs at each stage of the forward pass.
"""

from typing import Any, Optional

import torch
from torch import Tensor
from pytorch_lightning import LightningModule


def run_submodules_step_by_step(
    model: LightningModule,
    batch: dict[str, Tensor],
    recycling_steps: Optional[int] = None,
    num_sampling_steps: Optional[int] = None,
    diffusion_samples: Optional[int] = None,
    max_parallel_samples: Optional[int] = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run model submodules step by step and return intermediate outputs.

    This function mirrors the forward pass of Boltz1/Boltz2 but stores
    intermediate outputs at each stage for inspection.

    Parameters
    ----------
    model : LightningModule
        The Boltz model (Boltz1 or Boltz2).
    batch : dict[str, Tensor]
        The input batch from inspect_model_and_data.
    recycling_steps : Optional[int]
        Number of recycling iterations. If None, uses model's predict_args.
    num_sampling_steps : Optional[int]
        Number of diffusion sampling steps. If None, uses model's predict_args.
    diffusion_samples : Optional[int]
        Number of diffusion samples (multiplicity). If None, uses model's predict_args.
    max_parallel_samples : Optional[int]
        Maximum number of parallel samples. If None, uses model's predict_args.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict[str, Any]
        Dictionary containing intermediate outputs from each submodule:
        - 'input_embedder_output': Output from input embedder (s_inputs)
        - 's_init': Initial sequence embeddings
        - 'z_init': Initial pairwise embeddings
        - 'relative_position_encoding': Relative position encodings
        - 'msa_output': Output from MSA module
        - 'pairformer_output': (s, z) from pairformer
        - 'distogram_output': Distogram predictions
        - 'diffusion_conditioning_output': Conditioning for diffusion
        - 'structure_output': Structure predictions (if run)
        - 'confidence_output': Confidence predictions (if available)
    """
    # Get default parameters from model's predict_args if not provided
    predict_args = getattr(model, 'predict_args', {}) or {}
    if recycling_steps is None:
        recycling_steps = predict_args.get('recycling_steps', 0)
    if num_sampling_steps is None:
        num_sampling_steps = predict_args.get('sampling_steps', None)
    if diffusion_samples is None:
        diffusion_samples = predict_args.get('diffusion_samples', 1)
    if max_parallel_samples is None:
        max_parallel_samples = predict_args.get('max_parallel_samples', None)

    outputs = {}
    model.eval()

    if verbose:
        print("\n" + "=" * 70)
        print("RUNNING SUBMODULES STEP BY STEP")
        print("=" * 70)

    with torch.no_grad():
        # Step 1: Input Embedder
        if verbose:
            print("\n[1/8] Running input_embedder...")
        s_inputs = model.input_embedder(batch)
        outputs['input_embedder_output'] = s_inputs
        if verbose:
            print(f"      Output shape: {s_inputs.shape}")

        # Step 2: Initialize sequence embeddings
        if verbose:
            print("\n[2/8] Running s_init (sequence initialization)...")
        s_init = model.s_init(s_inputs)
        outputs['s_init'] = s_init
        if verbose:
            print(f"      Output shape: {s_init.shape}")

        # Step 3: Initialize pairwise embeddings
        if verbose:
            print("\n[3/8] Running z_init (pairwise initialization)...")
        z_init = (
            model.z_init_1(s_inputs)[:, :, None]
            + model.z_init_2(s_inputs)[:, None, :]
        )

        # Add relative position encoding
        if verbose:
            print("      Adding relative position encoding...")
        relative_position_encoding = model.rel_pos(batch)
        outputs['relative_position_encoding'] = relative_position_encoding
        z_init = z_init + relative_position_encoding

        # Add token bonds
        if verbose:
            print("      Adding token bonds...")
        z_init = z_init + model.token_bonds(batch["token_bonds"].float())

        # Add bond type features if available
        if hasattr(model, 'bond_type_feature') and model.bond_type_feature:
            if verbose:
                print("      Adding bond type features...")
            z_init = z_init + model.token_bonds_type(batch["type_bonds"].long())

        # Add contact conditioning
        if verbose:
            print("      Adding contact conditioning...")
        z_init = z_init + model.contact_conditioning(batch)

        outputs['z_init'] = z_init
        if verbose:
            print(f"      Output shape: {z_init.shape}")

        # Step 4: Initialize recycling tensors
        s = torch.zeros_like(s_init)
        z = torch.zeros_like(z_init)

        # Compute masks
        mask = batch["token_pad_mask"].float()
        pair_mask = mask[:, :, None] * mask[:, None, :]

        # Step 5: Run recycling iterations with MSA and Pairformer
        if verbose:
            print(f"\n[4/8] Running {recycling_steps + 1} recycling iteration(s)...")

        for i in range(recycling_steps + 1):
            if verbose:
                print(f"\n      --- Recycling iteration {i + 1}/{recycling_steps + 1} ---")

            # Apply recycling
            s = s_init + model.s_recycle(model.s_norm(s))
            z = z_init + model.z_recycle(model.z_norm(z))

            # Template module (if available)
            if hasattr(model, 'use_templates') and model.use_templates:
                if verbose:
                    print("      Running template_module...")
                template_module = model.template_module
                if hasattr(template_module, '_orig_mod'):
                    template_module = template_module._orig_mod
                z = z + template_module(
                    z, batch, pair_mask, use_kernels=getattr(model, 'use_kernels', False)
                )

            # MSA module
            if verbose:
                print("      Running msa_module...")
            msa_module = model.msa_module
            if hasattr(msa_module, '_orig_mod'):
                msa_module = msa_module._orig_mod
            msa_out = msa_module(
                z, s_inputs, batch, use_kernels=getattr(model, 'use_kernels', False)
            )
            z = z + msa_out
            outputs[f'msa_output_iter_{i}'] = msa_out

            # Pairformer module
            if verbose:
                print("      Running pairformer_module...")
            pairformer_module = model.pairformer_module
            if hasattr(pairformer_module, '_orig_mod'):
                pairformer_module = pairformer_module._orig_mod
            s, z = pairformer_module(
                s,
                z,
                mask=mask,
                pair_mask=pair_mask,
                use_kernels=getattr(model, 'use_kernels', False),
            )
            outputs[f'pairformer_output_iter_{i}'] = {'s': s.clone(), 'z': z.clone()}

        outputs['pairformer_output'] = {'s': s, 'z': z}
        if verbose:
            print(f"      Final s shape: {s.shape}")
            print(f"      Final z shape: {z.shape}")

        # Step 6: Distogram module
        if verbose:
            print("\n[5/8] Running distogram_module...")
        pdistogram = model.distogram_module(z)
        outputs['distogram_output'] = pdistogram
        if verbose:
            print(f"      Output shape: {pdistogram.shape}")

        # Step 7: Diffusion conditioning and structure prediction
        run_structure = getattr(model, 'run_trunk_and_structure', True)
        skip_structure = getattr(model, 'skip_run_structure', False)

        if run_structure and not skip_structure:
            if verbose:
                print("\n[6/8] Running diffusion_conditioning...")
            q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias = (
                model.diffusion_conditioning(
                    s_trunk=s,
                    z_trunk=z,
                    relative_position_encoding=relative_position_encoding,
                    feats=batch,
                )
            )
            diffusion_conditioning = {
                'q': q,
                'c': c,
                'to_keys': to_keys,
                'atom_enc_bias': atom_enc_bias,
                'atom_dec_bias': atom_dec_bias,
                'token_trans_bias': token_trans_bias,
            }
            outputs['diffusion_conditioning_output'] = diffusion_conditioning
            if verbose:
                print(f"      q shape: {q.shape}")
                print(f"      c shape: {c.shape}")

            # Structure module (sampling)
            if verbose:
                print(f"\n[7/8] Running structure_module.sample (diffusion_samples={diffusion_samples})...")
            with torch.autocast("cuda", enabled=False):
                struct_out = model.structure_module.sample(
                    s_trunk=s.float(),
                    s_inputs=s_inputs.float(),
                    feats=batch,
                    num_sampling_steps=num_sampling_steps,
                    atom_mask=batch["atom_pad_mask"].float(),
                    multiplicity=diffusion_samples,
                    max_parallel_samples=max_parallel_samples,
                    steering_args=getattr(model, 'steering_args', None),
                    diffusion_conditioning=diffusion_conditioning,
                )
            outputs['structure_output'] = struct_out
            if verbose:
                print(f"      sample_atom_coords shape: {struct_out['sample_atom_coords'].shape}")
        else:
            if verbose:
                print("\n[6/8] Skipping diffusion_conditioning (run_trunk_and_structure=False)")
                print("\n[7/8] Skipping structure_module.sample")

        # Step 8: Confidence module (if available)
        if hasattr(model, 'confidence_prediction') and model.confidence_prediction:
            if verbose:
                print("\n[8/8] Running confidence_module...")

            if 'structure_output' in outputs:
                x_pred = outputs['structure_output']['sample_atom_coords']
            else:
                # Use ground truth coords if no structure prediction
                x_pred = batch["coords"]
                if len(x_pred.shape) == 4:
                    x_pred = x_pred.squeeze(1)

            confidence_out = model.confidence_module(
                s_inputs=s_inputs.detach(),
                s=s.detach(),
                z=z.detach(),
                x_pred=x_pred.detach(),
                feats=batch,
                pred_distogram_logits=pdistogram[:, :, :, 0].detach(),
                multiplicity=diffusion_samples,
                run_sequentially=True,
                use_kernels=getattr(model, 'use_kernels', False),
            )
            outputs['confidence_output'] = confidence_out
            if verbose:
                print(f"      plddt shape: {confidence_out.get('plddt', 'N/A')}")
        else:
            if verbose:
                print("\n[8/8] Skipping confidence_module (not available)")

    if verbose:
        print("\n" + "=" * 70)
        print("SUBMODULE EXECUTION COMPLETE")
        print("=" * 70)
        print(f"\nAvailable outputs: {list(outputs.keys())}")
        print("=" * 70 + "\n")

    return outputs


def run_single_submodule(
    model: LightningModule,
    submodule_name: str,
    batch: dict[str, Tensor],
    **kwargs,
) -> Any:
    """Run a single submodule with appropriate inputs.

    This is a convenience function for running individual submodules.
    For submodules that require outputs from previous stages, you need
    to provide them via kwargs.

    Parameters
    ----------
    model : LightningModule
        The Boltz model.
    submodule_name : str
        Name of the submodule to run. Supported:
        - 'input_embedder'
        - 'msa_module'
        - 'pairformer_module'
        - 'distogram_module'
        - 'diffusion_conditioning'
        - 'structure_module'
        - 'confidence_module'
    batch : dict[str, Tensor]
        The input batch.
    **kwargs
        Additional arguments required by the submodule.

    Returns
    -------
    Any
        Output from the submodule.
    """
    model.eval()

    with torch.no_grad():
        if submodule_name == 'input_embedder':
            return model.input_embedder(batch)

        elif submodule_name == 'msa_module':
            # Requires: z, s_inputs
            z = kwargs.get('z')
            s_inputs = kwargs.get('s_inputs')
            if z is None or s_inputs is None:
                raise ValueError("msa_module requires 'z' and 's_inputs' in kwargs")
            msa_module = model.msa_module
            if hasattr(msa_module, '_orig_mod'):
                msa_module = msa_module._orig_mod
            return msa_module(
                z, s_inputs, batch, use_kernels=getattr(model, 'use_kernels', False)
            )

        elif submodule_name == 'pairformer_module':
            # Requires: s, z
            s = kwargs.get('s')
            z = kwargs.get('z')
            if s is None or z is None:
                raise ValueError("pairformer_module requires 's' and 'z' in kwargs")
            mask = batch["token_pad_mask"].float()
            pair_mask = mask[:, :, None] * mask[:, None, :]
            pairformer_module = model.pairformer_module
            if hasattr(pairformer_module, '_orig_mod'):
                pairformer_module = pairformer_module._orig_mod
            return pairformer_module(
                s, z, mask=mask, pair_mask=pair_mask,
                use_kernels=getattr(model, 'use_kernels', False)
            )

        elif submodule_name == 'distogram_module':
            # Requires: z
            z = kwargs.get('z')
            if z is None:
                raise ValueError("distogram_module requires 'z' in kwargs")
            return model.distogram_module(z)

        elif submodule_name == 'diffusion_conditioning':
            # Requires: s, z, relative_position_encoding
            s = kwargs.get('s')
            z = kwargs.get('z')
            rel_pos = kwargs.get('relative_position_encoding')
            if s is None or z is None:
                raise ValueError(
                    "diffusion_conditioning requires 's', 'z', and "
                    "'relative_position_encoding' in kwargs"
                )
            if rel_pos is None:
                rel_pos = model.rel_pos(batch)
            return model.diffusion_conditioning(
                s_trunk=s, z_trunk=z,
                relative_position_encoding=rel_pos, feats=batch
            )

        elif submodule_name == 'structure_module':
            # Requires: s, s_inputs, diffusion_conditioning
            s = kwargs.get('s')
            s_inputs = kwargs.get('s_inputs')
            diff_cond = kwargs.get('diffusion_conditioning')
            if s is None or s_inputs is None or diff_cond is None:
                raise ValueError(
                    "structure_module requires 's', 's_inputs', and "
                    "'diffusion_conditioning' in kwargs"
                )
            num_steps = kwargs.get('num_sampling_steps', None)
            with torch.autocast("cuda", enabled=False):
                return model.structure_module.sample(
                    s_trunk=s.float(),
                    s_inputs=s_inputs.float(),
                    feats=batch,
                    num_sampling_steps=num_steps,
                    atom_mask=batch["atom_pad_mask"].float(),
                    multiplicity=kwargs.get('multiplicity', 1),
                    max_parallel_samples=None,
                    steering_args=getattr(model, 'steering_args', None),
                    diffusion_conditioning=diff_cond,
                )

        elif submodule_name == 'confidence_module':
            # Requires: s_inputs, s, z, x_pred, pred_distogram_logits
            if not hasattr(model, 'confidence_module'):
                raise ValueError("Model does not have confidence_module")
            s_inputs = kwargs.get('s_inputs')
            s = kwargs.get('s')
            z = kwargs.get('z')
            x_pred = kwargs.get('x_pred')
            pred_distogram = kwargs.get('pred_distogram_logits')
            if any(v is None for v in [s_inputs, s, z, x_pred, pred_distogram]):
                raise ValueError(
                    "confidence_module requires 's_inputs', 's', 'z', 'x_pred', "
                    "and 'pred_distogram_logits' in kwargs"
                )
            return model.confidence_module(
                s_inputs=s_inputs.detach(),
                s=s.detach(),
                z=z.detach(),
                x_pred=x_pred.detach(),
                feats=batch,
                pred_distogram_logits=pred_distogram.detach(),
                multiplicity=kwargs.get('multiplicity', 1),
                run_sequentially=True,
                use_kernels=getattr(model, 'use_kernels', False),
            )

        else:
            raise ValueError(f"Unknown submodule: {submodule_name}")


def print_submodule_outputs_summary(outputs: dict[str, Any]) -> None:
    """Print a summary of submodule outputs.

    Parameters
    ----------
    outputs : dict[str, Any]
        Dictionary of outputs from run_submodules_step_by_step.
    """
    print("\n" + "=" * 70)
    print("SUBMODULE OUTPUTS SUMMARY")
    print("=" * 70)

    for key, value in outputs.items():
        print(f"\n{key}:")
        if isinstance(value, torch.Tensor):
            print(f"  Shape: {tuple(value.shape)}")
            print(f"  Dtype: {value.dtype}")
            print(f"  Device: {value.device}")
            print(f"  Min/Max: {value.min().item():.4f} / {value.max().item():.4f}")
        elif isinstance(value, dict):
            for k, v in value.items():
                if isinstance(v, torch.Tensor):
                    print(f"  {k}: Tensor{tuple(v.shape)}")
                else:
                    print(f"  {k}: {type(v).__name__}")
        else:
            print(f"  Type: {type(value).__name__}")

    print("\n" + "=" * 70)

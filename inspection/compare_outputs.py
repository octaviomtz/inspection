"""Compare outputs from step-by-step submodule execution vs full forward pass.

This module provides functions to verify that running submodules step by step
produces the same outputs as the model's forward pass.
"""

from typing import Any, Optional

import torch
from torch import Tensor
from pytorch_lightning import LightningModule

from inspection.run_submodules import run_submodules_step_by_step


def compare_tensors(
    tensor1: Tensor,
    tensor2: Tensor,
    name: str,
    rtol: float = 1e-4,
    atol: float = 1e-5,
) -> dict[str, Any]:
    """Compare two tensors and return comparison metrics.

    Parameters
    ----------
    tensor1 : Tensor
        First tensor (reference).
    tensor2 : Tensor
        Second tensor (to compare).
    name : str
        Name for reporting.
    rtol : float
        Relative tolerance.
    atol : float
        Absolute tolerance.

    Returns
    -------
    dict[str, Any]
        Comparison results including max_diff, mean_diff, allclose status.
    """
    if tensor1.shape != tensor2.shape:
        return {
            'name': name,
            'match': False,
            'error': f"Shape mismatch: {tensor1.shape} vs {tensor2.shape}",
        }

    diff = (tensor1 - tensor2).abs()
    max_diff = diff.max().item()
    mean_diff = diff.mean().item()
    is_close = torch.allclose(tensor1, tensor2, rtol=rtol, atol=atol)

    return {
        'name': name,
        'match': is_close,
        'max_diff': max_diff,
        'mean_diff': mean_diff,
        'rtol': rtol,
        'atol': atol,
        'shape': tuple(tensor1.shape),
    }


def run_model_forward(
    model: LightningModule,
    batch: dict[str, Tensor],
    recycling_steps: int,
    num_sampling_steps: Optional[int] = None,
    diffusion_samples: int = 1,
) -> dict[str, Tensor]:
    """Run the model's forward pass directly.

    Parameters
    ----------
    model : LightningModule
        The Boltz model.
    batch : dict[str, Tensor]
        Input batch.
    recycling_steps : int
        Number of recycling steps.
    num_sampling_steps : Optional[int]
        Number of diffusion sampling steps.
    diffusion_samples : int
        Number of diffusion samples.

    Returns
    -------
    dict[str, Tensor]
        Output from model forward pass.
    """
    model.eval()
    with torch.no_grad():
        out = model(
            feats=batch,
            recycling_steps=recycling_steps,
            num_sampling_steps=num_sampling_steps,
            diffusion_samples=diffusion_samples,
            run_confidence_sequentially=True,
        )
    return out


def compare_step_by_step_vs_forward(
    model: LightningModule,
    batch: dict[str, Tensor],
    recycling_steps: Optional[int] = None,
    num_sampling_steps: Optional[int] = None,
    verbose: bool = True,
    compare_structure: bool = False,
    seed: Optional[int] = None,
) -> dict[str, Any]:
    """Compare step-by-step submodule execution against full forward pass.

    This function runs both approaches and compares their outputs to verify
    correctness. By default, only deterministic outputs (s, z, pdistogram)
    are compared since diffusion sampling is stochastic.

    Parameters
    ----------
    model : LightningModule
        The Boltz model (Boltz1 or Boltz2).
    batch : dict[str, Tensor]
        Input batch from inspect_model_and_data.
    recycling_steps : Optional[int]
        Number of recycling steps. If None, uses model's predict_args.
    num_sampling_steps : Optional[int]
        Number of diffusion sampling steps. If None, uses model's predict_args.
    verbose : bool
        Print detailed comparison results.
    compare_structure : bool
        Whether to compare structure outputs (requires same random seed).
    seed : Optional[int]
        Random seed for reproducible structure comparison.

    Returns
    -------
    dict[str, Any]
        Comparison results with keys:
        - 'all_match': bool, whether all compared tensors match
        - 'comparisons': list of individual comparison results
        - 'forward_output': output from forward pass
        - 'step_by_step_output': output from step-by-step execution
    """
    # Get default parameters from model's predict_args
    predict_args = getattr(model, 'predict_args', {}) or {}
    if recycling_steps is None:
        recycling_steps = predict_args.get('recycling_steps', 0)
    if num_sampling_steps is None:
        num_sampling_steps = predict_args.get('sampling_steps', None)

    if verbose:
        print("\n" + "=" * 70)
        print("COMPARING STEP-BY-STEP VS FORWARD PASS")
        print("=" * 70)
        print(f"Recycling steps: {recycling_steps}")
        print(f"Sampling steps: {num_sampling_steps}")

    # Set seed if provided for reproducible structure comparison
    if seed is not None and compare_structure:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    # Run forward pass
    if verbose:
        print("\n[1/3] Running full forward pass...")
    forward_out = run_model_forward(
        model=model,
        batch=batch,
        recycling_steps=recycling_steps,
        num_sampling_steps=num_sampling_steps,
        diffusion_samples=1,
    )

    # Reset seed for step-by-step if needed
    if seed is not None and compare_structure:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)

    # Run step-by-step
    if verbose:
        print("\n[2/3] Running step-by-step submodule execution...")
    step_out = run_submodules_step_by_step(
        model=model,
        batch=batch,
        recycling_steps=recycling_steps,
        num_sampling_steps=num_sampling_steps,
        verbose=False,
    )

    # Compare outputs
    if verbose:
        print("\n[3/3] Comparing outputs...")

    comparisons = []

    # Compare s (sequence embeddings)
    if 'pairformer_output' in step_out and 's' in forward_out:
        comp = compare_tensors(
            forward_out['s'],
            step_out['pairformer_output']['s'],
            name='s (sequence embeddings)',
        )
        comparisons.append(comp)

    # Compare z (pair embeddings)
    if 'pairformer_output' in step_out and 'z' in forward_out:
        comp = compare_tensors(
            forward_out['z'],
            step_out['pairformer_output']['z'],
            name='z (pair embeddings)',
        )
        comparisons.append(comp)

    # Compare pdistogram
    if 'distogram_output' in step_out and 'pdistogram' in forward_out:
        comp = compare_tensors(
            forward_out['pdistogram'],
            step_out['distogram_output'],
            name='pdistogram (distogram predictions)',
        )
        comparisons.append(comp)

    # Compare structure outputs if requested
    if compare_structure and 'structure_output' in step_out:
        if 'sample_atom_coords' in forward_out and 'sample_atom_coords' in step_out['structure_output']:
            comp = compare_tensors(
                forward_out['sample_atom_coords'],
                step_out['structure_output']['sample_atom_coords'],
                name='sample_atom_coords (structure predictions)',
                rtol=1e-3,
                atol=1e-4,
            )
            comparisons.append(comp)

    # Compare confidence outputs if available
    if 'confidence_output' in step_out:
        if 'plddt' in forward_out and 'plddt' in step_out['confidence_output']:
            comp = compare_tensors(
                forward_out['plddt'],
                step_out['confidence_output']['plddt'],
                name='plddt (confidence)',
            )
            comparisons.append(comp)

        if 'pde' in forward_out and 'pde' in step_out['confidence_output']:
            comp = compare_tensors(
                forward_out['pde'],
                step_out['confidence_output']['pde'],
                name='pde (predicted distance error)',
            )
            comparisons.append(comp)

    # Determine overall match
    all_match = all(c.get('match', False) for c in comparisons)

    # Print results
    if verbose:
        print("\n" + "-" * 70)
        print("COMPARISON RESULTS")
        print("-" * 70)

        for comp in comparisons:
            status = "PASS" if comp.get('match', False) else "FAIL"
            print(f"\n{comp['name']}:")
            if 'error' in comp:
                print(f"  Status: {status}")
                print(f"  Error: {comp['error']}")
            else:
                print(f"  Status: {status}")
                print(f"  Shape: {comp['shape']}")
                print(f"  Max diff: {comp['max_diff']:.2e}")
                print(f"  Mean diff: {comp['mean_diff']:.2e}")
                print(f"  Tolerance: rtol={comp['rtol']}, atol={comp['atol']}")

        print("\n" + "-" * 70)
        overall = "ALL TESTS PASSED" if all_match else "SOME TESTS FAILED"
        print(f"OVERALL: {overall}")
        print("-" * 70 + "\n")

    return {
        'all_match': all_match,
        'comparisons': comparisons,
        'forward_output': forward_out,
        'step_by_step_output': step_out,
    }


def verify_submodule_outputs(
    model: LightningModule,
    batch: dict[str, Tensor],
    recycling_steps: Optional[int] = None,
    num_sampling_steps: Optional[int] = None,
) -> bool:
    """Quick verification that step-by-step outputs match forward pass.

    This is a convenience function that returns True if all deterministic
    outputs match, False otherwise.

    Parameters
    ----------
    model : LightningModule
        The Boltz model.
    batch : dict[str, Tensor]
        Input batch.
    recycling_steps : Optional[int]
        Number of recycling steps. If None, uses model's predict_args.
    num_sampling_steps : Optional[int]
        Number of sampling steps. If None, uses model's predict_args.

    Returns
    -------
    bool
        True if all deterministic outputs match within tolerance.
    """
    result = compare_step_by_step_vs_forward(
        model=model,
        batch=batch,
        recycling_steps=recycling_steps,
        num_sampling_steps=num_sampling_steps,
        verbose=False,
        compare_structure=False,
    )
    return result['all_match']

#!/usr/bin/env python
"""Example script demonstrating how to use the inspection utilities.

This script shows how to:
1. Access the captured batch and model after running boltz predict
2. Inspect model submodules
3. Register hooks to capture intermediate outputs

Usage:
    # First run boltz predict as normal:
    boltz predict examples/cycle_peptide_fake.yaml

    # Then in Python, you can access the captured data:
    from inspection import get_batch, get_model, print_batch_summary
"""

import sys
sys.path.insert(0, "/mnt/c/Users/octav/Documents/claude/proteinEBM/boltz")

from inspection import (
    get_batch,
    get_model,
    get_data_module,
    print_batch_summary,
    print_model_structure,
    list_submodules,
    get_submodule,
    print_module_info,
    register_forward_hooks,
    remove_hooks,
)


def example_after_prediction():
    """Example: Access data after boltz predict has run."""
    # Get the captured data (only works if boltz predict was run in same process)
    batch = get_batch()
    model = get_model()

    if batch is None:
        print("No batch captured. Run boltz predict first.")
        return

    print("\n=== Batch Summary ===")
    print_batch_summary(batch)

    print("\n=== Model Structure ===")
    print_model_structure(model, max_depth=2)

    # List top-level submodules
    print("\n=== Top-level Submodules ===")
    for name, module in list_submodules(model, max_depth=1):
        print(f"  {name}: {module.__class__.__name__}")

    # Get specific submodule
    print("\n=== Pairformer Module Info ===")
    try:
        pairformer = get_submodule(model, "pairformer_module")
        print_module_info(pairformer, "pairformer_module")
    except AttributeError:
        print("Pairformer module not found")


def example_standalone_data_loading():
    """Example: Load data standalone without running full prediction."""
    from pathlib import Path
    from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
    from boltz.data.types import Manifest

    # This is an example - adjust paths as needed
    # You would need to have processed data from a previous run
    print("This example requires processed data from a previous boltz predict run.")
    print("Paths would be: out_dir/processed/*, out_dir/msa/*, etc.")


def example_forward_hooks():
    """Example: Register hooks to capture intermediate outputs during forward pass."""
    import torch

    model = get_model()
    batch = get_batch()

    if model is None or batch is None:
        print("No model/batch captured. Run boltz predict first.")
        return

    # Register hooks on specific modules
    modules_to_hook = ["input_embedder", "pairformer_module", "structure_module"]
    storage, handles = register_forward_hooks(model, modules_to_hook)

    print(f"Registered {len(handles)} hooks")
    print("Run a forward pass to capture outputs...")

    # Note: To actually capture outputs, you'd need to run a forward pass:
    # with torch.no_grad():
    #     model(batch)  # This would populate 'storage'

    # Clean up hooks
    remove_hooks(handles)
    print("Hooks removed")


def list_batch_keys():
    """List all keys in the captured batch."""
    batch = get_batch()
    if batch is None:
        print("No batch captured.")
        return

    print("\n=== Batch Keys ===")
    for key in sorted(batch.keys()):
        value = batch[key]
        if hasattr(value, 'shape'):
            print(f"  {key}: Tensor{tuple(value.shape)}")
        elif isinstance(value, list):
            print(f"  {key}: list[{len(value)}]")
        else:
            print(f"  {key}: {type(value).__name__}")


if __name__ == "__main__":
    print("Inspection Example Script")
    print("=" * 60)
    print("\nThis script demonstrates the inspection utilities.")
    print("Run 'boltz predict' first, then import these functions.")
    print("\nAvailable functions:")
    print("  - example_after_prediction()")
    print("  - example_forward_hooks()")
    print("  - list_batch_keys()")

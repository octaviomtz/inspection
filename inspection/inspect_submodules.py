"""Interactive inspection of Boltz model submodules.

This module provides functions to load input data and inspect model submodules
during debugging. Set a breakpoint after calling `inspect_model_and_data` to
interactively explore the model and data.
"""

from typing import Optional, Any
import torch
from torch import nn
from pytorch_lightning import LightningModule, LightningDataModule


def inspect_model_and_data(
    model: LightningModule,
    data_module: LightningDataModule,
    device: Optional[torch.device] = None,
) -> dict[str, Any]:
    """Load input data and prepare model for inspection.

    Call this function and set a breakpoint on the return statement to
    interactively inspect the model submodules and input data.

    Parameters
    ----------
    model : LightningModule
        The Boltz model (Boltz1 or Boltz2).
    data_module : LightningDataModule
        The data module containing input data.
    device : Optional[torch.device]
        Device to load data to. Defaults to CPU.

    Returns
    -------
    dict[str, Any]
        Dictionary containing:
        - 'model': The model instance
        - 'batch': The first batch of input data
        - 'submodules': Dict of top-level submodule names to modules
        - 'batch_keys': List of keys in the batch

    Example
    -------
    Set a breakpoint after this call to inspect:

        >>> result = inspect_model_and_data(model, data_module)
        >>> # Breakpoint here - explore result['batch'], result['submodules'], etc.
        >>> batch = result['batch']
        >>> pairformer = result['submodules']['pairformer_module']
    """
    # Load the first batch
    dataloader = data_module.predict_dataloader()
    batch = next(iter(dataloader))

    # Transfer to device if specified
    if device is not None:
        batch = data_module.transfer_batch_to_device(batch, device, 0)

    # Collect top-level submodules
    submodules = {name: module for name, module in model.named_children()}

    # Get batch info
    batch_keys = list(batch.keys())
    batch_info = {}
    for key, value in batch.items():
        if isinstance(value, torch.Tensor):
            batch_info[key] = {
                'shape': tuple(value.shape),
                'dtype': str(value.dtype),
                'device': str(value.device),
            }
        elif isinstance(value, list):
            batch_info[key] = {'type': 'list', 'length': len(value)}
        else:
            batch_info[key] = {'type': type(value).__name__}

    # Model info
    model_info = {
        'class': model.__class__.__name__,
        'submodule_names': list(submodules.keys()),
        'total_params': sum(p.numel() for p in model.parameters()),
    }

    result = {
        'model': model,
        'batch': batch,
        'submodules': submodules,
        'batch_keys': batch_keys,
        'batch_info': batch_info,
        'model_info': model_info,
    }

    # Print summary for convenience
    print("\n" + "=" * 70)
    print("INSPECTION DATA READY")
    print("=" * 70)
    print(f"\nModel: {model_info['class']}")
    print(f"Total parameters: {model_info['total_params']:,}")
    print(f"\nSubmodules ({len(submodules)}):")
    for name in submodules.keys():
        print(f"  - {name}")
    print(f"\nBatch keys ({len(batch_keys)}):")
    for key in sorted(batch_keys)[:10]:
        info = batch_info[key]
        if 'shape' in info:
            print(f"  - {key}: Tensor{info['shape']}")
        else:
            print(f"  - {key}: {info.get('type', 'unknown')}")
    if len(batch_keys) > 10:
        print(f"  ... and {len(batch_keys) - 10} more")
    print("\n" + "=" * 70)
    print("Set a breakpoint here to inspect 'result' variable")
    print("=" * 70 + "\n")

    # === SET BREAKPOINT HERE ===
    # In VSCode, click in the gutter on the left of the 'return' line below
    # to set a breakpoint. Then you can inspect:
    #   - result['batch'] - the input data
    #   - result['submodules'] - dict of submodules
    #   - result['model'] - the full model
    return result  # <-- SET BREAKPOINT HERE


def get_submodule_by_path(model: nn.Module, path: str) -> nn.Module:
    """Get a nested submodule by dot-separated path.

    Parameters
    ----------
    model : nn.Module
        The root model.
    path : str
        Dot-separated path, e.g., "pairformer_module.blocks.0.attention".

    Returns
    -------
    nn.Module
        The requested submodule.
    """
    parts = path.split(".")
    current = model
    for part in parts:
        if part.isdigit():
            current = current[int(part)]
        else:
            current = getattr(current, part)
    return current


def list_all_submodules(model: nn.Module, max_depth: int = 3) -> list[str]:
    """List all submodule paths up to a certain depth.

    Parameters
    ----------
    model : nn.Module
        The model to inspect.
    max_depth : int
        Maximum depth to traverse.

    Returns
    -------
    list[str]
        List of dot-separated paths to submodules.
    """
    paths = []

    def _traverse(module: nn.Module, prefix: str, depth: int):
        if depth > max_depth:
            return
        for name, child in module.named_children():
            full_path = f"{prefix}.{name}" if prefix else name
            paths.append(full_path)
            _traverse(child, full_path, depth + 1)

    _traverse(model, "", 1)
    return paths


def run_submodule_forward(
    model: nn.Module,
    submodule_path: str,
    *args,
    **kwargs,
) -> Any:
    """Run forward pass on a specific submodule.

    Parameters
    ----------
    model : nn.Module
        The root model.
    submodule_path : str
        Path to the submodule.
    *args, **kwargs
        Arguments to pass to the submodule's forward method.

    Returns
    -------
    Any
        Output from the submodule.
    """
    submodule = get_submodule_by_path(model, submodule_path)
    with torch.no_grad():
        return submodule(*args, **kwargs)

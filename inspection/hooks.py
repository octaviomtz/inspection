"""Hooks for inspecting Boltz data and models during prediction.

This module provides a minimal interface to capture input data and model
before the Trainer runs, allowing inspection without modifying core code.
"""

from typing import Optional
import torch
from pytorch_lightning import LightningModule, LightningDataModule

from inspection.data_loader import load_batch_from_datamodule, print_batch_summary
from inspection.model_inspector import print_model_structure


# Global storage for inspection data
_inspection_data = {
    "batch": None,
    "model": None,
    "data_module": None,
}


def pre_predict_hook(
    model: LightningModule,
    data_module: LightningDataModule,
    device: Optional[torch.device] = None,
    load_batch: bool = True,
    batch_idx: int = 0,
    print_summary: bool = False,
) -> dict:
    """Hook called before trainer.predict() to capture data and model.

    Parameters
    ----------
    model : LightningModule
        The model being used for prediction.
    data_module : LightningDataModule
        The data module providing input data.
    device : Optional[torch.device]
        Device to load the batch to. If None, uses CPU.
    load_batch : bool
        Whether to load and store the first batch.
    batch_idx : int
        Which batch to load (default: 0).
    print_summary : bool
        Whether to print a summary of the data.

    Returns
    -------
    dict
        Dictionary containing 'batch', 'model', and 'data_module'.
    """
    _inspection_data["model"] = model
    _inspection_data["data_module"] = data_module

    if load_batch:
        batch = load_batch_from_datamodule(data_module, batch_idx=batch_idx, device=device)
        _inspection_data["batch"] = batch

        if print_summary:
            print("\n" + "=" * 60)
            print("INSPECTION HOOK: Data captured before prediction")
            print("=" * 60)
            print_batch_summary(batch)
            print_model_structure(model, max_depth=1)

    return _inspection_data


def get_inspection_data() -> dict:
    """Get the captured inspection data.

    Returns
    -------
    dict
        Dictionary containing 'batch', 'model', and 'data_module'.
    """
    return _inspection_data


def get_batch() -> Optional[dict]:
    """Get the captured batch.

    Returns
    -------
    Optional[dict]
        The captured batch, or None if not captured.
    """
    return _inspection_data.get("batch")


def get_model() -> Optional[LightningModule]:
    """Get the captured model.

    Returns
    -------
    Optional[LightningModule]
        The captured model, or None if not captured.
    """
    return _inspection_data.get("model")


def get_data_module() -> Optional[LightningDataModule]:
    """Get the captured data module.

    Returns
    -------
    Optional[LightningDataModule]
        The captured data module, or None if not captured.
    """
    return _inspection_data.get("data_module")


def clear_inspection_data() -> None:
    """Clear all captured inspection data."""
    _inspection_data["batch"] = None
    _inspection_data["model"] = None
    _inspection_data["data_module"] = None

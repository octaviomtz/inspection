"""Utilities for loading input data from Boltz DataModules."""

from typing import Optional, Iterator
import torch
from torch.utils.data import DataLoader
from pytorch_lightning import LightningDataModule


def get_dataloader(data_module: LightningDataModule) -> DataLoader:
    """Get the prediction dataloader from a DataModule.

    Parameters
    ----------
    data_module : LightningDataModule
        The data module (e.g., Boltz2InferenceDataModule).

    Returns
    -------
    DataLoader
        The prediction dataloader.
    """
    return data_module.predict_dataloader()


def load_batch_from_datamodule(
    data_module: LightningDataModule,
    batch_idx: int = 0,
    device: Optional[torch.device] = None,
) -> dict:
    """Load a specific batch from the DataModule.

    Parameters
    ----------
    data_module : LightningDataModule
        The data module (e.g., Boltz2InferenceDataModule).
    batch_idx : int
        The index of the batch to load (default: 0 for first batch).
    device : Optional[torch.device]
        Device to transfer the batch to. If None, keeps on CPU.

    Returns
    -------
    dict
        The batch dictionary containing all features.
    """
    dataloader = get_dataloader(data_module)

    for idx, batch in enumerate(dataloader):
        if idx == batch_idx:
            if device is not None:
                batch = data_module.transfer_batch_to_device(batch, device, 0)
            return batch

    raise IndexError(f"Batch index {batch_idx} out of range. DataLoader has fewer batches.")


def iterate_batches(
    data_module: LightningDataModule,
    device: Optional[torch.device] = None,
    max_batches: Optional[int] = None,
) -> Iterator[dict]:
    """Iterate over batches from the DataModule.

    Parameters
    ----------
    data_module : LightningDataModule
        The data module (e.g., Boltz2InferenceDataModule).
    device : Optional[torch.device]
        Device to transfer batches to. If None, keeps on CPU.
    max_batches : Optional[int]
        Maximum number of batches to yield. If None, yields all.

    Yields
    ------
    dict
        Batch dictionaries containing features.
    """
    dataloader = get_dataloader(data_module)

    for idx, batch in enumerate(dataloader):
        if max_batches is not None and idx >= max_batches:
            break
        if device is not None:
            batch = data_module.transfer_batch_to_device(batch, device, 0)
        yield batch


def describe_batch(batch: dict) -> dict:
    """Get a description of batch contents.

    Parameters
    ----------
    batch : dict
        The batch dictionary.

    Returns
    -------
    dict
        Dictionary with key -> (type, shape/len) mappings.
    """
    description = {}
    for key, value in batch.items():
        if isinstance(value, torch.Tensor):
            description[key] = ("Tensor", tuple(value.shape), str(value.dtype))
        elif isinstance(value, list):
            description[key] = ("list", len(value))
        else:
            description[key] = (type(value).__name__,)
    return description


def print_batch_summary(batch: dict) -> None:
    """Print a summary of batch contents.

    Parameters
    ----------
    batch : dict
        The batch dictionary.
    """
    print("=" * 60)
    print("BATCH SUMMARY")
    print("=" * 60)

    tensors = []
    others = []

    for key, value in sorted(batch.items()):
        if isinstance(value, torch.Tensor):
            tensors.append((key, value))
        else:
            others.append((key, value))

    print("\nTensor features:")
    print("-" * 60)
    for key, tensor in tensors:
        print(f"  {key:40} shape={tuple(tensor.shape)}, dtype={tensor.dtype}")

    print("\nOther features:")
    print("-" * 60)
    for key, value in others:
        if isinstance(value, list):
            print(f"  {key:40} list[{len(value)}]")
        else:
            print(f"  {key:40} {type(value).__name__}")

    print("=" * 60)

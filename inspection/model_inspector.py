"""Utilities for inspecting Boltz model submodules."""

from typing import Optional, Iterator
import torch
from torch import nn


def list_submodules(
    model: nn.Module,
    max_depth: Optional[int] = None,
    prefix: str = "",
) -> list[tuple[str, nn.Module]]:
    """List all submodules of a model.

    Parameters
    ----------
    model : nn.Module
        The model to inspect.
    max_depth : Optional[int]
        Maximum depth to traverse. None means unlimited.
    prefix : str
        Prefix for module names (used internally for recursion).

    Returns
    -------
    list[tuple[str, nn.Module]]
        List of (name, module) tuples.
    """
    result = []

    for name, module in model.named_children():
        full_name = f"{prefix}.{name}" if prefix else name
        result.append((full_name, module))

        current_depth = full_name.count(".") + 1
        if max_depth is None or current_depth < max_depth:
            result.extend(list_submodules(module, max_depth, full_name))

    return result


def get_submodule(model: nn.Module, path: str) -> nn.Module:
    """Get a submodule by its path.

    Parameters
    ----------
    model : nn.Module
        The model to search in.
    path : str
        Dot-separated path to the submodule (e.g., "pairformer_module.blocks.0").

    Returns
    -------
    nn.Module
        The requested submodule.

    Raises
    ------
    AttributeError
        If the path is invalid.
    """
    parts = path.split(".")
    current = model
    for part in parts:
        if part.isdigit():
            current = current[int(part)]
        else:
            current = getattr(current, part)
    return current


def print_model_structure(
    model: nn.Module,
    max_depth: int = 2,
    show_params: bool = True,
) -> None:
    """Print the model structure in a tree format.

    Parameters
    ----------
    model : nn.Module
        The model to inspect.
    max_depth : int
        Maximum depth to show (default: 2).
    show_params : bool
        Whether to show parameter counts.
    """
    print("=" * 70)
    print(f"MODEL STRUCTURE: {model.__class__.__name__}")
    print("=" * 70)

    def count_params(m: nn.Module) -> int:
        return sum(p.numel() for p in m.parameters())

    def print_tree(module: nn.Module, prefix: str, depth: int) -> None:
        children = list(module.named_children())
        for i, (name, child) in enumerate(children):
            is_last = i == len(children) - 1
            connector = "\\-- " if is_last else "|-- "
            child_prefix = "    " if is_last else "|   "

            param_str = ""
            if show_params:
                params = count_params(child)
                if params > 0:
                    if params >= 1e9:
                        param_str = f" ({params/1e9:.2f}B params)"
                    elif params >= 1e6:
                        param_str = f" ({params/1e6:.2f}M params)"
                    elif params >= 1e3:
                        param_str = f" ({params/1e3:.2f}K params)"
                    else:
                        param_str = f" ({params} params)"

            print(f"{prefix}{connector}{name}: {child.__class__.__name__}{param_str}")

            if depth < max_depth:
                print_tree(child, prefix + child_prefix, depth + 1)

    total_params = count_params(model)
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)

    print(f"\nTotal parameters: {total_params:,}")
    print(f"Trainable parameters: {trainable_params:,}")
    print()

    print_tree(model, "", 1)
    print("=" * 70)


def get_module_info(module: nn.Module) -> dict:
    """Get detailed information about a module.

    Parameters
    ----------
    module : nn.Module
        The module to inspect.

    Returns
    -------
    dict
        Dictionary with module information.
    """
    info = {
        "class": module.__class__.__name__,
        "total_params": sum(p.numel() for p in module.parameters()),
        "trainable_params": sum(p.numel() for p in module.parameters() if p.requires_grad),
        "buffers": list(module.named_buffers(recurse=False)),
        "direct_children": [(n, c.__class__.__name__) for n, c in module.named_children()],
    }

    # Get parameter shapes
    info["parameters"] = {
        name: tuple(param.shape)
        for name, param in module.named_parameters(recurse=False)
    }

    return info


def print_module_info(module: nn.Module, name: str = "") -> None:
    """Print detailed information about a module.

    Parameters
    ----------
    module : nn.Module
        The module to inspect.
    name : str
        Name to display (optional).
    """
    info = get_module_info(module)

    print("=" * 60)
    if name:
        print(f"MODULE: {name}")
    print(f"CLASS: {info['class']}")
    print("=" * 60)

    print(f"\nTotal parameters: {info['total_params']:,}")
    print(f"Trainable parameters: {info['trainable_params']:,}")

    if info["parameters"]:
        print("\nDirect parameters:")
        for pname, shape in info["parameters"].items():
            print(f"  {pname}: {shape}")

    if info["direct_children"]:
        print("\nDirect children:")
        for cname, cclass in info["direct_children"]:
            print(f"  {cname}: {cclass}")

    print("=" * 60)


def register_forward_hooks(
    model: nn.Module,
    module_names: Optional[list[str]] = None,
    storage: Optional[dict] = None,
) -> tuple[dict, list]:
    """Register forward hooks to capture intermediate outputs.

    Parameters
    ----------
    model : nn.Module
        The model to hook.
    module_names : Optional[list[str]]
        List of module paths to hook. If None, hooks top-level modules.
    storage : Optional[dict]
        Dictionary to store outputs in. Created if None.

    Returns
    -------
    tuple[dict, list]
        (storage dict, list of hook handles for removal)
    """
    if storage is None:
        storage = {}

    handles = []

    if module_names is None:
        module_names = [name for name, _ in model.named_children()]

    for name in module_names:
        try:
            module = get_submodule(model, name)
        except AttributeError:
            print(f"Warning: Module '{name}' not found, skipping.")
            continue

        def make_hook(n):
            def hook(mod, inp, out):
                storage[n] = {
                    "input": inp,
                    "output": out,
                }
            return hook

        handle = module.register_forward_hook(make_hook(name))
        handles.append(handle)

    return storage, handles


def remove_hooks(handles: list) -> None:
    """Remove registered hooks.

    Parameters
    ----------
    handles : list
        List of hook handles to remove.
    """
    for handle in handles:
        handle.remove()

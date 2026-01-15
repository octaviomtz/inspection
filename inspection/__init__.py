"""Inspection utilities for Boltz models."""

from inspection.data_loader import (
    load_batch_from_datamodule,
    get_dataloader,
    iterate_batches,
    describe_batch,
    print_batch_summary,
)
from inspection.model_inspector import (
    list_submodules,
    get_submodule,
    print_model_structure,
    get_module_info,
    print_module_info,
    register_forward_hooks,
    remove_hooks,
)
from inspection.hooks import (
    pre_predict_hook,
    get_inspection_data,
    get_batch,
    get_model,
    get_data_module,
    clear_inspection_data,
)
from inspection.inspect_submodules import (
    inspect_model_and_data,
    get_submodule_by_path,
    list_all_submodules,
    run_submodule_forward,
)
from inspection.run_submodules import (
    run_submodules_step_by_step,
    run_single_submodule,
    print_submodule_outputs_summary,
)
from inspection.write_submodule_outputs import (
    write_submodule_outputs,
)
from inspection.compare_outputs import (
    compare_step_by_step_vs_forward,
    verify_submodule_outputs,
    compare_tensors,
)

__all__ = [
    # Data loading
    "load_batch_from_datamodule",
    "get_dataloader",
    "iterate_batches",
    "describe_batch",
    "print_batch_summary",
    # Model inspection
    "list_submodules",
    "get_submodule",
    "print_model_structure",
    "get_module_info",
    "print_module_info",
    "register_forward_hooks",
    "remove_hooks",
    # Hooks
    "pre_predict_hook",
    "get_inspection_data",
    "get_batch",
    "get_model",
    "get_data_module",
    "clear_inspection_data",
    # Submodule inspection
    "inspect_model_and_data",
    "get_submodule_by_path",
    "list_all_submodules",
    "run_submodule_forward",
    # Step-by-step submodule execution
    "run_submodules_step_by_step",
    "run_single_submodule",
    "print_submodule_outputs_summary",
    # Write submodule outputs
    "write_submodule_outputs",
    # Output comparison
    "compare_step_by_step_vs_forward",
    "verify_submodule_outputs",
    "compare_tensors",
]

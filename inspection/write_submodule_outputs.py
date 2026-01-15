"""Write submodule outputs to disk in the same format as Trainer predictions.

This module provides functions to save the outputs from run_submodules_step_by_step
to disk, using the same format as BoltzWriter. This allows inspection of intermediate
results while maintaining compatibility with the standard output format.
"""

import json
from dataclasses import replace
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np
import torch
from torch import Tensor

from boltz.data.types import Coords, Interface, Record, Structure, StructureV2
from boltz.data.write.mmcif import to_mmcif
from boltz.data.write.pdb import to_pdb


def write_submodule_outputs(
    submodule_outputs: dict[str, Any],
    batch: dict[str, Tensor],
    output_dir: str | Path,
    data_dir: str | Path,
    output_format: Literal["pdb", "mmcif"] = "mmcif",
    boltz2: bool = False,
    write_embeddings: bool = True,
    write_confidence: bool = True,
    write_pae: bool = True,
    write_pde: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write submodule outputs to disk in the same format as Trainer predictions.

    This function takes the outputs from run_submodules_step_by_step and writes
    them to a "from_submodules" folder within the output directory, using the
    same format as BoltzWriter.

    Parameters
    ----------
    submodule_outputs : dict[str, Any]
        Dictionary of outputs from run_submodules_step_by_step. Expected keys:
        - 'structure_output': dict with 'sample_atom_coords'
        - 'pairformer_output': dict with 's' and 'z'
        - 'confidence_output': dict with plddt, pae, ptm, iptm, etc. (optional)
        - 'distogram_output': distogram predictions
    batch : dict[str, Tensor]
        The input batch containing 'record', 'atom_pad_mask', etc.
    output_dir : str | Path
        Base output directory. A "from_submodules" folder will be created inside.
    data_dir : str | Path
        Directory containing the structure .npz files (for loading structure info).
    output_format : Literal["pdb", "mmcif"]
        Output format for structure files (default: "mmcif").
    boltz2 : bool
        Whether using Boltz2 format (default: False).
    write_embeddings : bool
        Whether to write s and z embeddings (default: True).
    write_confidence : bool
        Whether to write confidence metrics (default: True).
    write_pae : bool
        Whether to write full PAE matrix (default: True).
    write_pde : bool
        Whether to write full PDE matrix (default: False).
    verbose : bool
        Print progress messages (default: True).

    Returns
    -------
    dict[str, Any]
        Dictionary containing information about written files:
        - 'output_dir': Path to the from_submodules directory
        - 'structure_files': List of written structure file paths
        - 'confidence_files': List of written confidence file paths
        - 'embedding_files': List of written embedding file paths
    """
    output_dir = Path(output_dir)
    data_dir = Path(data_dir)

    # Create the from_submodules folder
    submodules_dir = output_dir / "from_submodules"
    submodules_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("\n" + "=" * 70)
        print("WRITING SUBMODULE OUTPUTS")
        print("=" * 70)
        print(f"Output directory: {submodules_dir}")

    written_files = {
        'output_dir': submodules_dir,
        'structure_files': [],
        'confidence_files': [],
        'embedding_files': [],
    }

    # Get the records from the batch
    records: list[Record] = batch["record"]

    # Get structure coordinates
    if 'structure_output' not in submodule_outputs:
        if verbose:
            print("WARNING: No structure_output found, skipping structure files")
        return written_files

    # Get coordinates: shape [1, num_atoms, 3] from sample_atom_coords
    coords = submodule_outputs['structure_output']['sample_atom_coords']
    # Add batch dimension if needed: [batch, samples, num_atoms, 3]
    if len(coords.shape) == 3:
        coords = coords.unsqueeze(0)  # [1, num_atoms, 3] -> [1, 1, num_atoms, 3]

    pad_masks = batch["atom_pad_mask"]

    # Get confidence outputs if available
    confidence_output = submodule_outputs.get('confidence_output', {})
    has_confidence = bool(confidence_output) and write_confidence

    # Get embeddings from pairformer output
    pairformer_output = submodule_outputs.get('pairformer_output', {})
    s_embedding = pairformer_output.get('s')
    z_embedding = pairformer_output.get('z')

    # Prepare ranking (single sample, so rank is 0)
    num_samples = coords.shape[1] if len(coords.shape) == 4 else 1
    idx_to_rank = {i: i for i in range(num_samples)}

    # If we have confidence score, use it for ranking
    if has_confidence and 'plddt' in confidence_output:
        # Compute confidence score similar to predict_step
        complex_plddt = confidence_output.get('complex_plddt', torch.zeros(num_samples))
        iptm = confidence_output.get('iptm', torch.zeros(num_samples))
        ptm = confidence_output.get('ptm', torch.zeros(num_samples))

        # Use iptm if non-zero, else ptm
        if torch.is_tensor(iptm) and not torch.allclose(iptm, torch.zeros_like(iptm)):
            tm_score = iptm
        else:
            tm_score = ptm

        confidence_score = (4 * complex_plddt + tm_score) / 5

        if confidence_score.numel() > 1:
            argsort = torch.argsort(confidence_score, descending=True)
            idx_to_rank = {idx.item(): rank for rank, idx in enumerate(argsort)}

    # Iterate over records
    for record_idx, record in enumerate(records):
        if verbose:
            print(f"\nProcessing record: {record.id}")

        # Load the structure
        structure_path = data_dir / f"{record.id}.npz"
        if not structure_path.exists():
            if verbose:
                print(f"  WARNING: Structure file not found: {structure_path}")
            continue

        if boltz2:
            structure: StructureV2 = StructureV2.load(structure_path)
        else:
            structure: Structure = Structure.load(structure_path)

        # Compute chain map with masked removed
        chain_map = {}
        for i, mask in enumerate(structure.mask):
            if mask:
                chain_map[len(chain_map)] = i

        # Remove masked chains
        structure = structure.remove_invalid_chains()

        # Create structure directory
        struct_dir = submodules_dir / record.id
        struct_dir.mkdir(exist_ok=True)

        # Get coordinates for this record
        coord = coords[record_idx] if coords.shape[0] > 1 else coords[0]
        pad_mask = pad_masks[record_idx] if pad_masks.shape[0] > 1 else pad_masks[0]

        for model_idx in range(num_samples):
            # Get model coord
            if len(coord.shape) == 3:
                model_coord = coord[model_idx]
            else:
                model_coord = coord

            # Unpad
            coord_unpad = model_coord[pad_mask.bool()]
            coord_unpad = coord_unpad.cpu().numpy()

            # Get plddts
            plddts = None
            if has_confidence and 'plddt' in confidence_output:
                plddt_tensor = confidence_output['plddt']
                if len(plddt_tensor.shape) > 1 and plddt_tensor.shape[0] > model_idx:
                    plddts = plddt_tensor[model_idx]
                elif len(plddt_tensor.shape) == 1:
                    plddts = plddt_tensor

            # New atom table
            atoms = structure.atoms.copy()
            atoms["coords"] = coord_unpad
            atoms["is_present"] = True
            if boltz2:
                coord_unpad_typed = [(x,) for x in coord_unpad]
                coord_unpad_typed = np.array(coord_unpad_typed, dtype=Coords)

            # New residue table
            residues = structure.residues.copy()
            residues["is_present"] = True

            # Update the structure
            interfaces = np.array([], dtype=Interface)
            if boltz2:
                new_structure: StructureV2 = replace(
                    structure,
                    atoms=atoms,
                    residues=residues,
                    interfaces=interfaces,
                    coords=coord_unpad_typed,
                )
            else:
                new_structure: Structure = replace(
                    structure,
                    atoms=atoms,
                    residues=residues,
                    interfaces=interfaces,
                )

            # Update chain info
            chain_info = []
            for chain in new_structure.chains:
                old_chain_idx = chain_map.get(chain["asym_id"], chain["asym_id"])
                if old_chain_idx < len(record.chains):
                    old_chain_info = record.chains[old_chain_idx]
                    new_chain_info = replace(
                        old_chain_info,
                        chain_id=int(chain["asym_id"]),
                        valid=True,
                    )
                    chain_info.append(new_chain_info)

            # Create path name
            rank = idx_to_rank.get(model_idx, model_idx)
            outname = f"{record.id}_model_{rank}"

            # Save the structure
            if output_format == "pdb":
                path = struct_dir / f"{outname}.pdb"
                with path.open("w") as f:
                    f.write(to_pdb(new_structure, plddts=plddts, boltz2=boltz2))
            elif output_format == "mmcif":
                path = struct_dir / f"{outname}.cif"
                with path.open("w") as f:
                    f.write(to_mmcif(new_structure, plddts=plddts, boltz2=boltz2))
            else:
                path = struct_dir / f"{outname}.npz"
                from dataclasses import asdict
                np.savez_compressed(path, **asdict(new_structure))

            written_files['structure_files'].append(str(path))
            if verbose:
                print(f"  Wrote structure: {path.name}")

            # Save confidence summary
            if has_confidence and 'plddt' in confidence_output:
                conf_path = struct_dir / f"confidence_{record.id}_model_{rank}.json"
                confidence_summary_dict = {}

                # Compute confidence score
                complex_plddt_val = confidence_output.get('complex_plddt')
                iptm_val = confidence_output.get('iptm')
                ptm_val = confidence_output.get('ptm')

                if complex_plddt_val is not None and ptm_val is not None:
                    if iptm_val is not None and not torch.allclose(iptm_val, torch.zeros_like(iptm_val)):
                        tm_val = iptm_val
                    else:
                        tm_val = ptm_val

                    if complex_plddt_val.numel() > model_idx:
                        c_plddt = complex_plddt_val[model_idx].item() if complex_plddt_val.dim() > 0 else complex_plddt_val.item()
                    else:
                        c_plddt = complex_plddt_val.item()

                    if tm_val.numel() > model_idx:
                        t_val = tm_val[model_idx].item() if tm_val.dim() > 0 else tm_val.item()
                    else:
                        t_val = tm_val.item()

                    confidence_summary_dict["confidence_score"] = (4 * c_plddt + t_val) / 5

                for key in [
                    "ptm",
                    "iptm",
                    "ligand_iptm",
                    "protein_iptm",
                    "complex_plddt",
                    "complex_iplddt",
                    "complex_pde",
                    "complex_ipde",
                ]:
                    if key in confidence_output:
                        val = confidence_output[key]
                        if torch.is_tensor(val):
                            if val.numel() > model_idx:
                                confidence_summary_dict[key] = val[model_idx].item() if val.dim() > 0 else val.item()
                            else:
                                confidence_summary_dict[key] = val.item()

                # Handle pair_chains_iptm
                if "pair_chains_iptm" in confidence_output:
                    pair_iptm = confidence_output["pair_chains_iptm"]
                    confidence_summary_dict["chains_ptm"] = {
                        str(idx): pair_iptm[idx][idx][model_idx].item() if pair_iptm[idx][idx].numel() > model_idx else pair_iptm[idx][idx].item()
                        for idx in pair_iptm
                    }
                    confidence_summary_dict["pair_chains_iptm"] = {
                        str(idx1): {
                            str(idx2): pair_iptm[idx1][idx2][model_idx].item() if pair_iptm[idx1][idx2].numel() > model_idx else pair_iptm[idx1][idx2].item()
                            for idx2 in pair_iptm[idx1]
                        }
                        for idx1 in pair_iptm
                    }

                with conf_path.open("w") as f:
                    f.write(json.dumps(confidence_summary_dict, indent=4))

                written_files['confidence_files'].append(str(conf_path))
                if verbose:
                    print(f"  Wrote confidence: {conf_path.name}")

                # Save plddt
                if plddts is not None:
                    plddt_path = struct_dir / f"plddt_{record.id}_model_{rank}.npz"
                    plddt_np = plddts.cpu().numpy() if torch.is_tensor(plddts) else plddts
                    np.savez_compressed(plddt_path, plddt=plddt_np)
                    written_files['confidence_files'].append(str(plddt_path))
                    if verbose:
                        print(f"  Wrote plddt: {plddt_path.name}")

            # Save pae
            if write_pae and 'pae' in confidence_output:
                pae = confidence_output['pae']
                if pae.numel() > 0:
                    if pae.dim() > 2 and pae.shape[0] > model_idx:
                        pae_model = pae[model_idx]
                    else:
                        pae_model = pae
                    pae_path = struct_dir / f"pae_{record.id}_model_{rank}.npz"
                    np.savez_compressed(pae_path, pae=pae_model.cpu().numpy())
                    written_files['confidence_files'].append(str(pae_path))
                    if verbose:
                        print(f"  Wrote pae: {pae_path.name}")

            # Save pde
            if write_pde and 'pde' in confidence_output:
                pde = confidence_output['pde']
                if pde.numel() > 0:
                    if pde.dim() > 2 and pde.shape[0] > model_idx:
                        pde_model = pde[model_idx]
                    else:
                        pde_model = pde
                    pde_path = struct_dir / f"pde_{record.id}_model_{rank}.npz"
                    np.savez_compressed(pde_path, pde=pde_model.cpu().numpy())
                    written_files['confidence_files'].append(str(pde_path))
                    if verbose:
                        print(f"  Wrote pde: {pde_path.name}")

        # Save embeddings (once per record, not per model)
        if write_embeddings and s_embedding is not None and z_embedding is not None:
            s_np = s_embedding.cpu().numpy()
            z_np = z_embedding.cpu().numpy()

            emb_path = struct_dir / f"embeddings_{record.id}.npz"
            np.savez_compressed(emb_path, s=s_np, z=z_np)
            written_files['embedding_files'].append(str(emb_path))
            if verbose:
                print(f"  Wrote embeddings: {emb_path.name}")

    if verbose:
        print("\n" + "=" * 70)
        print("WRITING COMPLETE")
        print(f"  Structure files: {len(written_files['structure_files'])}")
        print(f"  Confidence files: {len(written_files['confidence_files'])}")
        print(f"  Embedding files: {len(written_files['embedding_files'])}")
        print("=" * 70 + "\n")

    return written_files

from abc import ABC, abstractmethod
from typing import Optional, Dict, Any, Set, List, Union

import torch
import numpy as np
from boltz.data import const
from boltz.model.potentials.schedules import (
    ParameterSchedule,
    ExponentialInterpolation,
    PiecewiseStepFunction,
)
from boltz.model.loss.diffusionv2 import weighted_rigid_align


class Potential(ABC):
    def __init__(
        self,
        parameters: Optional[
            Dict[str, Union[ParameterSchedule, float, int, bool]]
        ] = None,
    ):
        self.parameters = parameters

    def compute(self, coords, feats, parameters):
        index, args, com_args, ref_args, operator_args = self.compute_args(
            feats, parameters
        )

        if index.shape[1] == 0:
            return torch.zeros(coords.shape[:-2], device=coords.device)

        if com_args is not None:
            com_index, atom_pad_mask = com_args
            unpad_com_index = com_index[atom_pad_mask]
            unpad_coords = coords[..., atom_pad_mask, :]
            coords = torch.zeros(
                (*unpad_coords.shape[:-2], unpad_com_index.max() + 1, 3),
                device=coords.device,
            ).scatter_reduce(
                -2,
                unpad_com_index.unsqueeze(-1).expand_as(unpad_coords),
                unpad_coords,
                "mean",
            )
        else:
            com_index, atom_pad_mask = None, None

        if ref_args is not None:
            ref_coords, ref_mask, ref_atom_index, ref_token_index = ref_args
            coords = coords[..., ref_atom_index, :]
        else:
            ref_coords, ref_mask, ref_atom_index, ref_token_index = (
                None,
                None,
                None,
                None,
            )

        if operator_args is not None:
            negation_mask, union_index = operator_args
        else:
            negation_mask, union_index = None, None

        value = self.compute_variable(
            coords,
            index,
            ref_coords=ref_coords,
            ref_mask=ref_mask,
            compute_gradient=False,
        )
        energy = self.compute_function(
            value, *args, negation_mask=negation_mask, compute_derivative=False
        )

        if union_index is not None:
            neg_exp_energy = torch.exp(-1 * parameters["union_lambda"] * energy)
            Z = torch.zeros(
                (*energy.shape[:-1], union_index.max() + 1), device=union_index.device
            ).scatter_reduce(
                -1,
                union_index.expand_as(neg_exp_energy),
                neg_exp_energy,
                "sum",
            )
            softmax_energy = neg_exp_energy / Z[..., union_index]
            softmax_energy[Z[..., union_index] == 0] = 0
            return (energy * softmax_energy).sum(dim=-1)

        return energy.sum(dim=tuple(range(1, energy.dim())))

    def compute_gradient(self, coords, feats, parameters):
        index, args, com_args, ref_args, operator_args = self.compute_args(
            feats, parameters
        )
        if index.shape[1] == 0:
            return torch.zeros_like(coords)

        if com_args is not None:
            com_index, atom_pad_mask = com_args
            unpad_coords = coords[..., atom_pad_mask, :]
            unpad_com_index = com_index[atom_pad_mask]
            coords = torch.zeros(
                (*unpad_coords.shape[:-2], unpad_com_index.max() + 1, 3),
                device=coords.device,
            ).scatter_reduce(
                -2,
                unpad_com_index.unsqueeze(-1).expand_as(unpad_coords),
                unpad_coords,
                "mean",
            )
            com_counts = torch.bincount(com_index[atom_pad_mask])
        else:
            com_index, atom_pad_mask = None, None

        if ref_args is not None:
            ref_coords, ref_mask, ref_atom_index, ref_token_index = ref_args
            coords = coords[..., ref_atom_index, :]
        else:
            ref_coords, ref_mask, ref_atom_index, ref_token_index = (
                None,
                None,
                None,
                None,
            )

        if operator_args is not None:
            negation_mask, union_index = operator_args
        else:
            negation_mask, union_index = None, None

        value, grad_value = self.compute_variable(
            coords,
            index,
            ref_coords=ref_coords,
            ref_mask=ref_mask,
            compute_gradient=True,
        )
        energy, dEnergy = self.compute_function(
            value, 
            *args, negation_mask=negation_mask, compute_derivative=True
        )
        if union_index is not None:
            neg_exp_energy = torch.exp(-1 * parameters["union_lambda"] * energy)
            Z = torch.zeros(
                (*energy.shape[:-1], union_index.max() + 1), device=union_index.device
            ).scatter_reduce(
                -1,
                union_index.expand_as(energy),
                neg_exp_energy,
                "sum",
            )
            softmax_energy = neg_exp_energy / Z[..., union_index]
            softmax_energy[Z[..., union_index] == 0] = 0
            f = torch.zeros(
                (*energy.shape[:-1], union_index.max() + 1), device=union_index.device
            ).scatter_reduce(
                -1,
                union_index.expand_as(energy),
                energy * softmax_energy,
                "sum",
            )
            dSoftmax = (
                dEnergy
                * softmax_energy
                * (1 + parameters["union_lambda"] * (energy - f[..., union_index]))
            )
            prod = dSoftmax.tile(grad_value.shape[-3]).unsqueeze(
                -1
            ) * grad_value.flatten(start_dim=-3, end_dim=-2)
            if prod.dim() > 3:
                prod = prod.sum(dim=list(range(1, prod.dim() - 2)))
            grad_atom = torch.zeros_like(coords).scatter_reduce(
                -2,
                index.flatten(start_dim=0, end_dim=1)
                .unsqueeze(-1)
                .expand((*coords.shape[:-2], -1, 3)),
                prod,
                "sum",
            )
        else:
            prod = dEnergy.tile(grad_value.shape[-3]).unsqueeze(
                -1
            ) * grad_value.flatten(start_dim=-3, end_dim=-2)
            if prod.dim() > 3:
                prod = prod.sum(dim=list(range(1, prod.dim() - 2)))
            grad_atom = torch.zeros_like(coords).scatter_reduce(
                -2,
                index.flatten(start_dim=0, end_dim=1)
                .unsqueeze(-1)
                .expand((*coords.shape[:-2], -1, 3)),  # 9 x 516 x 3
                prod,
                "sum",
            )

        if com_index is not None:
            grad_atom = grad_atom[..., com_index, :]
        elif ref_token_index is not None:
            grad_atom = grad_atom[..., ref_token_index, :]

        return grad_atom

    def compute_parameters(self, t):
        if self.parameters is None:
            return None
        parameters = {
            name: parameter
            if not isinstance(parameter, ParameterSchedule)
            else parameter.compute(t)
            for name, parameter in self.parameters.items()
        }
        return parameters

    @abstractmethod
    def compute_function(
        self, value, *args, negation_mask=None, compute_derivative=False
    ):
        raise NotImplementedError

    @abstractmethod
    def compute_variable(self, coords, index, compute_gradient=False):
        raise NotImplementedError

    @abstractmethod
    def compute_args(self, t, feats, **parameters):
        raise NotImplementedError

    def get_reference_coords(self, feats, parameters):
        return None, None


class FlatBottomPotential(Potential):
    def compute_function(
        self,
        value,
        k,
        lower_bounds,
        upper_bounds,
        negation_mask=None,
        compute_derivative=False,
    ):
        if lower_bounds is None:
            lower_bounds = torch.full_like(value, float("-inf"))
        if upper_bounds is None:
            upper_bounds = torch.full_like(value, float("inf"))
        lower_bounds = lower_bounds.expand_as(value).clone()
        upper_bounds = upper_bounds.expand_as(value).clone()

        if negation_mask is not None:
            unbounded_below_mask = torch.isneginf(lower_bounds)
            unbounded_above_mask = torch.isposinf(upper_bounds)
            unbounded_mask = unbounded_below_mask + unbounded_above_mask
            assert torch.all(unbounded_mask + negation_mask)
            lower_bounds[~unbounded_above_mask * ~negation_mask] = upper_bounds[
                ~unbounded_above_mask * ~negation_mask
            ]
            upper_bounds[~unbounded_above_mask * ~negation_mask] = float("inf")
            upper_bounds[~unbounded_below_mask * ~negation_mask] = lower_bounds[
                ~unbounded_below_mask * ~negation_mask
            ]
            lower_bounds[~unbounded_below_mask * ~negation_mask] = float("-inf")

        neg_overflow_mask = value < lower_bounds
        pos_overflow_mask = value > upper_bounds

        energy = torch.zeros_like(value)
        energy[neg_overflow_mask] = (k * (lower_bounds - value))[neg_overflow_mask]
        energy[pos_overflow_mask] = (k * (value - upper_bounds))[pos_overflow_mask]
        if not compute_derivative:
            return energy

        dEnergy = torch.zeros_like(value)
        dEnergy[neg_overflow_mask] = (
            -1 * k.expand_as(neg_overflow_mask)[neg_overflow_mask]
        )
        dEnergy[pos_overflow_mask] = (
            1 * k.expand_as(pos_overflow_mask)[pos_overflow_mask]
        )

        return energy, dEnergy


class ReferencePotential(Potential):
    def compute_variable(
        self, coords, index, ref_coords, ref_mask, compute_gradient=False
    ):
        aligned_ref_coords = weighted_rigid_align(
            ref_coords.float(),
            coords[:, index].float(),
            ref_mask,
            ref_mask,
        )

        r = coords[:, index] - aligned_ref_coords
        r_norm = torch.linalg.norm(r, dim=-1)

        if not compute_gradient:
            return r_norm

        r_hat = r / r_norm.unsqueeze(-1)
        grad = (r_hat * ref_mask.unsqueeze(-1)).unsqueeze(1)
        return r_norm, grad


class DistancePotential(Potential):
    def compute_variable(
        self, coords, index, ref_coords=None, ref_mask=None, compute_gradient=False
    ):
        r_ij = coords.index_select(-2, index[0]) - coords.index_select(-2, index[1])
        r_ij_norm = torch.linalg.norm(r_ij, dim=-1)
        r_hat_ij = r_ij / r_ij_norm.unsqueeze(-1)

        if not compute_gradient:
            return r_ij_norm

        grad_i = r_hat_ij
        grad_j = -1 * r_hat_ij
        grad = torch.stack((grad_i, grad_j), dim=1)
        return r_ij_norm, grad


class DihedralPotential(Potential):
    def compute_variable(
        self, coords, index, ref_coords=None, ref_mask=None, compute_gradient=False
    ):
        r_ij = coords.index_select(-2, index[0]) - coords.index_select(-2, index[1])
        r_kj = coords.index_select(-2, index[2]) - coords.index_select(-2, index[1])
        r_kl = coords.index_select(-2, index[2]) - coords.index_select(-2, index[3])

        n_ijk = torch.cross(r_ij, r_kj, dim=-1)
        n_jkl = torch.cross(r_kj, r_kl, dim=-1)

        r_kj_norm = torch.linalg.norm(r_kj, dim=-1)
        n_ijk_norm = torch.linalg.norm(n_ijk, dim=-1)
        n_jkl_norm = torch.linalg.norm(n_jkl, dim=-1)

        sign_phi = torch.sign(
            r_kj.unsqueeze(-2) @ torch.cross(n_ijk, n_jkl, dim=-1).unsqueeze(-1)
        ).squeeze(-1, -2)
        phi = sign_phi * torch.arccos(
            torch.clamp(
                (n_ijk.unsqueeze(-2) @ n_jkl.unsqueeze(-1)).squeeze(-1, -2)
                / (n_ijk_norm * n_jkl_norm),
                -1 + 1e-8,
                1 - 1e-8,
            )
        )

        if not compute_gradient:
            return phi

        a = (
            (r_ij.unsqueeze(-2) @ r_kj.unsqueeze(-1)).squeeze(-1, -2) / (r_kj_norm**2)
        ).unsqueeze(-1)
        b = (
            (r_kl.unsqueeze(-2) @ r_kj.unsqueeze(-1)).squeeze(-1, -2) / (r_kj_norm**2)
        ).unsqueeze(-1)

        grad_i = n_ijk * (r_kj_norm / n_ijk_norm**2).unsqueeze(-1)
        grad_l = -1 * n_jkl * (r_kj_norm / n_jkl_norm**2).unsqueeze(-1)
        grad_j = (a - 1) * grad_i - b * grad_l
        grad_k = (b - 1) * grad_l - a * grad_i
        grad = torch.stack((grad_i, grad_j, grad_k, grad_l), dim=1)
        return phi, grad


class AbsDihedralPotential(DihedralPotential):
    def compute_variable(
        self, coords, index, ref_coords=None, ref_mask=None, compute_gradient=False
    ):
        if not compute_gradient:
            phi = super().compute_variable(
                coords, index, compute_gradient=compute_gradient
            )
            phi = torch.abs(phi)
            return phi

        phi, grad = super().compute_variable(
            coords, index, compute_gradient=compute_gradient
        )
        grad[(phi < 0)[..., None, :, None].expand_as(grad)] *= -1
        phi = torch.abs(phi)

        return phi, grad


class PoseBustersPotential(FlatBottomPotential, DistancePotential):
    def compute_args(self, feats, parameters):
        pair_index = feats["rdkit_bounds_index"][0]
        lower_bounds = feats["rdkit_lower_bounds"][0].clone()
        upper_bounds = feats["rdkit_upper_bounds"][0].clone()
        bond_mask = feats["rdkit_bounds_bond_mask"][0]
        angle_mask = feats["rdkit_bounds_angle_mask"][0]

        lower_bounds[bond_mask * ~angle_mask] *= 1.0 - parameters["bond_buffer"]
        upper_bounds[bond_mask * ~angle_mask] *= 1.0 + parameters["bond_buffer"]
        lower_bounds[~bond_mask * angle_mask] *= 1.0 - parameters["angle_buffer"]
        upper_bounds[~bond_mask * angle_mask] *= 1.0 + parameters["angle_buffer"]
        lower_bounds[bond_mask * angle_mask] *= 1.0 - min(
            parameters["bond_buffer"], parameters["angle_buffer"]
        )
        upper_bounds[bond_mask * angle_mask] *= 1.0 + min(
            parameters["bond_buffer"], parameters["angle_buffer"]
        )
        lower_bounds[~bond_mask * ~angle_mask] *= 1.0 - parameters["clash_buffer"]
        upper_bounds[~bond_mask * ~angle_mask] = float("inf")

        vdw_radii = torch.zeros(
            const.num_elements, dtype=torch.float32, device=pair_index.device
        )
        vdw_radii[1:119] = torch.tensor(
            const.vdw_radii, dtype=torch.float32, device=pair_index.device
        )
        atom_vdw_radii = (
            feats["ref_element"].float() @ vdw_radii.unsqueeze(-1)
        ).squeeze(-1)[0]
        bond_cutoffs = 0.35 + atom_vdw_radii[pair_index].mean(dim=0)
        lower_bounds[~bond_mask] = torch.max(lower_bounds[~bond_mask], bond_cutoffs[~bond_mask])
        upper_bounds[bond_mask] = torch.min(upper_bounds[bond_mask], bond_cutoffs[bond_mask])

        k = torch.ones_like(lower_bounds)

        return pair_index, (k, lower_bounds, upper_bounds), None, None, None


class ConnectionsPotential(FlatBottomPotential, DistancePotential):
    def compute_args(self, feats, parameters):
        pair_index = feats["connected_atom_index"][0]
        lower_bounds = None
        upper_bounds = torch.full(
            (pair_index.shape[1],), parameters["buffer"], device=pair_index.device
        )
        k = torch.ones_like(upper_bounds)

        return pair_index, (k, lower_bounds, upper_bounds), None, None, None


class VDWOverlapPotential(FlatBottomPotential, DistancePotential):
    def compute_args(self, feats, parameters):
        atom_chain_id = (
            torch.bmm(
                feats["atom_to_token"].float(), feats["asym_id"].unsqueeze(-1).float()
            )
            .squeeze(-1)
            .long()
        )[0]
        atom_pad_mask = feats["atom_pad_mask"][0].bool()
        chain_sizes = torch.bincount(atom_chain_id[atom_pad_mask])
        single_ion_mask = (chain_sizes > 1)[atom_chain_id]

        vdw_radii = torch.zeros(
            const.num_elements, dtype=torch.float32, device=atom_chain_id.device
        )
        vdw_radii[1:119] = torch.tensor(
            const.vdw_radii, dtype=torch.float32, device=atom_chain_id.device
        )
        atom_vdw_radii = (
            feats["ref_element"].float() @ vdw_radii.unsqueeze(-1)
        ).squeeze(-1)[0]

        pair_index = torch.triu_indices(
            atom_chain_id.shape[0],
            atom_chain_id.shape[0],
            1,
            device=atom_chain_id.device,
        )

        pair_pad_mask = atom_pad_mask[pair_index].all(dim=0)
        pair_ion_mask = single_ion_mask[pair_index[0]] * single_ion_mask[pair_index[1]]

        num_chains = atom_chain_id.max() + 1
        connected_chain_index = feats["connected_chain_index"][0]
        connected_chain_matrix = torch.eye(
            num_chains, device=atom_chain_id.device, dtype=torch.bool
        )
        connected_chain_matrix[connected_chain_index[0], connected_chain_index[1]] = (
            True
        )
        connected_chain_matrix[connected_chain_index[1], connected_chain_index[0]] = (
            True
        )
        connected_chain_mask = connected_chain_matrix[
            atom_chain_id[pair_index[0]], atom_chain_id[pair_index[1]]
        ]

        pair_index = pair_index[
            :, pair_pad_mask * pair_ion_mask * ~connected_chain_mask
        ]

        lower_bounds = atom_vdw_radii[pair_index].sum(dim=0) * (
            1.0 - parameters["buffer"]
        )
        upper_bounds = None
        k = torch.ones_like(lower_bounds)

        return pair_index, (k, lower_bounds, upper_bounds), None, None, None


class SymmetricChainCOMPotential(FlatBottomPotential, DistancePotential):
    def compute_args(self, feats, parameters):
        atom_chain_id = (
            torch.bmm(
                feats["atom_to_token"].float(), feats["asym_id"].unsqueeze(-1).float()
            )
            .squeeze(-1)
            .long()
        )[0]
        atom_pad_mask = feats["atom_pad_mask"][0].bool()
        chain_sizes = torch.bincount(atom_chain_id[atom_pad_mask])
        single_ion_mask = chain_sizes > 1

        pair_index = feats["symmetric_chain_index"][0]
        pair_ion_mask = single_ion_mask[pair_index[0]] * single_ion_mask[pair_index[1]]
        pair_index = pair_index[:, pair_ion_mask]
        lower_bounds = torch.full(
            (pair_index.shape[1],),
            parameters["buffer"],
            dtype=torch.float32,
            device=pair_index.device,
        )
        upper_bounds = None
        k = torch.ones_like(lower_bounds)

        return (
            pair_index,
            (k, lower_bounds, upper_bounds),
            (atom_chain_id, atom_pad_mask),
            None,
            None,
        )


class StereoBondPotential(FlatBottomPotential, AbsDihedralPotential):
    def compute_args(self, feats, parameters):
        stereo_bond_index = feats["stereo_bond_index"][0]
        stereo_bond_orientations = feats["stereo_bond_orientations"][0].bool()

        lower_bounds = torch.zeros(
            stereo_bond_orientations.shape, device=stereo_bond_orientations.device
        )
        upper_bounds = torch.zeros(
            stereo_bond_orientations.shape, device=stereo_bond_orientations.device
        )
        lower_bounds[stereo_bond_orientations] = torch.pi - parameters["buffer"]
        upper_bounds[stereo_bond_orientations] = float("inf")
        lower_bounds[~stereo_bond_orientations] = float("-inf")
        upper_bounds[~stereo_bond_orientations] = parameters["buffer"]

        k = torch.ones_like(lower_bounds)

        return stereo_bond_index, (k, lower_bounds, upper_bounds), None, None, None


class ChiralAtomPotential(FlatBottomPotential, DihedralPotential):
    def compute_args(self, feats, parameters):
        chiral_atom_index = feats["chiral_atom_index"][0]
        chiral_atom_orientations = feats["chiral_atom_orientations"][0].bool()

        lower_bounds = torch.zeros(
            chiral_atom_orientations.shape, device=chiral_atom_orientations.device
        )
        upper_bounds = torch.zeros(
            chiral_atom_orientations.shape, device=chiral_atom_orientations.device
        )
        lower_bounds[chiral_atom_orientations] = parameters["buffer"]
        upper_bounds[chiral_atom_orientations] = float("inf")
        upper_bounds[~chiral_atom_orientations] = -1 * parameters["buffer"]
        lower_bounds[~chiral_atom_orientations] = float("-inf")

        k = torch.ones_like(lower_bounds)
        return chiral_atom_index, (k, lower_bounds, upper_bounds), None, None, None


class PlanarBondPotential(FlatBottomPotential, AbsDihedralPotential):
    def compute_args(self, feats, parameters):
        double_bond_index = feats["planar_bond_index"][0].T
        double_bond_improper_index = torch.tensor(
            [
                [1, 2, 3, 0],
                [4, 5, 0, 3],
            ],
            device=double_bond_index.device,
        ).T
        improper_index = (
            double_bond_index[:, double_bond_improper_index]
            .swapaxes(0, 1)
            .flatten(start_dim=1)
        )
        lower_bounds = None
        upper_bounds = torch.full(
            (improper_index.shape[1],),
            parameters["buffer"],
            device=improper_index.device,
        )
        k = torch.ones_like(upper_bounds)

        return improper_index, (k, lower_bounds, upper_bounds), None, None, None


class TemplateReferencePotential(FlatBottomPotential, ReferencePotential):
    def compute_args(self, feats, parameters):
        if "template_mask_cb" not in feats or "template_force" not in feats:
            return torch.empty([1, 0]), None, None, None, None

        template_mask = feats["template_mask_cb"][feats["template_force"]]
        if template_mask.shape[0] == 0:
            return torch.empty([1, 0]), None, None, None, None

        ref_coords = feats["template_cb"][feats["template_force"]].clone()
        ref_mask = feats["template_mask_cb"][feats["template_force"]].clone()
        ref_atom_index = (
            torch.bmm(
                feats["token_to_rep_atom"].float(),
                torch.arange(
                    feats["atom_pad_mask"].shape[1],
                    device=feats["atom_pad_mask"].device,
                    dtype=torch.float32,
                )[None, :, None],
            )
            .squeeze(-1)
            .long()
        )[0]
        ref_token_index = (
            torch.bmm(
                feats["atom_to_token"].float(),
                feats["token_index"].unsqueeze(-1).float(),
            )
            .squeeze(-1)
            .long()
        )[0]

        index = torch.arange(
            template_mask.shape[-1], dtype=torch.long, device=template_mask.device
        )[None]
        upper_bounds = torch.full(
            template_mask.shape, float("inf"), device=index.device, dtype=torch.float32
        )
        ref_idxs = torch.argwhere(template_mask).T
        upper_bounds[ref_idxs.unbind()] = feats["template_force_threshold"][
            feats["template_force"]
        ][ref_idxs[0]]

        lower_bounds = None
        k = torch.ones_like(upper_bounds)
        return (
            index,
            (k, lower_bounds, upper_bounds),
            None,
            (ref_coords, ref_mask, ref_atom_index, ref_token_index),
            None,
        )


class ContactPotentital(FlatBottomPotential, DistancePotential):
    def compute_args(self, feats, parameters):
        index = feats["contact_pair_index"][0]
        union_index = feats["contact_union_index"][0]
        negation_mask = feats["contact_negation_mask"][0]
        lower_bounds = None
        upper_bounds = feats["contact_thresholds"][0].clone()
        k = torch.ones_like(upper_bounds)
        return (
            index,
            (k, lower_bounds, upper_bounds),
            None,
            None,
            (negation_mask, union_index),
        )


class CDR3ConformationPotential(FlatBottomPotential, DihedralPotential):
    """Potential to control CDR3 loop conformations via backbone psi dihedral angles.

    This potential allows steering CDR3 regions toward specific conformations by
    constraining psi angles (N-CA-C-N') to target ranges. Common CDR3 conformations:

    - extended: psi ≈ 140° (2.44 rad) - elongated loop
    - compact: psi ≈ -40° (-0.70 rad) - tight turn
    - kinked: psi ≈ 0° - sharp bend (often at apex)
    - custom: user-specified angle ranges

    Features expected in feats:
    - cdr3_dihedral_index: [4, N] tensor of atom indices (N, CA, C, N') for psi dihedrals
    - cdr3_target_lower: [N] tensor of lower angle bounds (radians)
    - cdr3_target_upper: [N] tensor of upper angle bounds (radians)
    """

    def compute_args(self, feats, parameters):
        # Get device from feats for consistent tensor placement
        device = feats["atom_pad_mask"].device

        # Check if CDR3 features are present
        if "cdr3_dihedral_index" not in feats:
            return torch.empty([4, 0], dtype=torch.long, device=device), (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        dihedral_index = feats["cdr3_dihedral_index"][0]  # [4, N_dihedrals]

        if dihedral_index.shape[1] == 0:
            return dihedral_index, (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        # Get target angle bounds
        lower_bounds = feats["cdr3_target_lower"][0].clone()  # [N_dihedrals]
        upper_bounds = feats["cdr3_target_upper"][0].clone()  # [N_dihedrals]

        # Apply buffer from parameters (allows some flexibility)
        if "buffer" in parameters:
            lower_bounds = lower_bounds - parameters["buffer"]
            upper_bounds = upper_bounds + parameters["buffer"]

        # Spring constant for energy penalty
        k = torch.ones_like(lower_bounds)

        return dihedral_index, (k, lower_bounds, upper_bounds), None, None, None


class AntigenOrientationPotential(FlatBottomPotential, DistancePotential):
    """Potential to optimize antigen orientation relative to CDR loops.

    This potential encourages the antigen to orient such that its residues
    are in contact with CDR loops from heavy and light chains. It works by:
    1. Computing distances between antigen CA atoms and CDR CA atoms
    2. Applying a flat-bottom energy penalty for distances > contact_threshold
    3. Providing gradients that push antigen atoms toward nearest CDR atoms

    Features expected in feats:
    - antigen_atom_index: [N_antigen] tensor of antigen CA atom indices
    - cdr_atom_index: [N_cdr] tensor of CDR CA atom indices
    - antigen_orientation_threshold: scalar contact threshold (Angstrom)
    """

    def compute_args(self, feats, parameters):
        """Extract antigen and CDR atom indices and compute pair indices.

        Returns pair indices for all antigen-CDR atom pairs along with
        distance thresholds for the flat-bottom potential.
        """
        device = feats["atom_pad_mask"].device

        # Check if antigen orientation features are present
        if "antigen_atom_index" not in feats:
            return torch.empty([2, 0], dtype=torch.long, device=device), (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        antigen_atom_index = feats["antigen_atom_index"][0]  # [N_antigen]
        cdr_atom_index = feats["cdr_atom_index"][0]  # [N_cdr]

        if antigen_atom_index.shape[0] == 0 or cdr_atom_index.shape[0] == 0:
            return torch.empty([2, 0], dtype=torch.long, device=device), (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        # Get contact threshold
        threshold = feats["antigen_orientation_threshold"][0].item()

        # Create pair indices: for each antigen atom, find distance to nearest CDR atom
        # We create all antigen-CDR pairs and use union_index to group by antigen atom
        n_antigen = antigen_atom_index.shape[0]
        n_cdr = cdr_atom_index.shape[0]

        # Create all pairs: antigen_atom_i <-> cdr_atom_j
        antigen_expanded = antigen_atom_index.unsqueeze(1).expand(-1, n_cdr).flatten()
        cdr_expanded = cdr_atom_index.unsqueeze(0).expand(n_antigen, -1).flatten()

        pair_index = torch.stack([antigen_expanded, cdr_expanded], dim=0)  # [2, N_antigen * N_cdr]

        # Union index groups pairs by antigen atom for soft-min computation
        union_index = torch.arange(n_antigen, device=device).unsqueeze(1).expand(-1, n_cdr).flatten()

        # Upper bounds = contact_threshold for all pairs
        upper_bounds = torch.full(
            (pair_index.shape[1],), threshold, dtype=torch.float32, device=device
        )
        lower_bounds = None

        # Spring constant
        k = torch.ones_like(upper_bounds)

        # No negation (we want distance < threshold)
        negation_mask = torch.zeros(pair_index.shape[1], dtype=torch.bool, device=device)

        return (
            pair_index,
            (k, lower_bounds, upper_bounds),
            None,
            None,
            (negation_mask, union_index),
        )


class EmbeddingInterfacePotential(FlatBottomPotential, DistancePotential):
    """Potential that uses pair embeddings to weight distance-based interface steering.

    Instead of treating all CDR-antigen pairs equally, this potential extracts
    pair embedding norms from the model's pair representation z to weight the
    distance-based contact potential. Pairs with stronger embeddings (higher norms)
    get more steering force.

    Features expected in feats:
    - embedding_interface_antigen_atom_idx: [N_antigen] tensor of antigen CA atom indices
    - embedding_interface_cdr_atom_idx: [N_cdr] tensor of CDR CA atom indices
    - embedding_interface_antigen_token_idx: [N_antigen] tensor of antigen token indices
    - embedding_interface_cdr_token_idx: [N_cdr] tensor of CDR token indices
    - embedding_interface_threshold: scalar contact threshold (Angstrom)
    """

    def __init__(self, parameters):
        super().__init__(parameters)
        self._embedding_weights = None

    def set_embedding_weights(self, z_trunk, feats):
        """Extract CDR-antigen pair embedding norms as importance weights.

        Parameters
        ----------
        z_trunk : Tensor
            [batch, N_tokens, N_tokens, token_z] pair embeddings from trunk
        feats : dict
            Feature dictionary containing token indices
        """
        if "embedding_interface_antigen_token_idx" not in feats:
            self._embedding_weights = None
            return

        ag_idx = feats["embedding_interface_antigen_token_idx"][0]  # [N_ag]
        cdr_idx = feats["embedding_interface_cdr_token_idx"][0]     # [N_cdr]

        if ag_idx.shape[0] == 0 or cdr_idx.shape[0] == 0:
            self._embedding_weights = None
            return

        # Extract CDR-antigen sub-matrix from pair embeddings
        # z_trunk: [batch, N_tokens, N_tokens, token_z]
        # We want z[batch, ag_i, cdr_j, :] for all ag_i, cdr_j pairs
        z_ag = z_trunk[:, ag_idx]            # [batch, N_ag, N_tokens, token_z]
        z_interface = z_ag[:, :, cdr_idx]    # [batch, N_ag, N_cdr, token_z]

        # Compute embedding norms as importance weights: [N_ag, N_cdr]
        weights = z_interface.norm(dim=-1).mean(dim=0)  # average over batch

        # Normalize to [0, 1] range
        w_min = weights.min()
        w_max = weights.max()
        weights = (weights - w_min) / (w_max - w_min + 1e-8)

        # Flatten to match pair ordering (antigen_i, cdr_j)
        self._embedding_weights = weights.flatten()  # [N_ag * N_cdr]

    def compute_args(self, feats, parameters):
        """Extract antigen and CDR atom indices and compute pair indices.

        Same pair construction as AntigenOrientationPotential but with
        spring constants weighted by pair embedding magnitudes.
        """
        device = feats["atom_pad_mask"].device

        if "embedding_interface_antigen_atom_idx" not in feats:
            return torch.empty([2, 0], dtype=torch.long, device=device), (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        antigen_atom_index = feats["embedding_interface_antigen_atom_idx"][0]  # [N_antigen]
        cdr_atom_index = feats["embedding_interface_cdr_atom_idx"][0]  # [N_cdr]

        if antigen_atom_index.shape[0] == 0 or cdr_atom_index.shape[0] == 0:
            return torch.empty([2, 0], dtype=torch.long, device=device), (
                torch.empty([0], dtype=torch.float32, device=device),
                None,
                None,
            ), None, None, None

        # Get contact threshold
        threshold = feats["embedding_interface_threshold"][0].item()

        # Create pair indices: all antigen-CDR pairs
        n_antigen = antigen_atom_index.shape[0]
        n_cdr = cdr_atom_index.shape[0]

        antigen_expanded = antigen_atom_index.unsqueeze(1).expand(-1, n_cdr).flatten()
        cdr_expanded = cdr_atom_index.unsqueeze(0).expand(n_antigen, -1).flatten()

        pair_index = torch.stack([antigen_expanded, cdr_expanded], dim=0)  # [2, N_antigen * N_cdr]

        # Union index groups pairs by antigen atom for soft-min computation
        union_index = torch.arange(n_antigen, device=device).unsqueeze(1).expand(-1, n_cdr).flatten()

        # Upper bounds = contact_threshold for all pairs
        upper_bounds = torch.full(
            (pair_index.shape[1],), threshold, dtype=torch.float32, device=device
        )
        lower_bounds = None

        # Spring constant weighted by embedding magnitudes
        k = torch.ones_like(upper_bounds)
        if self._embedding_weights is not None:
            embed_w = self._embedding_weights.to(device)
            if embed_w.shape[0] == k.shape[0]:
                k = k * embed_w

        # No negation (we want distance < threshold)
        negation_mask = torch.zeros(pair_index.shape[1], dtype=torch.bool, device=device)

        return (
            pair_index,
            (k, lower_bounds, upper_bounds),
            None,
            None,
            (negation_mask, union_index),
        )


def get_potentials(steering_args, boltz2=False):
    potentials = []
    if steering_args["fk_steering"] or steering_args["physical_guidance_update"]:
        potentials.extend(
            [
                SymmetricChainCOMPotential(
                    parameters={
                        "guidance_interval": 4,
                        "guidance_weight": 0.5
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 0.5,
                        "buffer": ExponentialInterpolation(
                            start=1.0, end=5.0, alpha=-2.0
                        ),
                    }
                ),
                VDWOverlapPotential(
                    parameters={
                        "guidance_interval": 5,
                        "guidance_weight": (
                            PiecewiseStepFunction(thresholds=[0.4], values=[0.125, 0.0])
                            if steering_args["physical_guidance_update"]
                            else 0.0
                        ),
                        "resampling_weight": PiecewiseStepFunction(
                            thresholds=[0.6], values=[0.01, 0.0]
                        ),
                        "buffer": 0.225,
                    }
                ),
                ConnectionsPotential(
                    parameters={
                        "guidance_interval": 1,
                        "guidance_weight": 0.15
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 1.0,
                        "buffer": 2.0,
                    }
                ),
                PoseBustersPotential(
                    parameters={
                        "guidance_interval": 1,
                        "guidance_weight": 0.01
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 0.1,
                        "bond_buffer": 0.125,
                        "angle_buffer": 0.125,
                        "clash_buffer": 0.10,
                    }
                ),
                ChiralAtomPotential(
                    parameters={
                        "guidance_interval": 1,
                        "guidance_weight": 0.1
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 1.0,
                        "buffer": 0.52360,
                    }
                ),
                StereoBondPotential(
                    parameters={
                        "guidance_interval": 1,
                        "guidance_weight": 0.05
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 1.0,
                        "buffer": 0.52360,
                    }
                ),
                PlanarBondPotential(
                    parameters={
                        "guidance_interval": 1,
                        "guidance_weight": 0.05
                        if steering_args["physical_guidance_update"]
                        else 0.0,
                        "resampling_weight": 1.0,
                        "buffer": 0.26180,
                    }
                ),
            ]
        )
    if boltz2 and (
        steering_args["fk_steering"] or steering_args["contact_guidance_update"]
    ):
        potentials.extend(
            [
                ContactPotentital(
                    parameters={
                        "guidance_interval": 4,
                        "guidance_weight": (
                            PiecewiseStepFunction(
                                thresholds=[0.25, 0.75], values=[0.0, 0.5, 1.0]
                            )
                            if steering_args["contact_guidance_update"]
                            else 0.0
                        ),
                        "resampling_weight": 1.0,
                        "union_lambda": ExponentialInterpolation(
                            start=8.0, end=0.0, alpha=-2.0
                        ),
                    }
                ),
                TemplateReferencePotential(
                    parameters={
                        "guidance_interval": 2,
                        "guidance_weight": 0.1
                        if steering_args["contact_guidance_update"]
                        else 0.0,
                        "resampling_weight": 1.0,
                    }
                ),
            ]
        )
    # Add CDR3 conformation potential if enabled
    if boltz2 and steering_args.get("cdr3_steering", False):
        potentials.append(
            CDR3ConformationPotential(
                parameters={
                    "guidance_interval": 2,
                    "guidance_weight": (
                        PiecewiseStepFunction(
                            thresholds=[0.3, 0.7], values=[0.5, 1.0, 0.5]
                        )
                        if steering_args["physical_guidance_update"]
                        else 0.0
                    ),
                    "resampling_weight": 1.0,
                    "buffer": 0.35,  # ~20 degrees flexibility
                }
            )
        )
    # Add antigen orientation potential if enabled
    if boltz2 and steering_args.get("antigen_steering", False):
        potentials.append(
            AntigenOrientationPotential(
                parameters={
                    "guidance_interval": 2,
                    "guidance_weight": PiecewiseStepFunction(
                        thresholds=[0.3, 0.7],
                        values=[1.5, 1.0, 0.3]  # Strong early, weaker late
                    ),
                    "resampling_weight": PiecewiseStepFunction(
                        thresholds=[0.5],
                        values=[1.0, 0.5]  # Heavy resampling early
                    ),
                    "union_lambda": ExponentialInterpolation(
                        start=8.0, end=0.0, alpha=-2.0
                    ),
                }
            )
        )
    # Add embedding interface potential if enabled
    if boltz2 and steering_args.get("embedding_interface_steering", False):
        potentials.append(
            EmbeddingInterfacePotential(
                parameters={
                    "guidance_interval": 2,
                    "guidance_weight": PiecewiseStepFunction(
                        thresholds=[0.3, 0.7],
                        values=[1.5, 1.0, 0.3]  # Strong early, weaker late
                    ),
                    "resampling_weight": PiecewiseStepFunction(
                        thresholds=[0.5],
                        values=[1.0, 0.5]  # Heavy resampling early
                    ),
                    "union_lambda": ExponentialInterpolation(
                        start=8.0, end=0.0, alpha=-2.0
                    ),
                }
            )
        )
    return potentials

"""A+.2: Composite re-ranking for FK particle steering.

Instead of relying solely on model confidence (iPTM/pLDDT) for model selection,
this module combines confidence, FK steering energy, and interface pLDDT into a
composite score. FK steering energy correlates with docking quality better than
raw confidence for steered predictions.
"""

from __future__ import annotations

from typing import Optional

import torch


def composite_rerank(
    confidence_scores: torch.Tensor,
    fk_energies: Optional[torch.Tensor] = None,
    interface_plddts: Optional[torch.Tensor] = None,
    alpha: float = 0.5,
    beta: float = 0.3,
    gamma: float = 0.2,
) -> torch.Tensor:
    """Compute composite re-ranking scores for FK-steered predictions.

    Parameters
    ----------
    confidence_scores : torch.Tensor
        Shape [N], model confidence scores (higher = better).
    fk_energies : torch.Tensor, optional
        Shape [N], FK potential energies (lower = better).
    interface_plddts : torch.Tensor, optional
        Shape [N], interface pLDDT scores (higher = better).
    alpha : float
        Weight for confidence score component.
    beta : float
        Weight for FK energy component (negated: lower energy = higher score).
    gamma : float
        Weight for interface pLDDT component.

    Returns
    -------
    torch.Tensor
        Shape [N], composite scores (higher = better).
    """
    score = alpha * confidence_scores

    if fk_energies is not None:
        # Normalize energies to [0, 1] range before combining
        e_min = fk_energies.min()
        e_max = fk_energies.max()
        if e_max - e_min > 1e-8:
            normalized_energies = (fk_energies - e_min) / (e_max - e_min)
        else:
            normalized_energies = torch.zeros_like(fk_energies)
        # Negate: lower energy = higher score
        score = score - beta * normalized_energies

    if interface_plddts is not None:
        score = score + gamma * interface_plddts

    return score


def rerank_predictions(
    predictions: list[dict],
    alpha: float = 0.5,
    beta: float = 0.3,
    gamma: float = 0.2,
) -> list[dict]:
    """Re-rank a list of prediction dicts by composite score.

    Each dict should have:
        - 'confidence_score': float
        - 'fk_energy': float (optional, lower = better)
        - 'interface_plddt': float (optional)

    Parameters
    ----------
    predictions : list[dict]
        List of prediction result dictionaries.
    alpha, beta, gamma : float
        Weights for confidence, energy, and interface pLDDT.

    Returns
    -------
    list[dict]
        Sorted predictions with 'composite_score' added, best first.
    """
    if not predictions:
        return predictions

    conf = torch.tensor([p["confidence_score"] for p in predictions])
    energies = None
    if all("fk_energy" in p for p in predictions):
        energies = torch.tensor([p["fk_energy"] for p in predictions])
    plddts = None
    if all("interface_plddt" in p for p in predictions):
        plddts = torch.tensor([p["interface_plddt"] for p in predictions])

    scores = composite_rerank(conf, energies, plddts, alpha, beta, gamma)

    for p, s in zip(predictions, scores):
        p["composite_score"] = s.item()

    return sorted(predictions, key=lambda p: p["composite_score"], reverse=True)

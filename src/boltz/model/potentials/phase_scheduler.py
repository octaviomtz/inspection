"""
Phase scheduler for Strategy G+: Progressive refinement with contact-triggered phase transitions.

Implements adaptive phase transitions based on contact score monitoring and smooth β scheduling.
"""

from typing import Optional, Tuple
from collections import deque

import torch
import numpy as np


class AdaptivePhaseScheduler:
    """
    Manages adaptive phase transitions for Strategy G+ based on contact score improvement.

    Three phases:
    - Phase 1 (Exploration): High CDR3 flexibility, weak contact guidance
    - Phase 2 (Transition): Neutral CDR3, moderate guidance
    - Phase 3 (Optimization): CDR3 optimization, strong guidance
    """

    def __init__(
        self,
        improvement_threshold: float = 0.05,
        improvement_window: int = 10,
        initial_phase: int = 1,
        verbose: bool = False,
    ):
        """
        Initialize the phase scheduler.

        Parameters
        ----------
        improvement_threshold : float
            Minimum improvement rate (5% by default) to stay in current phase.
            When improvement drops below this for N consecutive steps, advance phase.
        improvement_window : int
            Window size for computing rolling average of contact scores.
        initial_phase : int
            Starting phase (1, 2, or 3).
        verbose : bool
            Print phase transitions.
        """
        self.improvement_threshold = improvement_threshold
        self.improvement_window = improvement_window
        self.verbose = verbose

        self.current_phase = initial_phase
        self.contact_history = deque(maxlen=improvement_window)
        self.phase_start_step = 0
        self.steps_in_current_phase = 0
        self.consecutive_low_improvement = 0

    def get_phase_parameters(
        self,
        current_contact_score: float,
        timestep: int,
        num_sampling_steps: int,
    ) -> Tuple[int, float, float]:
        """
        Get current phase and associated parameters based on contact score.

        Parameters
        ----------
        current_contact_score : float
            Current contact score between CDR and antigen.
        timestep : int
            Current diffusion step (0 to num_sampling_steps).
        num_sampling_steps : int
            Total number of diffusion steps.

        Returns
        -------
        phase : int
            Current phase (1, 2, or 3).
        guidance_weight : float
            Weight for contact-based guidance (0.0 to 1.0).
        cdr3_beta : float
            β value for CDR3-specific pair scaling.
        """
        # Track contact score history
        self.contact_history.append(current_contact_score)
        self.steps_in_current_phase += 1

        # Check for phase advancement
        if len(self.contact_history) >= self.improvement_window and self.current_phase < 3:
            improvement = self._compute_improvement()

            if improvement < self.improvement_threshold:
                self.consecutive_low_improvement += 1
            else:
                self.consecutive_low_improvement = 0

            # Advance phase after consecutive low-improvement steps
            if self.consecutive_low_improvement >= 2 and self.current_phase < 3:
                self.current_phase += 1
                self.phase_start_step = timestep
                self.steps_in_current_phase = 0
                self.consecutive_low_improvement = 0
                if self.verbose:
                    print(f"Advancing to phase {self.current_phase} at step {timestep}")

        # Compute phase progress (0 to 1 within current phase)
        phase_progress = self._compute_phase_progress(timestep, num_sampling_steps)

        # Get phase-specific parameters with smooth transitions
        guidance_weight = self._smooth_guidance_weight(phase_progress)
        cdr3_beta = self._smooth_cdr3_beta(phase_progress)

        return self.current_phase, guidance_weight, cdr3_beta

    def _compute_improvement(self) -> float:
        """
        Compute contact score improvement rate over the window.

        Returns
        -------
        float
            Improvement rate: (current - oldest) / oldest, clamped to [0, 1].
        """
        if len(self.contact_history) < 2:
            return 1.0  # No improvement data yet, stay in phase

        oldest = self.contact_history[0]
        current = self.contact_history[-1]

        if oldest <= 0:
            return 1.0  # Avoid division by zero

        improvement = (current - oldest) / abs(oldest)
        return max(0.0, improvement)  # Clamp to [0, ∞)

    def _compute_phase_progress(self, timestep: int, num_sampling_steps: int) -> float:
        """
        Compute progress within current phase for smooth transitions.

        Returns
        -------
        float
            Phase progress in [0, 1], where 0 = phase start, 1 = phase end.
        """
        # Estimate phase duration as ~1/3 of total steps
        phase_duration = max(1, num_sampling_steps // 3)
        progress = self.steps_in_current_phase / phase_duration
        return min(1.0, progress)

    def _smooth_guidance_weight(self, phase_progress: float) -> float:
        """
        Compute smooth guidance weight based on current phase and progress.

        Parameters
        ----------
        phase_progress : float
            Progress within current phase [0, 1].

        Returns
        -------
        float
            Guidance weight [0.1, 0.5] with smooth ramps between phases.
        """
        if self.current_phase == 1:
            # Phase 1: Start with 0.1, gradually increase to 0.25
            return 0.1 + 0.15 * phase_progress
        elif self.current_phase == 2:
            # Phase 2: 0.25 to 0.4 (transition)
            return 0.25 + 0.15 * phase_progress
        else:  # Phase 3
            # Phase 3: 0.4 to 0.5 (optimization)
            return 0.4 + 0.1 * phase_progress

    def _smooth_cdr3_beta(self, phase_progress: float) -> float:
        """
        Compute smooth CDR3 β value based on current phase and progress.

        Parameters
        ----------
        phase_progress : float
            Progress within current phase [0, 1].

        Returns
        -------
        float
            CDR3 β value [-0.5, 0.0] following phase progression.
        """
        if self.current_phase == 1:
            # Phase 1: Exploratory β, from -0.5 to -0.2
            return -0.5 + 0.3 * phase_progress
        elif self.current_phase == 2:
            # Phase 2: Transition β, from -0.2 to 0.0
            return -0.2 + 0.2 * phase_progress
        else:  # Phase 3
            # Phase 3: Neutral to slight optimization, stay near 0
            return 0.0  # Could add small positive β for final optimization

    def reset(self):
        """Reset scheduler for a new prediction."""
        self.current_phase = 1
        self.contact_history.clear()
        self.phase_start_step = 0
        self.steps_in_current_phase = 0
        self.consecutive_low_improvement = 0


def smooth_linear_interpolation(
    x: float,
    x_min: float = 0.0,
    x_max: float = 1.0,
    y_min: float = 0.0,
    y_max: float = 1.0,
) -> float:
    """
    Smooth linear interpolation between two values.

    Parameters
    ----------
    x : float
        Input value.
    x_min, x_max : float
        Input range.
    y_min, y_max : float
        Output range.

    Returns
    -------
    float
        Interpolated value.
    """
    if x_max == x_min:
        return y_min
    t = (x - x_min) / (x_max - x_min)
    t = max(0.0, min(1.0, t))  # Clamp to [0, 1]
    return y_min + t * (y_max - y_min)

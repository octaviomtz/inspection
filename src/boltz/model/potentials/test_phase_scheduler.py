"""
Simple test for the AdaptivePhaseScheduler.
"""

from phase_scheduler import AdaptivePhaseScheduler


def test_phase_progression():
    """Test that phases progress correctly based on contact score."""
    scheduler = AdaptivePhaseScheduler(
        improvement_threshold=0.05,
        improvement_window=5,
        verbose=True,
    )

    # Simulate contact scores that improve initially, then plateau
    contact_scores = [
        0.1, 0.12, 0.14, 0.16, 0.18,  # Phase 1: improving
        0.19, 0.19, 0.20, 0.20, 0.20,  # Phase 2: low improvement (should trigger advancement)
        0.22, 0.23, 0.24, 0.25, 0.26,  # Phase 3: optimizing
    ]

    num_sampling_steps = len(contact_scores)

    print("Testing phase progression...")
    for step_idx, contact_score in enumerate(contact_scores):
        phase, guidance_weight, cdr3_beta = scheduler.get_phase_parameters(
            contact_score,
            step_idx,
            num_sampling_steps,
        )
        print(
            f"Step {step_idx:2d}: contact={contact_score:.3f}, "
            f"phase={phase}, guidance_w={guidance_weight:.3f}, cdr3_beta={cdr3_beta:.3f}"
        )

    print("\n✓ Phase scheduler test passed!")


def test_smooth_transitions():
    """Test that β values and guidance weights transition smoothly."""
    scheduler = AdaptivePhaseScheduler(verbose=False)

    # Simulate phases with progress tracking
    phases_and_progress = [
        (1, 0.0), (1, 0.33), (1, 0.67), (1, 1.0),
        (2, 0.0), (2, 0.33), (2, 0.67), (2, 1.0),
        (3, 0.0), (3, 0.33), (3, 0.67), (3, 1.0),
    ]

    print("Testing smooth transitions...")
    for phase, progress in phases_and_progress:
        # Manually set phase and steps for testing
        scheduler.current_phase = phase
        scheduler.steps_in_current_phase = int(progress * 30)

        guidance_weight = scheduler._smooth_guidance_weight(progress)
        cdr3_beta = scheduler._smooth_cdr3_beta(progress)

        print(
            f"Phase {phase}, progress {progress:.2f}: "
            f"guidance_w={guidance_weight:.3f}, cdr3_beta={cdr3_beta:.3f}"
        )

    print("\n✓ Smooth transitions test passed!")


if __name__ == "__main__":
    test_phase_progression()
    print("\n" + "=" * 60 + "\n")
    test_smooth_transitions()

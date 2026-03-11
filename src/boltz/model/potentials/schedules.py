import math
from abc import ABC


class ParameterSchedule(ABC):
    def compute(self, t):
        raise NotImplementedError


class ExponentialInterpolation(ParameterSchedule):
    def __init__(self, start, end, alpha):
        self.start = start
        self.end = end
        self.alpha = alpha

    def compute(self, t):
        if self.alpha != 0:
            return self.start + (self.end - self.start) * (
                math.exp(self.alpha * t) - 1
            ) / (math.exp(self.alpha) - 1)
        else:
            return self.start + (self.end - self.start) * t


class PiecewiseStepFunction(ParameterSchedule):
    def __init__(self, thresholds, values):
        self.thresholds = thresholds
        self.values = values

    def compute(self, t):
        assert len(self.thresholds) > 0
        assert len(self.values) == len(self.thresholds) + 1

        idx = 0
        while idx < len(self.thresholds) and t > self.thresholds[idx]:
            idx += 1
        return self.values[idx]


class CosineDecaySchedule(ParameterSchedule):
    """Cosine decay from start_val to end_val over a time range.

    t is steering_t (1.0 at beginning of diffusion, 0.0 at end).
    Returns start_val when t >= t_start, end_val when t <= t_end,
    and a smooth cosine interpolation in between.
    """

    def __init__(self, start_val, end_val, t_start, t_end):
        self.start_val = start_val
        self.end_val = end_val
        self.t_start = t_start
        self.t_end = t_end

    def compute(self, t):
        if t >= self.t_start:
            return self.start_val
        elif t <= self.t_end:
            return self.end_val
        else:
            progress = (self.t_start - t) / (self.t_start - self.t_end)
            return self.end_val + 0.5 * (self.start_val - self.end_val) * (
                1 + math.cos(math.pi * progress)
            )


class LinearRampSchedule(ParameterSchedule):
    """Linear ramp from start_val to end_val over a time range.

    t is steering_t (1.0 at beginning of diffusion, 0.0 at end).
    """

    def __init__(self, start_val, end_val, t_start, t_end):
        self.start_val = start_val
        self.end_val = end_val
        self.t_start = t_start
        self.t_end = t_end

    def compute(self, t):
        if t >= self.t_start:
            return self.start_val
        elif t <= self.t_end:
            return self.end_val
        else:
            progress = (self.t_start - t) / (self.t_start - self.t_end)
            return self.start_val + (self.end_val - self.start_val) * progress

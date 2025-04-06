import numpy as np
# from src.Reservoir import Reservoir
from src.OnlineMeanVar import OnlineMeanVar

class Red:
    """
    Represents read event distribution attributes, using an online mean-variance tracker.
    """

    def __init__(self, initlen: int = 30, skip_calc : bool = False, reservoir_size: int = 100):
        self.data_density = 0.0
        self.n_datapoints = 0
        self.contained_datapoints = 0
        self.n_segments = 0
        self.contained_segments = 0
        self.n_reads = 0

        # Additional RED attributes
        self.signal_stats = OnlineMeanVar(initlen, skip_calc)
        self.dwell_time_stats = OnlineMeanVar(initlen, skip_calc) # TODO is this allowed, AFAIK segment lengths are negative binomial distributed instead of normally distributed

        # reservoir sampler
        # self.signal_reservoir = Reservoir(reservoir_size)
        # self.dwell_time_reservoir = Reservoir(reservoir_size)

    def append(self, xs: np.ndarray):
        """Appends new values to the RED model and updates mean/variance."""
        self.signal_stats.append(xs)
        # self.signal_reservoir.add(xs)

        dwell_time = np.array([len(xs)])
        self.dwell_time_stats.append(dwell_time)
        # self.dwell_time_reservoir.add(dwell_time)

    def get_signal_mean_stdev(self) -> tuple[float, float]:
        """Computes and retrieves the signal mean and standard deviation safely."""
        return self.signal_stats.get_mean_stdev()
    
    def get_dwell_time_mean_stdev(self) -> tuple[float, float]:
        """Computes and retrieves the dwell time mean and standard deviation safely."""
        return self.dwell_time_stats.get_mean_stdev()

    def add_reads(self, n: int = 1):
        """Increments the read count."""
        self.n_reads += n

    def add_segments(self, n: int = 1):
        """Increments the segment count."""
        self.n_segments += n

    def add_datapoints(self, n: int):
        """Increments the total datapoints count."""
        self.n_datapoints += n

    def add_contained_datapoints(self, n: int):
        """Increments the contained datapoints count."""
        self.contained_datapoints += n

    def add_contained_segments(self, n: int):
        """Increments the contained segments count."""
        self.contained_segments += n

    def add_data_density(self, density: float):
        """Increments the data density value."""
        self.data_density += density

    def set_data_density(self, density: float):
        """Safely sets the data density value."""
        self.data_density = density

    def expected_model_density(self) -> float:
        """Computes the expected model density based on standard deviation."""
        _, stdev = self.signal_stats.get_mean_stdev()
        return 1 / (2 * np.sqrt(np.pi) * stdev) if stdev else np.nan  # Avoid division by zero

    def __str__(self):
        """String representation for debugging."""
        signal_mean, signal_std = self.signal_stats.get_mean_stdev()
        segment_len_mean, segment_len_std = self.dwell_time_stats.get_mean_stdev() # TODO add to output
        return (
            f"{signal_mean:.8f}\t{signal_std:.8f}\t"
            f"{segment_len_mean:.8f}\t{segment_len_std:.8f}\t"
            f"{self.data_density:.8f}\t{self.expected_model_density():.8f}\t"
            f"{self.n_datapoints:.0f}\t{self.contained_datapoints:.0f}\t"
            f"{self.n_segments:.0f}\t{self.contained_segments:.0f}\t"
            f"{self.n_reads:.0f}"
        )

    def magnipore_string(self) -> str:
        """String representation for Magnipore output."""
        return (
            f"{self.n_datapoints:.0f}\t{self.contained_datapoints:.0f}\t"
            f"{self.n_segments:.0f}\t{self.contained_segments:.0f}\t"
            f"{self.n_reads:.0f}"
        )

    # def get_samples(self) -> np.ndarray:
    #     """
    #     Retrieves the current samples from the reservoir.

    #     Returns
    #     -------
    #     np.ndarray
    #         An array containing the samples currently held in the reservoir.
    #     """
    #     return self.signal_reservoir.samples()
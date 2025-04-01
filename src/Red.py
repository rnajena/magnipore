import numpy as np
from src.Reservoir import Reservoir

class Red:
    """
    Represents read event distribution attributes, using an online mean-variance tracker.
    """

    def __init__(self, initlen: int = 30, skip_calc : bool = False, reservoir_size : int = 100):
        # Online mean-variance tracking
        self.k = 0      # Shift value (computed from first batch)
        self.n = 0        # Total count of values
        self.ex = 0.0     # Sum of (X - k)
        self.ex2 = 0.0    # Sum of (X - k)^2
        self.initvec = np.zeros(0, dtype=np.float32)         # Shared buffer for initial values
        self.INITLEN = initlen                # Minimum number of values before shift computation
        self.SKIP_CALC = skip_calc

        # Additional RED attributes
        self.data_density = 0.0
        self.n_datapoints = 0
        self.contained_datapoints = 0
        self.n_segments = 0
        self.contained_segments = 0
        self.n_reads = 0
        self.mean = None
        self.std = None
        self.var = None

        # reservoir sampler
        self.reservoir = Reservoir(100)

    def append(self, xs: np.ndarray):
        """Appends new values to the RED model and updates mean/variance."""
        if not self.initvec.size:
            self.initvec = xs
        else:
            self.initvec = np.append(self.initvec, xs) # More efficient than `np.append()`

        # Process data if buffer is full
        if len(self.initvec) >= self.INITLEN:
            self._flush()

        self.reservoir.add(xs)

    def _flush(self):
        """Processes buffered values and updates running statistics safely."""
        if not self.n:  # First flush: Compute k
            self.k = np.mean(self.initvec)

        self.n += self.initvec.size
        diff = self.initvec - self.k
        self.ex += np.sum(diff)
        # self.ex2 += np.sum(diff * diff)
        self.ex2 += np.einsum('i,i->', diff, diff)  # Equivalent to np.sum(diff**2), but faster
        self.initvec = np.zeros(0, dtype=np.float32)  # Reset buffer

    def get_mean_stdev(self) -> tuple[float, float]:
        """Computes and retrieves the mean and standard deviation safely."""
        if self.SKIP_CALC:
            return self.mean, self.std
        
        self._flush()
        if self.n < 2:
            return self.k, 0.0  # Avoid division by zero
        
        self.var = (self.ex2 - (self.ex ** 2) / self.n) / (self.n - 1)
        self.mean = (self.ex / self.n) + self.k
        self.std = np.sqrt(self.var)
        return self.mean, self.std

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
        _, stdev = self.get_mean_stdev()
        return 1 / (2 * np.sqrt(np.pi) * stdev) if stdev else np.nan  # Avoid division by zero

    def __str__(self):
        """String representation for debugging."""
        self.mean, self.std = self.get_mean_stdev()
        return (
            f"{self.mean:.8f}\t{self.std:.8f}\t"
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

    def get_samples(self) -> np.ndarray:
        """
        Retrieves the current samples from the reservoir.

        Returns
        -------
        np.ndarray
            An array containing the samples currently held in the reservoir.
        """
        return self.reservoir.samples()
from numpy import float32, mean, sum, dot, zeros, append
from math import sqrt

class OnlineMeanVar:
    __slots__ = (
        "k", "n", "ex", "ex2", "initvec", "INITLEN", "SKIP_CALC", "mean", "std"
    )
  
    def __init__(self, initlen: int = 20, skip_calc : bool = False):
        # Online mean-variance tracking
        self.k : float = 0                # Shift value (computed from first batch)
        self.n : int = 0                  # Total count of values
        self.ex : float = 0.0             # Sum of (X - k)
        self.ex2 : float = 0.0            # Sum of (X - k)^2
        self.initvec = zeros(0, float32)
        self.INITLEN : int = initlen      # Minimum number of values before shift computation
        self.SKIP_CALC : bool = skip_calc
        self.mean : float | None = None
        self.std : float | None = None

    def append(self, xs):
        """Appends new values to the RED model and updates mean/variance."""
        if not self.initvec.size:
            self.initvec = xs
        else:
            self.initvec = append(self.initvec, xs)

        # Process data if buffer is full
        if len(self.initvec) >= self.INITLEN or self.n > 0:
            self._flush()

    def _flush(self):
        """Processes buffered values and updates running statistics safely."""
        if not self.n:  # First flush: Compute k
            self.k = float(mean(self.initvec))

        self.n += self.initvec.size
        diff = self.initvec - self.k
        self.ex += sum(diff) # float(sum(diff, dtype=float64))
        self.ex2 += dot(diff, diff)  # dot avoids temporary array
        self.initvec = zeros(0, float32)  # Reset buffer

    def get_mean_stdev(self) -> tuple[float, float]:
        """Computes and retrieves the mean and standard deviation safely."""
        if self.SKIP_CALC:
            return self.mean, self.std
        
        self._flush()
        if self.n < 2:
            return self.k, 0.0  # Avoid division by zero
        
        var = (self.ex2 - (self.ex ** 2) / self.n) / (self.n - 1)
        self.mean = (self.ex / self.n) + self.k
        self.std = sqrt(var)
        return self.mean, self.std
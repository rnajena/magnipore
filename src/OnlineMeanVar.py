import numpy as np

class OnlineMeanVar:
  
    def __init__(self, initlen: int = 30, skip_calc : bool = False):
        # Online mean-variance tracking
        self.k = 0      # Shift value (computed from first batch)
        self.n = 0        # Total count of values
        self.ex = 0.0     # Sum of (X - k)
        self.ex2 = 0.0    # Sum of (X - k)^2
        self.initvec = np.zeros(0, dtype=np.float32)         # Shared buffer for initial values
        self.INITLEN = initlen                # Minimum number of values before shift computation
        self.SKIP_CALC = skip_calc
        self.mean = None
        self.std = None

    def append(self, xs: np.ndarray):
        """Appends new values to the RED model and updates mean/variance."""
        if not self.initvec.size:
            self.initvec = xs
        else:
            self.initvec = np.append(self.initvec, xs) # More efficient than `np.append()`

        # Process data if buffer is full
        if len(self.initvec) >= self.INITLEN:
            self._flush()

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
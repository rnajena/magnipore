# from numpy import sqrt, pi, nan
from numpy import mean, std
from src.Reservoir import Reservoir
# from fitter import Fitter
# from src.OnlineMeanVar import OnlineMeanVar
import pickle

class Red:
    """
    Represents read event distribution attributes, using an online mean-variance tracker.
    """
    __slots__ = [
        'data_density', 'n_datapoints', 'contained_datapoints',
        'n_segments', 'contained_segments', 'n_reads',
        'signal_reservoir' # , 'dwell_time_reservoir'
    ]

    # TODO try k = 1000 (should reduce variance by around root(10) compared to k=100)
    # helpful for multiple modes - better separation of modes
    def __init__(self, initlen: int | tuple = 1000, skip_calc : bool = False): #, reservoir_size: int = 100):
        self.data_density = 0.0
        self.n_datapoints = 0
        self.contained_datapoints = 0
        self.n_segments = 0
        self.contained_segments = 0
        self.n_reads = 0

        # Additional RED attributes
        # self.signal_stats = OnlineMeanVar(initlen, skip_calc)
        # self.dwell_time_stats = OnlineMeanVar(initlen, skip_calc) # TODO is this allowed?, AFAIK segment lengths are negative binomial distributed instead of normally distributed
        #! This already leads to memory issues with E.coli sized genome!

        # reservoir sampler
        self.signal_reservoir = Reservoir(initlen, dtype="float16") # print(np.float16().dtype.num)  # Output: 23
        # self.signal_reservoir = Reservoir(initlen, dtype='float32') # maybe use float16? to reduce memory usage
        # self.dwell_time_reservoir = Reservoir(initlen, dtype="uint16") # print(np.uint16().dtype.num)  # Output: 4

    def append(self, xs):
        """Appends new values to the RED model and updates mean/variance."""
        # self.signal_stats.append(xs)
        # self.dwell_time_stats.append(xs.shape)
        try:
            # max length possible is 65535 (which is already way too large for a ONT segment - either a stuck pore or missegmentation), above will be skipped
            # self.dwell_time_reservoir.add(xs.shape)
            self.signal_reservoir.add(xs)
        except OverflowError:
            pass # currently just do nothing if overflow occurs

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

    # def get_signal_mean_stdev(self) -> tuple[float, float]:
    #     """Computes and retrieves the signal mean and standard deviation safely."""
    #     return self.signal_stats.get_mean_stdev()

    def get_signal_mean_stdev(self) -> tuple[float, float]:
        """Computes and retrieves the signal mean and standard deviation safely."""
        samples = self.signal_reservoir.samples()
        return mean(samples), std(samples)
    
    # def get_dwell_time_mean_stdev(self) -> tuple[float, float]:
    #     """Computes and retrieves the dwell time mean and standard deviation safely."""
    #     return self.dwell_time_stats.get_mean_stdev()

    # def expected_model_density(self) -> float:
    #     """Computes the expected model density based on standard deviation."""
    #     _, stdev = self.get_signal_mean_stdev()
    #     return 1 / (2 * sqrt(pi) * stdev) if stdev else nan  # Avoid division by zero

    # def __str__(self):
    #     """String representation for debugging."""
    #     signal_mean, signal_std = self.get_signal_mean_stdev()
    #     # segment_len_mean, segment_len_std = self.dwell_time_stats.get_mean_stdev()
    #     return (
    #         f"{signal_mean:.8f}\t{signal_std:.8f}\t"
    #         # f"{segment_len_mean:.8f}\t{segment_len_std:.8f}\t"
    #         f"{self.data_density:.8f}\t{self.expected_model_density():.8f}\t"
    #         f"{self.n_datapoints:.0f}\t{self.contained_datapoints:.0f}\t"
    #         f"{self.n_segments:.0f}\t{self.contained_segments:.0f}\t"
    #         f"{self.n_reads:.0f}"
    #     )

    # def magnipore_string(self) -> str:
    #     """String representation for Magnipore output."""
    #     return (
    #         f"{self.n_datapoints:.0f}\t{self.contained_datapoints:.0f}\t"
    #         f"{self.n_segments:.0f}\t{self.contained_segments:.0f}\t"
    #         f"{self.n_reads:.0f}"
    #     )

    def get_signals(self):
        """
        Retrieves the current samples from the reservoir.

        Returns
        -------
        np.ndarray
            An array containing the samples currently held in the reservoir.
        """
        return self.signal_reservoir.samples()
    
    # def get_dwell_times(self):
    #     """
    #     Retrieves the current samples from the reservoir.

    #     Returns
    #     -------
    #     np.ndarray
    #         An array containing the samples currently held in the reservoir.
    #     """
    #     return self.dwell_time_reservoir.samples()
    
    # def fit(self):
    #     """
    #     Fits a distribution to the samples in the reservoir.

    #     Returns
    #     -------
    #     dict
    #         A dictionary containing the fitted distribution parameters.
    #     """
    #     f = Fitter(self.get_samples(), distributions=['norm', 'lognorm', 'gamma'])
    #     f.fit()
    #     return f.fitted_param['norm']
    
    # '_fit', 'alpha', 'anglit', 'arcsine', 'argus', 'beta', 'betaprime', 'bradford', 'burr', 'burr12', 'cauchy', 'chi', 'chi2', 'cosine', 'crystalball', 'dgamma', 'dpareto_lognorm', 'dweibull', 'erlang', 'expon', 'exponnorm', 'exponpow', 'exponweib', 'f', 'fatiguelife', 'fisk', 'foldcauchy', 'foldnorm', 'gamma', 'gausshyper', 'genexpon', 'genextreme', 'gengamma', 'genhalflogistic', 'genhyperbolic', 'geninvgauss', 'genlogistic', 'gennorm', 'genpareto', 'gibrat', 'gompertz', 'gumbel_l', 'gumbel_r', 'halfcauchy', 'halfgennorm', 'halflogistic', 'halfnorm', 'hypsecant', 'invgamma', 'invgauss', 'invweibull', 'irwinhall', 'jf_skew_t', 'johnsonsb', 'johnsonsu', 'kappa3', 'kappa4', 'ksone', 'kstwo', 'kstwobign', 'landau', 'laplace', 'laplace_asymmetric', 'levy', 'levy_l', 'levy_stable', 'loggamma', 'logistic', 'loglaplace', 'lognorm', 'loguniform', 'lomax', 'maxwell', 'mielke', 'moyal', 'multivariate_normal', 'nakagami', 'ncf', 'nct', 'ncx2', 'norm', 'norminvgauss', 'pareto', 'pearson3', 'powerlaw', 'powerlognorm', 'powernorm', 'rayleigh', 'rdist', 'recipinvgauss', 'reciprocal', 'rel_breitwigner', 'rice', 'rv_continuous', 'rv_histogram', 'semicircular', 'skewcauchy', 'skewnorm', 'studentized_range', 't', 'trapezoid', 'trapz', 'triang', 'truncexpon', 'truncnorm', 'truncpareto', 'truncweibull_min', 'tukeylambda', 'uniform', 'vonmises', 'vonmises_fisher', 'vonmises_line', 'wald', 'weibull_max', 'weibull_min', 'wrapcauchy'


    def __getstate__(self):
        """
        Custom method for pickling the Red object.
        Returns a dictionary of the object's state.
        """
        state = {slot: getattr(self, slot) for slot in self.__slots__}
        # Ensure Reservoir objects are picklable
        state['signal_reservoir'] = pickle.dumps(self.signal_reservoir)
        # state['dwell_time_reservoir'] = pickle.dumps(self.dwell_time_reservoir)
        return state

    def __setstate__(self, state):
        """
        Custom method for unpickling the Red object.
        Restores the object's state from the dictionary.
        """
        # Restore Reservoir objects from their pickled state
        state['signal_reservoir'] = pickle.loads(state['signal_reservoir'])
        # state['dwell_time_reservoir'] = pickle.loads(state['dwell_time_reservoir'])
        for slot, value in state.items():
            setattr(self, slot, value)
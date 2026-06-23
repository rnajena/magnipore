# adapted from https://github.com/rnajena/ONT-Nuc-Stats/blob/master/python-online-algorithms/reservoir.py
import numpy as np

# BUG Currently entries are initialized to np.zeros. This is a bug for many applications. The inner array should be dynamically sized then.

class Reservoir:
  """
  Implements reservoir sampling using algorithm L.

  Parameters:
  k (int): Number of samples to hold.
  """
  __slots__ = ['rs', 'k', 'd', 'W', 'cnt', 'next', 'dtype']

  def __init__(self, k : int | tuple, dtype : str = "float32"):
    self.rs = None  # Delay initialization of the reservoir array
    self.k : int | tuple = k
    self.d = k[0] if type(self.k) is tuple else k
    self.W : float = self._genW()
    self.cnt : int = 0
    self.next : int = 0
    self.dtype : str = dtype  # Store dtype for later initialization

  def setK(self, k : int | tuple):
    """
    Sets the size of the reservoir to k. If the reservoir is already full,
    it will not change the size but will reset the internal state.
    
    Parameters:
    k (int | tuple): New size of the reservoir.
    """
    self.k = k

  def add(self, xs : np.ndarray):
    """
    Adds elements from the iterable `xs` to the reservoir. If the reservoir
    is not yet full, it adds elements directly. Once full, it replaces elements
    probabilistically according to reservoir sampling algorithm L. Updates
    internal counters and calculates the next target index for possible replacement.
    
    Parameters:
    xs (iterable): An iterable of elements to be added to the reservoir.
    """
    if self.rs is None:
      # Initialize the reservoir array when data is first added to reduce memory usage
      self.rs = np.zeros(self.k, self.dtype)
      
    for x in xs:
      self.cnt += 1
      # always add elements if array not filled yet
      if (self.cnt < self.d):
        self.rs[self.cnt-1] = x
      # always add, but prepare the next jump target
      elif (self.cnt == self.d):
        self.rs[self.cnt-1] = x
        self.next = self.cnt + self._nextjump()
      # reached a jump target
      # PERF It would be possible to implement this even more efficient, by not looping over x in xs.
      elif (self.cnt >= self.next):
        self.rs[int(np.random.random()*self.d)] = x
        self.W *= self._genW()
        self.next = self.cnt + self._nextjump()

  def samples(self) -> np.ndarray:
    """
    Returns the current samples in the reservoir.

    If the reservoir is full, it will return k samples, otherwise it will return
    the number of samples that have been added so far.

    Returns:
      np.ndarray: An array containing the current samples in the reservoir.
    """
    if self.rs is None:
      return np.array([], dtype=self.dtype)  # Return an empty array if no data has been added
    return self.rs[:min(self.d, self.cnt)]
    
  def _nextjump(self) -> int:
    """
    Calculates the next jump length L in the reservoir sampling algorithm.

    Returns an integer value L, which is the number of elements to skip
    in the input stream before the next element is added to the reservoir.
    This value is calculated as floor(log(random())/log(1-W))+1, where W is
    the current value of W in the reservoir sampling algorithm.
    """
    return np.floor(np.log(np.random.random())/np.log(1-self.W))+1
    
  def _genW(self) -> float:
    """
    Generates a random value for W in the reservoir sampling algorithm.

    Based on the value of k, this function returns a random value for W, which is
    used to calculate the next jump length L. The value of W is calculated as
    exp(log(random())/k).

    Returns:
      float: A random value for W.
    """
    return np.exp(np.log(np.random.random()) / self.d)
  
  def __getstate__(self):
    """
    Custom method for pickling the Reservoir object.
    Returns a dictionary of the object's state.
    """
    state = {slot: getattr(self, slot) for slot in self.__slots__}
    # Serialize the NumPy array as bytes if it exists
    if self.rs is not None:
      state['rs'] = self.rs.tobytes()
      state['rs_dtype'] = self.dtype    # Store dtype for reconstruction
      state['rs_shape'] = self.rs.shape # Store shape for reconstruction
    else:
      state['rs'] = None
    return state

  def __setstate__(self, state):
    """
    Custom method for unpickling the Reservoir object.
    Restores the object's state from the dictionary.
    """
    # Restore all attributes in __slots__, setting defaults if missing
    for slot in self.__slots__:
      setattr(self, slot, state.get(slot, None))

    # Handle the 'rs' attribute specifically
    if state.get('rs') is not None:
      try:
        # Attempt to reconstruct 'rs' using the expected keys
        self.rs = np.frombuffer(
          state['rs'], 
          dtype=np.dtype(state.get('rs_dtype', 'float32'))
        ).reshape(state['rs_shape'])
      except KeyError:
        # Fallback for older versions or missing keys
        self.rs = None
      except TypeError:
        # Handle cases where 'rs' is not bytes-like
        if isinstance(state['rs'], list):
          self.rs = np.array(state['rs'], dtype="float32").reshape((self.d,))
        else:
          raise
    else:
      self.rs = None
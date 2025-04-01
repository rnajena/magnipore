# adapted from https://github.com/rnajena/ONT-Nuc-Stats/blob/master/python-online-algorithms/reservoir.py
import numpy as np

# BUG Currently entries are initialized to np.zeros. This is a bug for many applications. The inner array should be dynamically sized then.

class Reservoir:
  """
  Implements reservoir sampling using algorithm L.

  Parameters:
  k (int): Number of samples to hold.
  """
  def __init__(self, k : int = 100):
    assert(k>0)
    self.rs = np.zeros(k)
    self.k = k
    self.W = self._genW()
    self.cnt = 0
    self.next = 0

  def add(self, xs : np.ndarray):
    """
    Adds elements from the iterable `xs` to the reservoir. If the reservoir
    is not yet full, it adds elements directly. Once full, it replaces elements
    probabilistically according to reservoir sampling algorithm L. Updates
    internal counters and calculates the next target index for possible replacement.
    
    Parameters:
    xs (iterable): An iterable of elements to be added to the reservoir.
    """
    for x in xs:
      self.cnt += 1
      # always add elements if array not filled yet
      if (self.cnt < self.k):
        self.rs[self.cnt-1] = x
      # always add, but prepare the next jump target
      elif (self.cnt == self.k):
        self.rs[self.cnt-1] = x
        self.next = self.cnt + self._nextjump()
      # reached a jump target
      # PERF It would be possible to implement this even more efficient, by not looping over x in xs.
      elif (self.cnt >= self.next):
        self.rs[int(np.random.random()*self.k)] = x
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
    r = min(self.k, self.cnt)
    return self.rs[:r]
    
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
    return np.exp(np.log(np.random.random()) / self.k)
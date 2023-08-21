from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

from magnipore import magnipore, nanosherlock, Helper

__all__ = [
    "magnipore",
    "nanosherlock",
    "Helper",
]
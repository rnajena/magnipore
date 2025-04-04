try:
  from src._version import version as __version__
  __version_str__ = f"v{__version__}"
except ImportError:
  __version__ = "unknown"
  __version_str__ = "unknown"
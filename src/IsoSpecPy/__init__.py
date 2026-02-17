# Try to use nanobind version first, fall back to CFFI if not available
try:
    from .IsoSpecPy_nb import *

    _backend = "nanobind"
except ImportError:
    from .IsoSpecPy import *

    _backend = "cffi"

import importlib.metadata

__version__ = importlib.metadata.version("IsoSpecPy")

# Old, deprecated name, for compatibility with 1.9.X only
IsoLayered = IsoTotalProb


# For backward compatibility with 1.0.X:
class CompatIsoWrapper(object):
    _initialized = False

    def __getattr__(self, name):
        # Lazy import to avoid circular dependencies
        if not self._initialized:
            try:
                from .IsoSpecPyOld import IsoSpec, IsoSpecify, IsoPlot

                self.IsoSpec = IsoSpec
                self.IsoSpecify = IsoSpecify
                self.IsoPlot = IsoPlot
                self._initialized = True
            except ImportError:
                # If old API is not available, that's okay for nanobind-only usage
                pass

        if not hasattr(self, name):
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{name}'"
            )
        return getattr(self, name)


IsoSpecPy = CompatIsoWrapper()

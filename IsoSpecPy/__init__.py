from .IsoSpecPy import *
import importlib.metadata

__version__ = importlib.metadata.version("IsoSpecPy")

# Old, deprecated name, for compatibility with 1.9.X only
IsoLayered = IsoTotalProb

# For backward compatibility with 1.0.X:
class CompatIsoWrapper(object):
    def __init__(self):
        from .IsoSpecPyOld import IsoSpec, IsoSpecify, IsoPlot
        self.IsoSpec = IsoSpec
        self.IsoSpecify = IsoSpecify
        self.IsoPlot = IsoPlot


IsoSpecPy = CompatIsoWrapper()


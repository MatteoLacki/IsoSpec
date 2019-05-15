from .IsoSpecPy import *


# For backward compatibility with 1.0.X:
class CompatIsoWrapper(object):
    def __init__(self):
        from .IsoSpecPyOld import IsoSpec, IsoSpecify, IsoPlot
        self.IsoSpec = IsoSpec
        self.IsoSpecify = IsoSpecify
        self.IsoPlot = IsoPlot


IsoSpecPy = CompatIsoWrapper()


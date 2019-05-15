from .IsoSpecPyNew import *


# For backward compatibility with 1.0.X:
class CompatIsoWrapper(object):
    def __init__(self):
        from .IsoSpecPy import IsoSpec, IsoSpecify, IsoPlot
        self.IsoSpec = IsoSpec
        self.IsoSpecify = IsoSpecify
        self.IsoPlot = IsoPlot


IsoSpecPy = CompatIsoWrapper()


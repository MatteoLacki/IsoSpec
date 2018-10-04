# Calculates the isotopic distribution of water in several ways

import IsoSpecPy
from math import exp

try:
    if IsoSpecPy.__version__[:3] != '1.9':
        raise AttributeError
except AttributeError:
    print "This file is meant to be used with IsoSpecPy version 1.9.X. You seem to have a different version installed on your system."
    import sys
    sys.exit(-1)

i = IsoSpecPy.IsoOrderedGenerator(formula="H2O1", get_confs=True)

print "Calculating isotopic distribution of water. Here's a list of all configurations, in a guaranteed order of nonincreasing probability:"

for mass, log_prob, conf in i:
    print "Mass:", mass
    print "log(probability):", log_prob
    print "probability:", exp(log_prob)
    print "Number of Protium atoms:", conf[0][0]
    print "Number of Deuterium atoms", conf[0][1]
    print "Number of O16 atoms:", conf[1][0]
    print "Number of O17 atoms:", conf[1][1]
    print "Number of O18 atoms:", conf[1][2]

print
print "Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?"

hydrogen_probs = (0.5, 0.5)
oxygen_probs = (0.5, 0.3, 0.2)
hydrogen_masses = (1.00782503207, 2.0141017778)
oxygen_masses = (15.99491461956, 16.99913170, 17.9991610)

i = IsoSpecPy.IsoOrderedGenerator(dimNumber = 2, isotopeNumbers = (2, 3), atomCounts = (2, 1), isotopeMasses = (hydrogen_masses, oxygen_masses), isotopeProbabilities = (hydrogen_probs, oxygen_probs), get_confs=True)


mass, log_prob, conf = i.__iter__().next()

print "The first configuration has the following parameters:"
print "Mass:", mass
print "log(probability):", log_prob
print "probability:", exp(log_prob)
print "Number of Protium atoms:", conf[0][0]
print "Number of Deuterium atoms", conf[0][1]
print "Number of O16 atoms:", conf[1][0]
print "Number of O17 atoms:", conf[1][1]
print "Number of O18 atoms:", conf[1][2]






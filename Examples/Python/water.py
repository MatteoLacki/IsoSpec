# Calculates the isotopic distribution of water in several ways

import IsoSpecPy
from math import exp

try:
    if IsoSpecPy.__version__[:3] != '1.9':
        raise AttributeError
except AttributeError:
    print("This file is meant to be used with IsoSpecPy version 1.9.X. You seem to have a different version installed on your system.")
    import sys
    sys.exit(-1)

i = IsoSpecPy.IsoLayeredGenerator(formula="H2O1", prob_to_cover = 0.999, get_confs=True)

print("Calculating isotopic distribution of water. Here's a list of configurations necessary to cover at least 0.999 of total probability:")

for mass, log_prob, conf in i:
    print("")
    print("Mass: " + str(mass))
    print("log(probability): " + str(log_prob))
    print("probability: " + str(exp(log_prob)))
    print("Number of Protium atoms: " + str(conf[0][0]))
    print("Number of Deuterium atoms: " + str(conf[0][1]))
    print("Number of O16 atoms: " + str(conf[1][0]))
    print("Number of O17 atoms: " + str(conf[1][1]))
    print("Number of O18 atoms: " + str(conf[1][2]))

print("")
print("Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?")

hydrogen_probs = (0.5, 0.5)
oxygen_probs = (0.5, 0.3, 0.2)
hydrogen_masses = (1.00782503207, 2.0141017778)
oxygen_masses = (15.99491461956, 16.99913170, 17.9991610)

i = IsoSpecPy.IsoLayeredGenerator(atomCounts = (2, 1), isotopeMasses = (hydrogen_masses, oxygen_masses), isotopeProbabilities = (hydrogen_probs, oxygen_probs), prob_to_cover = 0.999, get_confs=True)


mass, log_prob, conf = next(i.__iter__())

print("The first configuration has the following parameters:")
print("Mass: " + str(mass))
print("log(probability): " + str(log_prob))
print("probability: " + str(exp(log_prob)))
print("Number of Protium atoms: " + str(conf[0][0]))
print("Number of Deuterium atoms: " + str(conf[0][1]))
print("Number of O16 atoms: " + str(conf[1][0]))
print("Number of O17 atoms: " + str(conf[1][1]))
print("Number of O18 atoms: " + str(conf[1][2]))







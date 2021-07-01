# Calculates the isotopic distribution of water in several ways

'''
We shall walk through a set of configurations which covers at least 99.9% of the total
probability. For water we could obviously go through the entire spectrum (100%), but for
larger molecules the entire spectrum has far too many configurations. Luckily, most of the
likelihood is concentrated in a relatively small set of most probable isotopologues
- and this method allows one to quickly calculate such a set, parametrising on the
percentage of coverage of total probability space required.

This is usually better than just calculating all isotopes with probabilities above a
iven threshold, as it allows one to directly parametrise on the accuracy of the
simplified spectrum - that is, the L1 distance between the full and simplified spectrum
will be less than (in this example) 0.001.

If for some reason one would wish to just calculate a set of peaks with probabilities
above a given threshold - it is possible using the IsoThresholdGenerator class.

Note: the returned set will usually contain a bit more configurations than necessary
to achieve the desired coverage. These configurations need to be computed anyway, however
it is possible to discard them using the optional trim argument.
'''

import IsoSpecPy
from math import exp

try:
    if IsoSpecPy.__version__[:2] != '2.':
        raise AttributeError
except AttributeError:
    print("This file is meant to be used with IsoSpecPy version 2.0.X. You seem to have a different version installed on your system.")
    import sys
    sys.exit(-1)

i = IsoSpecPy.IsoTotalProb(formula="H2O1", prob_to_cover = 0.999, get_confs=True)

print("Calculating isotopic distribution of water. Here's a list of configurations necessary to cover at least 0.999 of total probability:")

for mass, prob, conf in i:
    print("")
    print("Mass: " + str(mass))
    print("probability: " + str(prob))
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

i = IsoSpecPy.IsoTotalProb(atomCounts = (2, 1), isotopeMasses = (hydrogen_masses, oxygen_masses), isotopeProbabilities = (hydrogen_probs, oxygen_probs), prob_to_cover = 0.999, get_confs=True)



print("The first configuration has the following parameters:")
print("Mass: " + str(i.masses[0]))
print("probability: " + str(i.probs[0]))
conf = i.confs[0]
conf_hydrogen = conf[0]
conf_oxygen = conf[1]
print("Number of Protium atoms: " + str(conf_hydrogen[0]))
print("Number of Deuterium atoms: " + str(conf_hydrogen[1]))
print("Number of O16 atoms: " + str(conf_oxygen[0]))
print("Number of O17 atoms: " + str(conf_oxygen[1]))
print("Number of O18 atoms: " + str(conf_oxygen[2]))







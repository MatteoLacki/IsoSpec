# Calculates the isotopic distribution of water in several ways

from IsoSpecPy import IsoSpecPy
from math import exp

i = IsoSpecPy.IsoSpec.IsoFromFormula("H2O1", 0.9)

print "The isotopologue set containing at least 0.9 probability has", len(i), "element(s)"

confs = i.getConfs()

print "The first configuration has the following parameters:"
print "Mass:", confs[0][0]
print "log-prob:", confs[0][1] 
print "probability:", exp(confs[0][1])
print "Protium atoms:", confs[0][2][0][0]
print "Deuterium atoms", confs[0][2][0][1]
print "O16 atoms:", confs[0][2][1][0]
print "O17 atoms:", confs[0][2][1][1]
print "O18 atoms:", confs[0][2][1][2]

print
print "Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?"

hydrogen_probs = (0.5, 0.5)
oxygen_probs = (0.5, 0.3, 0.2)
hydrogen_masses = (1.00782503207, 2.0141017778)
oxygen_masses = (15.99491461956, 16.99913170, 17.9991610)
atom_counts = (2, 1)

i = IsoSpecPy.IsoSpec(atom_counts, (hydrogen_masses, oxygen_masses), (hydrogen_probs, oxygen_probs), 0.9)

print "The isotopologue set containing at least 0.9 probability has", len(i), "element(s)"

confs = i.getConfs()

print "The first configuration has the following parameters:"
print "Mass:", confs[0][0]
print "log-prob:", confs[0][1]
print "probability:", exp(confs[0][1])
print "Protium atoms:", confs[0][2][0][0]
print "Deuterium atoms", confs[0][2][0][1]
print "O16 atoms:", confs[0][2][1][0]
print "O17 atoms:", confs[0][2][1][1]
print "O18 atoms:", confs[0][2][1][2]





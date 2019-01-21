# Calculates the isotopic distribution of isotopically labelled glucose, where one of the carbons is replaced with 14C
# We assume that 95% of the molecules are radiolabelled.

import IsoSpecPy
from IsoSpecPy import PeriodicTbl
from math import exp

try:
    if IsoSpecPy.__version__[:3] != '1.9':
        raise AttributeError
except AttributeError:
    print("This file is meant to be used with IsoSpecPy version 1.9.X. You seem to have a different version installed on your system.")
    import sys
    sys.exit(-1)



# As one of the atoms has a nonstandard isotopic distribution we can't do the construction from chemical formula
# We have to manually grab isotopic masses and probabilities instead, the construct an artificial "element" which is carbon
# with one added isotope (14C) and shifted isotopic distribution

# Formula of radiolabeled glucose is C5H12O6(14C)1
# First, deal with standard elements
normal_carbon_masses = PeriodicTbl.symbol_to_masses["C"]
normal_carbon_probs = PeriodicTbl.symbol_to_probs["C"]
hydrogen_masses = PeriodicTbl.symbol_to_masses["H"]
hydrogen_probs = PeriodicTbl.symbol_to_probs["H"]
oxygen_masses = PeriodicTbl.symbol_to_masses["O"]
oxygen_probs = PeriodicTbl.symbol_to_probs["O"]

# Now, the radiolabelled carbon
# 14C isn't normally considered, here we add an extra isotope to the standard ones
radiolabelled_carbon_masses = PeriodicTbl.symbol_to_masses["C"] + (14.003241989,)

# Assuming that the labelling was only 95% efficient, that is only 95% 
# of the molecules have standard C replaced with 14C. Non-replaced molecules have standard
# isotopic abundance (realtive to each other)
radiolabelled_carbon_probs = (0.05*normal_carbon_probs[0], 0.05*normal_carbon_probs[1], 0.95)

atom_counts = (5, 12, 6, 1) # 5 normal carbons, 12 H's, 6 O's and one radiolabelled carbon

i = IsoSpecPy.IsoLayeredGenerator(atomCounts = atom_counts, 
                                  isotopeMasses = (normal_carbon_masses, hydrogen_masses, oxygen_masses, radiolabelled_carbon_masses), 
                                  isotopeProbabilities = (normal_carbon_probs, hydrogen_probs, oxygen_probs, radiolabelled_carbon_probs), 
                                  prob_to_cover = 0.99, 
                                  get_confs=True)


print("The list of configurations that, taken together, cover at least 99% of the probability space is:")

for mass, prob, conf in i:
    print("=====================================================")

    print("Mass: " + str(mass))
    print("probability: " + str(prob))
    
    print("Number of 12C atoms: " + str(conf[0][0] + conf[3][0])) # Counting the normal and unsuccesfully radiolabelled atoms
    print("Number of 13C atoms: " + str(conf[0][1] + conf[3][1]))
    print("Number of 14C atoms: " + str(conf[3][2])) # Note: conf[0][2] will raise IndexError, normal carbons don't consider 14C
    print("Number of Protium atoms: " + str(conf[1][0]))
    print("Number of Deuterium atoms: " + str(conf[1][1]))
    print("Number of O16 atoms: " + str(conf[2][0]))
    print("Number of O17 atoms: " + str(conf[2][1]))
    print("Number of O18 atoms: " + str(conf[2][2]))







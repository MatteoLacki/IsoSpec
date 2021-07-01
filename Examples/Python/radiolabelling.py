# Calculates the isotopic distribution of isotopically labelled glucose, where two of the carbons are replaced with 14C
# We assume that radiolabelling is 95% efficient (that is, there is a 95% chance for each of the radiolabel atoms to be 14C, and
# a 5% chance of them being either 12C or 13C, with the probability of each of those proportional to their standard isotopic abundances)

import IsoSpecPy
from IsoSpecPy import PeriodicTbl
from math import exp

try:
    if IsoSpecPy.__version__[:2] != '2.':
        raise AttributeError
except AttributeError:
    print("This file is meant to be used with IsoSpecPy version 2.0.X. You seem to have a different version installed on your system.")
    import sys
    sys.exit(-1)



# Now, the radiolabelled carbon
# 14C isn't normally considered in the isotopic distribution, here we add an extra isotope to the standard ones
radiolabelled_carbon_masses = PeriodicTbl.symbol_to_masses["C"] + (14.003241989,)

# Assuming that the labelling was only 95% efficient, that is only 95% 
# of the radiolabel atoms have standard C replaced with 14C. Non-replaced atoms have standard
# isotopic abundance (realtive to each other)
normal_carbon_probs = PeriodicTbl.symbol_to_probs["C"]
radiolabelled_carbon_probs = (0.05*normal_carbon_probs[0], 0.05*normal_carbon_probs[1], 0.95)

i = IsoSpecPy.IsoTotalProb(formula = "C4H12O6", # The formula for glucose, sans the radiolabel atoms
                                  # Here we specify additional "elements" which occur *in addition* to those from the formula
                                  atomCounts = (2,),
                                  isotopeMasses = (radiolabelled_carbon_masses,), 
                                  isotopeProbabilities = (radiolabelled_carbon_probs,), 
                                  # And the rest of parameters for configuration
                                  prob_to_cover = 0.99,
                                  get_confs=True)

# Radiolabelling (or isotopic labelling) with more than one element looks like this: 
# Let's say we wanted to have glucose with one 14C carbon, and two deuteriums, all with 95% probability
# Then it would be:
#i = IsoSpecPy.IsoLayeredGenerator(formula = "C5H10O6", # The formula for glucose, sans the radiolabel atoms
#                                  atomCounts = (1, 2),   
#                                  isotopeMasses = (radiolabelled_carbon_masses, PeriodicTbl.symbol_to_masses["H"]),
#                                  isotopeProbabilities = (radiolabelled_carbon_probs, (0.05, 0.95)),
#                                  # And the rest of parameters for configuration
#                                  prob_to_cover = 0.99,
#                                  get_confs=True)



print("The list of configurations that, taken together, cover at least 99% of the probability space is:")

for mass, prob, conf in i:
    print("=====================================================")

    print("Mass: " + str(mass))
    print("probability: " + str(prob))
    
    # The configuration fragments corresponding to the elements from the formula is first, the maually specified elements follow in order
    print("Number of 12C atoms: " + str(conf[0][0] + conf[3][0])) # Counting the normal and unsuccesfully radiolabelled atoms
    print("Number of 13C atoms: " + str(conf[0][1] + conf[3][1]))
    print("Number of 14C atoms: " + str(conf[3][2])) # Note: conf[0][2] will raise IndexError, normal carbons don't consider 14C
    print("Number of Protium atoms: " + str(conf[1][0]))
    print("Number of Deuterium atoms: " + str(conf[1][1]))
    print("Number of O16 atoms: " + str(conf[2][0]))
    print("Number of O17 atoms: " + str(conf[2][1]))
    print("Number of O18 atoms: " + str(conf[2][2]))







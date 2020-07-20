from IsoSpecPy.IsoSpecPy import IsoTotalProb, IsoThreshold
from IsoSpecPy.PeriodicTbl import symbol_to_masses


# Suppose we want to simulate the spectrum of
# C100H200Ag46
# except that the isotope counts of silver in each molecule
# are fixed: each molecule contains exactly 44 atoms of 107Ag,
# and 2 atoms of 109Ag. The carbon and hydrogen follow their
# natural distributions.

# Here's how to do that:

# We start by defining two artificial elements, each representing
# one isotope of silver. The element is going to be monoisotopic,
# with proper mass retrieved from the periodic table.
# The new elements will be passed to Iso* methods using the
# atomCounts, isotopeMasses, isotopeProbabilities parameters,
# and the rest of the formula - using standard "formula"
# argument.

# defining extra elements
# 44 107Ag atoms and 2 109Ag
atomCounts = [44, 2]

# Retrieving the isotope masses from IsoSpec's builtin periodic table
# (one could input them manually too)
M107Ag, M109Ag = symbol_to_masses['Ag']

# The masses are passed as list of masses for each new element:
# we have two elements, each has one isotope, hence the following
# list structure
isotopeMasses = [[M107Ag], [M109Ag]]

# The single isotope of each new element occurs with 1.0 probability
isotopeProbabilities = [[1.0], [1.0]]


# peaks above .001 threshold
print("Here is the list of peaks above 0.001 intensity threshold for C100H200[107Ag]44[109Ag]2")
MZI_thr = IsoThreshold(.001,
                       'C100H200', # the remaining part of the formula
                       get_confs=True,
                       atomCounts=atomCounts,
                       isotopeMasses=isotopeMasses,
                       isotopeProbabilities=isotopeProbabilities)

print("Masses list:", list(MZI_thr.masses))
print("Probabilities list:", list(MZI_thr.probs))
print("Configurations list:", list(MZI_thr.confs))

print()
print()

# Alternatively: peaks that make 99.9% of the isotopic distribution
print("Here is the list of peaks that comprise at least 99.9% of total intensity for C100H200[107Ag]44[109Ag]2")
MZI_opt = IsoTotalProb(.999,
                       'C100H200', # the remaining part of the formula
                       get_confs=True,
                       atomCounts=atomCounts,
                       isotopeMasses=isotopeMasses,
                       isotopeProbabilities=isotopeProbabilities)

print("Masses list:", list(MZI_opt.masses))
print("Probabilities list:", list(MZI_opt.probs))
print("Configurations list:", list(MZI_opt.confs))

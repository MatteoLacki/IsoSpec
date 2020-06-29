from IsoSpecPy.IsoSpecPy import IsoParamsFromDict, IsoThreshold, IsoTotalProb, Iso
from IsoSpecPy.PeriodicTbl import symbol_to_masses

# defining extra elements
# 44 107Ag atoms and 2 109Ag
atomCounts = [44, 2]

# exact IUPAC masses of 107Ag and 109Ag
M107Ag, M109Ag = symbol_to_masses['Ag']

# isotope masses must be passed in as lists of lists of masses
# here each new element has only one isotope, so each inner list
# has only one element
isotopeMasses = [[M107Ag] , [M109Ag]]

# each new element has no mass indeterminacy
isotopeProbabilities = [[1], [1]]


# peaks above .001 threshold
MZI_thr = IsoThreshold(.001, 
                       'C100H200', # the remaining part of the formula
                       get_confs=True,
                       atomCounts=atomCounts,
                       isotopeMasses=isotopeMasses,
                       isotopeProbabilities=isotopeProbabilities)

list(MZI_thr.masses)
list(MZI_thr.probs)
list(MZI_thr.confs)

# peaks that make 99.9% of the isotopic distribution
MZI_opt = IsoTotalProb(.999,
                       'C100H200', # the remaining part of the formula
                       get_confs=True,
                       atomCounts=atomCounts,
                       isotopeMasses=isotopeMasses,
                       isotopeProbabilities=isotopeProbabilities)

list(MZI_opt.masses)
list(MZI_opt.probs)
list(MZI_opt.confs)




# Fasta with diffs.
MZI_thr = IsoThreshold(.001, 'C100H200', 
                       get_confs=True,
                       atomCounts=[10],
                       isotopeMasses=[[10,12]],
                       isotopeProbabilities=[[.8,.2]])

%load_ext autoreload
%autoreload 2
from IsoSpecPy.IsoSpecPy import IsoTotalProb

# isotopic distribution of protein with fasta sequence AAAPPQAAC
MZI_opt = IsoTotalProb(.999,
                       get_confs=True,
                       fasta='AAAPPQAAC')

# isotopic distribution of protein with fasta sequence AAAPPQAAC,
# with carbamidomethylation of cystein
carbamidomethylation = 'C2H3N1O1'
MZI_opt2 = IsoTotalProb(.999,
                        formula=carbamidomethylation,
                        get_confs=True,
                        fasta='AAAPPQAAC')
list(MZI_opt2.masses)


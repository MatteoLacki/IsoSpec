from IsoSpecPy.IsoSpecPy import IsoTotalProb

# isotopic distribution of protein with fasta sequence AAAPPGQAAC
MZI_opt = IsoTotalProb(.999,
                       fasta='AAAPPGQAAC')

# isotopic distribution of protein with fasta sequence AAAPPGQAAC,
# with thioglycinie modification
MZI_opt2 = IsoTotalProb(.999,
                        formula='O-1S1',
                        fasta='AAAPPGQAAC')
print("Masses and probabilities in AAAPPGQAAC with thioglycine modification:")
print([(m, p) for (m, p) in zip(MZI_opt2.masses, MZI_opt2.probs)])


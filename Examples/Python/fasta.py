from IsoSpecPy.IsoSpecPy import IsoTotalProb, ParsePeptideSequence

# isotopic distribution of protein with fasta sequence AAAPPGQAAC
MZI_opt = IsoTotalProb(.999,
                       fasta='AAAPPGQAAC')

# isotopic distribution of protein with fasta sequence AAAPPGQAAC,
# with an arbitrary elemental offset (thioglycine modification: -O +S).
# This is the general escape hatch for an offset that isn't a named Unimod
# entry -- you have to know the formula and get the sign right by hand.
MZI_opt2 = IsoTotalProb(.999,
                        formula='O-1S1',
                        fasta='AAAPPGQAAC')
print("Masses and probabilities in AAAPPGQAAC with thioglycine modification:")
print([(m, p) for (m, p) in zip(MZI_opt2.masses, MZI_opt2.probs)])

# For a *named* modification, [UNIMOD:<id>] notation is resolved directly
# against a table of ~980 Unimod entries in the packaged CSV, no manual
# formula arithmetic needed. peptide_sequence= is the recommended spelling
# (fasta= is now just an alias of it, kept for backward compatibility).
# Placement matches this monorepo's SAGE search-engine fork's own peptide
# output: [UNIMOD:<id>]-SEQUENCE for an N-terminal mod, X[UNIMOD:<id>] for
# an internal one, SEQUENCE-[UNIMOD:<id>] for a C-terminal one.
MZI_opt3 = IsoTotalProb(.999,
                        get_confs=True,
                        peptide_sequence='AAAPPGQAAC[UNIMOD:4]')  # Carbamidomethyl
print("Masses and probabilities in AAAPPGQAAC[UNIMOD:4] (Carbamidomethyl):")
print([(m, p) for (m, p) in zip(MZI_opt3.masses, MZI_opt3.probs)])

# ...or just the composition, no envelope:
print(dict(ParsePeptideSequence('AAAPPGQAAC[UNIMOD:4]')))


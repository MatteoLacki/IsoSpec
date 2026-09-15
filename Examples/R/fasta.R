library(IsoSpecR)

# A peptide sequence, with a named modification given as an arbitrary
# elemental offset (thioglycine: -O +S) added on top -- the general escape
# hatch for an offset that isn't a named Unimod entry. IsoSpecify only ever
# takes a raw composition vector, so combine the two by simple vector
# addition, matching Python's formula=+fasta= summation:
base <- RParsePeptideSequence("AAAPPGQAAC")
thioglycine_offset <- c(O = -1, S = 1)

composition <- base
composition["O"] <- composition["O"] + thioglycine_offset["O"]
composition["S"] <- composition["S"] + thioglycine_offset["S"]

res <- IsoSpecify(molecule = composition, stopCondition = 0.999)
cat("Masses and probabilities in AAAPPGQAAC with thioglycine modification:\n")
print(res)

# For a *named* modification, [UNIMOD:<id>] notation is resolved directly
# against a table of ~980 Unimod entries in the packaged CSV -- no
# manual formula arithmetic needed. Placement matches this monorepo's SAGE
# search-engine fork's own peptide output: "[UNIMOD:<id>]-SEQUENCE" for an
# N-terminal mod, "X[UNIMOD:<id>]" for an internal one,
# "SEQUENCE-[UNIMOD:<id>]" for a C-terminal one.
composition2 <- RParsePeptideSequence("AAAPPGQAAC[UNIMOD:4]")  # Carbamidomethyl
res2 <- IsoSpecify(molecule = composition2, stopCondition = 0.999)
cat("Masses and probabilities in AAAPPGQAAC[UNIMOD:4] (Carbamidomethyl):\n")
print(res2)

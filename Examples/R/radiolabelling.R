library(IsoSpecR)


# We will model a glucose molecule where one of the carbons has been replaced with 
# (14C) with 95% probability.


# First, we define an alternate isotopic table for elements that contains the readiocarbon: 14C with 95% probability, 
# and either 12C or 13C with their natural relative probabilities otherwise. First, we load up standard isotopic data:
data(isotopicData)
isotopes = isotopicData$IsoSpec

# And we add an artificial element X, with 3 isotopes, that will represent the radiocarbon. 
radiolabel = data.frame(
    element = c('X', 'X', 'X'),
    isotope = c('X12', 'X13', 'X14'),
    # Grab the masses of 12C and 13C from IsoSpec builtin data, provide mass of 14C manually
    mass = c(isotopes[isotopes$isotope == 'C12', 'mass'], isotopes[isotopes$isotope == 'C13', 'mass'], 14.003241989),
    # Isotope abundances are: 95% 14C, and the remaining 5% is split between 12C and 13C proportionally to their natural abundances
    abundance = c(isotopes[isotopes$isotope == 'C12', 'abundance'] * 0.05, isotopes[isotopes$isotope == 'C13', 'abundance'] * 0.05, 0.95),
    # unused parameter, kept for backward compatibility, ignore that
    ratioC = c(NA, NA, NA)
)

isotopes = rbind(isotopes, radiolabel)



radioglucose = c(C=5, H=12, O=6, X=1) # A glucose molecule: one of the carbons is replaced with an artificial element X (the readiocarbon)

p = .99 # joint probability of the output

IsoSpecify(radioglucose, p, isotopes = isotopes, showCounts = TRUE)


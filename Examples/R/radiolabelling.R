library(IsoSpecR)

### DOES NOT WORK YET

# We will model a glucose molecule where one of the carbons has been replaced with 
# (14C) with 95% probability.


# First, we define an alternate isotopic table for elements that contains the readiocarbon: 14C with 95% probability, 
# and either 12C or 13C with their standard relative probabilities otherwise. First, we load up standard isotopic data:
data(isotopicData)
isotopes = isotopicData$IsoSpec

# And we add an artificial element X, with 3 isotopes that will represent the radiocarbon. 
# This part is a bit messy...
levels(isotopes[[2]]) = c(levels(isotopes[[2]]), "X12", "X13", "X14")
isotopes = rbind(isotopes, c("X", "X12", 12.0, 0.05*0.989211941850467, 0))
isotopes = rbind(isotopes, c("X", "X13", 13.0033548352, 0.05*0.0107880581495331, 0))
isotopes = rbind(isotopes, c("X", "X14", 14.003241989, 0.95, 0))



radioglucose = c(C=5, H=12, O=6, X=1) # A glucose molecule: one of the carbons is replaced with an artificial element X

p = .99 # joint probability of the output

IsoSpecify(radioglucose, p, isotopes = isotopes)


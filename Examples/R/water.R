library(IsoSpecR)

# A water molecule:
water <- c(H=2,O=1)

# Desired joint probability p of the p-optimal set of isotopologues (90%): 
p <- .9

# Get raw information:
resRaw <- IsoSpecify( molecule=water, stopCondition=.99 )
print(resRaw)

# While not doing any high-throughput analysis, you might consider to present reported data in a fancy way. 
# ATTENTION: this will CORRUPT THE LINEAR RUNTIME claim, as setting "fancy = TRUE" involves sorting the results. 
res <- IsoSpecify( molecule=water, stopCondition=.99 )

# ATTENTION: while turned on, the algorithm's time complexity is nlog(n) instead of linear.


print('The first configuration has the following parameters:')
print('Mass:');res$mass
print('log(probability):');res$logProb
print('probability:');res$prob
print('Number of Protium atoms:');res$H1
print('Number of Deuterium atoms:');res$H2
print('Number of O16 atoms:');res$O16
print('Number of O17 atoms:');res$O17
print('Number of O18 atoms:');res$O18

print("Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?")
print('In R, we have to provide an additional parameter to the algorithm: a data.frame containing the new isotopic ratios.')
modifiedIsotopes <- data.frame(
	element = c('H', 'H', 'O', 'O', 'O'),
	isotope = c('H1', 'H2', 'O16', 'O17', 'O18'),
	mass  	= c(1.00782503207, 2.0141017778,15.99491461956, 16.99913170, 17.9991610),
	abundance = c(0.5, 0.5,0.5, 0.3, 0.2)
)

modRes <- IsoSpecify( molecule=water, stopCondition=.99, isotopes=modifiedIsotopes )

print('The number of configuration must be bigger, the probability being less concentrated on any isotope.')
modRes

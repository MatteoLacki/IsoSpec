from IsoSpecPy import IsoBinned

# Here we will compute a binned isotopic spectrum of a molecule
# C275H440N80O90S6, with resolution of 0.01 Da, and print the resulting bins
# The bins are centered starting at 0.0, which means that the first bin covers
# the range [-0.005 Da, 0.005 Da), the second one covers [0.005 Da, 0.015 Da), etc.
# The spectrum will be computed until cumulative probability of 0.9999 is reached.
formula = "C275H440N80O90S6"
resolution = 0.01  # in Da
bin_centers_start = 0.0  # in Da
precision = 0.9999
binned_spectrum = IsoBinned(bin_width=resolution,
                            bin_middle=bin_centers_start,
                            formula=formula,
                            target_total_prob=precision)

from matplotlib import pyplot as plt
plt.bar(list(binned_spectrum.masses),
        list(binned_spectrum.probs),
        width=resolution,
        align='center')
plt.xlabel('Mass (Da)')
plt.ylabel('Probability')
plt.title(f'Binned isotopic spectrum of {formula}')
plt.show()

# Now, let's compute the same at 1.0 Da resolution, but suppose you want the bins
# to be centered at 0.5, i.e., the first bin should cover the range [0.0 Da, 1.0 Da),
# the second one [1.0 Da, 2.0 Da), etc.
resolution = 1.0  # in Da
bin_centers_start = 0.5  # in Da
binned_spectrum = IsoBinned(bin_width=resolution,
                            bin_middle=bin_centers_start,
                            formula=formula,
                            target_total_prob=precision)

plt.bar(list(binned_spectrum.masses),
        list(binned_spectrum.probs),
        width=resolution,
        align='center')
plt.xlabel('Mass (Da)')
plt.ylabel('Probability')
plt.title(f'Binned isotopic spectrum of {formula} at 1 Da resolution')
plt.show()
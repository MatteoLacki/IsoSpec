## What is IsoSpec? 

IsoSpec is a **fine structure isotopic calculator**.

It has been presented in a paper in Analytical Chemistry in 2017.

Given a molecular formula and some probability threshold *0 < P < 1*, it will provide you with the smallest possible set of <a href="http://goldbook.iupac.org/I03351.html" target="_self">isotopologues</a> that are jointly *P* probable. 

### *Example*
Say, that you happen to be interested in the top 50% probable isotopologues of <a href="http://www.rcsb.org/pdb/explore.do?structureId=2zp6" target="_self">*Bovine Insulin*</a>. We happen to know, that its molecular formula is C<sub>254</sub>H<sub>377</sub>N<sub>65</sub>O<sub>75</sub>S<sub>6</sub>. 

In that case, IsoSpec would provide you with:

|   |            mass|  probability| <sup>1</sup>H| <sup>2</sup>H| <sup>12</sup>C| <sup>13</sup>C| <sup>14</sup>N| <sup>15</sup>N| <sup>16</sup>O| <sup>17</sup>O| <sup>18</sup>O| <sup>32</sup>S| <sup>33</sup>S| <sup>34</sup>S| <sup>36</sup>S|
|:--|---------------:|------------:|-------------:|-------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|
|1  | 5731.6075806688| 0.1123023514|           377|             0|            252|              2|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|2  | 5732.6109355040| 0.1028778936|           377|             0|            251|              3|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|3  | 5730.6042258336| 0.0814037470|           377|             0|            253|              1|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|4  | 5733.6142903392| 0.0704027660|           377|             0|            250|              4|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|5  | 5734.6176451744| 0.0383896060|           377|             0|            249|              5|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|7  | 5733.6033765247| 0.0301636876|           377|             0|            252|              2|             65|              0|             75|              0|              0|              5|              0|              1|              0|
|6  | 5729.6008709984| 0.0293871014|           377|             0|            254|              0|             65|              0|             75|              0|              0|              6|              0|              0|              0|
|8  | 5734.6067313599| 0.0276323390|           377|             0|            251|              3|             65|              0|             75|              0|              0|              5|              0|              1|              0|
|9  | 5732.6046155640| 0.0266824062|           377|             0|            252|              2|             64|              1|             75|              0|              0|              6|              0|              0|              0|


Now, how this could be of any use to **you**? Well, we did suppose you are rather into ...

### Mass Spectrometry


Say **you** hoped to know, if *Bovine Insulin* is present in the sample you analyzed with your <a href="https://en.wikipedia.org/wiki/Mass_spectrometry" target="_self">**mass spectrometer**</a>.
How would you do that? 

You would use our software to generate the isotopologues and then try to find them in the experimental spectrum.





## How to install IsoSpec?

IsoSpec is written in **C++** has bindings to **Python** (IsoSpecPy), **R** (IsoSpecR), and **C**. 

### Python
**IMPORTANT**: please note that the Python package is standalone, in the
sense that it does not need the C/C++ library to be installed separately.

Requirements:
* Python (v2 and v3 are 
* C++11 compiler: clang++ (v3.3 or later) or g++ (v4.7 or later)
* setuptools
* cffi (this, if not present, will be automatically downloaded and
          installed by the setup script, however you may prefer to use 
          your distribution's package manager to install that)


```
pip install IsoSpecPy
```

alternatively, you can download our package from here and then

```
cd IsoSpecPy
sudo python setup.py install
```

Again, clang++ is the preferred compiler and will be used if found by the 
setup script. If you want to override the behaviour (if you have clang++ 
installed, but still want to use g++) you will have to replace the last 
command with:

```
ISO_USE_DEFAULT_CXX=TRUE sudo python setup.py install
```

### R

Requirements:
* R (>= 3.2.1)

The package is hosted at <a href="https://cran.r-project.org/web/packages/IsoSpecR/index.html" target="_self">**CRAN**</a>. 
This means that it can be automatically downloaded. Just start an R console (or R studio) and run

```R
    install.packages('IsoSpecR')
```

Then, follow the instructions. For Windows users, this will result in downloading a precompiled version of the package.

The package can be also directly downloaded from this page. If you use either Linux of Mac OSX, then simply:

1. Download the package.
2. Move to the folder containing the IsoSpecR folder. 
3. Run in terminal

```
	R CMD build IsoSpecR 
	R CMD INSTALL IsoSpecR_1.0.tar.gz  
```

All necessary packages should download automatically. 


### C/C++

Requirements:

* C++11 compiler: clang++ (v3.3 or later) or g++ (v4.7 or later)
* Make (GNU make is OK)


Note: clang++ is the default (and preferred) compiler as it produces 
faster code (in our tests, the difference in speed is about 20%). 
If you'd like to use g++ instead, please edit the first line of 
IsoSpec++/Makefile appropriately.

Next, execute the following commands:

```
cd IsoSpec++
make
```

You may copy the resulting .so file to a convenient location.
If you wish to develop software using this library, you will also 
have to place the header files (*.hpp/*.h) somewhere your C/C++
compiler can find them.

## Here are some example sessions:

### Python

```python
# Calculates the isotopic distribution of water in several ways

from IsoSpecPy import IsoSpecPy
from math import exp

i = IsoSpecPy.IsoSpec.IsoFromFormula("H2O1", 0.9)

print "The isotopologue set containing at least 0.9 probability has", len(i), "element(s)"

confs = i.getConfs()

print "The first configuration has the following parameters:"
print "Mass:", confs[0][0]
print "log(probability):", confs[0][1] 
print "probability:", exp(confs[0][1])
print "Number of Protium atoms:", confs[0][2][0][0]
print "Number of Deuterium atoms", confs[0][2][0][1]
print "Number of O16 atoms:", confs[0][2][1][0]
print "Number of O17 atoms:", confs[0][2][1][1]
print "Number of O18 atoms:", confs[0][2][1][2]

print
print "Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?"

hydrogen_probs = (0.5, 0.5)
oxygen_probs = (0.5, 0.3, 0.2)
hydrogen_masses = (1.00782503207, 2.0141017778)
oxygen_masses = (15.99491461956, 16.99913170, 17.9991610)
atom_counts = (2, 1)

i = IsoSpecPy.IsoSpec(atom_counts, (hydrogen_masses, oxygen_masses), (hydrogen_probs, oxygen_probs), 0.9)

print "The isotopologue set containing at least 0.9 probability has", len(i), "element(s)"

confs = i.getConfs()

print "The first configuration has the following parameters:"
print "Mass:", confs[0][0]
print "log-prob:", confs[0][1]
print "probability:", exp(confs[0][1])
print "Number of Protium atoms:", confs[0][2][0][0]
print "Number of Deuterium atoms", confs[0][2][0][1]
print "Number of O16 atoms:", confs[0][2][1][0]
print "Number of O17 atoms:", confs[0][2][1][1]
print "Number of O18 atoms:", confs[0][2][1][2]
```

### R

```R
library(IsoSpecR)

# A water molecule:
water <- c(H=2,O=1)

# Desired joint probability p of the p-optimal set of isotopologues (90%): 
p <- .9

# The fancy representation of the results is on.
# ATTENTION: while turned on, the algorithm's time complexity is nlog(n) instead of linear.
res <- IsoSpecify( molecule=water, stopCondition=.99, fancy=TRUE )

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
print('In R, we have to preper additional parameter for the algorithm: a data.frame containing the new isotopic ratios.')
modifiedIsotopes <- data.frame(
	element = c('H', 'H', 'O', 'O', 'O'),
	isotope = c('H1', 'H2', 'O16', 'O17', 'O18'),
	mass  	= c(1.00782503207, 2.0141017778,15.99491461956, 16.99913170, 17.9991610),
	abundance = c(0.5, 0.5,0.5, 0.3, 0.2)
)

modRes <- IsoSpecify( molecule=water, stopCondition=.99, fancy=TRUE, isotopes=modifiedIsotopes )

print('The number of configuration must be bigger, the probability being less concentrated on any isotope.')
modRes
```

### C++

```C++

#include <iostream>
#include "../../IsoSpec++/isoSpec++.hpp"


int main()
{
    IsoSpec* iso = IsoSpec::IsoFromFormula("H2O1", 0.9);

    iso->processConfigurationsUntilCutoff();

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << iso->getNoVisitedConfs() << " element(s)" << std::endl;

    std::tuple<double*,double*,int*,int> product = iso->getCurrentProduct();

    double* masses = std::get<0>(product);
    double* logprobs = std::get<1>(product);
    int* configs = std::get<2>(product);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << masses[0] << std::endl;
    std::cout << "log-prob: " << logprobs[0] << std::endl;
    std::cout << "probability: " << exp(logprobs[0]) << std::endl;

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;

    delete iso;
    delete masses;
    delete logprobs;
    delete configs;


    std::cout << "Now what if both isotopes of hydrogen were equally probable, while prob. of O16 was 50%, O17 at 30% and O18 at 20%?" << std::endl;

    const int elementNumber = 2;
    const int isotopeNumbers[2] = {2,3};

    const int atomCounts[2] = {2,1};


    const double hydrogen_masses[2] = {1.00782503207, 2.0141017778};
    const double oxygen_masses[3] = {15.99491461956, 16.99913170, 17.9991610};

    const double* isotope_masses[2] = {hydrogen_masses, oxygen_masses};

    const double hydrogen_probs[2] = {0.5, 0.5};
    const double oxygen_probs[3] = {0.5, 0.3, 0.2};

    const double* probs[2] = {hydrogen_probs, oxygen_probs};



    iso = new IsoSpecLayered(elementNumber, isotopeNumbers, atomCounts, isotope_masses, probs, 0.9);

    iso->processConfigurationsUntilCutoff();

    std::cout <<  "The isotopologue set containing at least 0.9 probability has " << iso->getNoVisitedConfs() << " element(s)" << std::endl;

    product = iso->getCurrentProduct();

    masses = std::get<0>(product);
    logprobs = std::get<1>(product);
    configs = std::get<2>(product);

    std::cout << "The first configuration has the following parameters: " << std::endl;
    std::cout << "Mass: " << masses[0] << std::endl;
    std::cout << "log-prob: " << logprobs[0] << std::endl;
    std::cout << "probability: " << exp(logprobs[0]) << std::endl;

    // Successive isotopologues are ordered by the appearance in the formula of the element, then by nucleon number, and concatenated into one array
    std::cout << "Protium atoms: " << configs[0] << std::endl;
    std::cout << "Deuterium atoms " << configs[1] << std::endl;
    std::cout << "O16 atoms: " << configs[2] << std::endl;
    std::cout << "O17 atoms: " << configs[3] << std::endl;
    std::cout << "O18 atoms: " << configs[4] << std::endl;

    delete iso;
    delete masses;
    delete logprobs;
    delete configs;
}
```

## Interested in IsoSpec?

Contact us! Mail to:

* matteo.lacki@gmail.com 

or 

* mist@mimuw.edu.pl

and we might help you in including our software in your project!

Please note, that we are giving you IsoSpec for free.
Please, be kind-hearted, and if you are to use our software in your research, do not forget to cite us:

Łącki, M. K., Startek, M., Valkenborg, D., & Gambin, A. (2017). IsoSpec: Hyperfast Fine Structure Calculator. Analytical Chemistry, 89(6), 3272-3277.

More information on the paper is <a href="http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b01459" target="_self">**HERE**</a>

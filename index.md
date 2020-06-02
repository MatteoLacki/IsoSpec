## What is IsoSpec? 

IsoSpec is a **fine structure isotopic calculator**.
It can reveal for you a given fraction of the isotopic distribution of mass for a given molecule based only on its molecular composition.
Mass is not unique due to the presence of isotopes that occur randomly in Nature, albeit with well studied frequencies.
Say you want to reveal *P = 50 %* of the masses of <a href="http://www.rcsb.org/pdb/explore.do?structureId=2zp6" target="_self">*Bovine Insulin*</a>.
We happen to know, that its molecular formula is C<sub>254</sub>H<sub>377</sub>N<sub>65</sub>O<sub>75</sub>S<sub>6</sub>. 

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

Before getting into details, a small break for a commercial!

### Using IsoSpec? Cite us! Citations are almost like money. Or are they? Well, if you do want to pay us, we will not be angry, no, we will gracefully accept every donation too.

<a href="http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b01459" target="_self">
IsoSpec: Hyperfast Fine Structure Calculator
Mateusz K. Łącki, Michał Startek, Dirk Valkenborg, and Anna Gambin
Analytical Chemistry 2017 89 (6), 3272-3277
DOI: 10.1021/acs.analchem.6b01459
</a>

## How to install IsoSpec?

IsoSpec is written in **C++** and has bindings to **Python** (IsoSpecPy), **R** (IsoSpecR), and **C**. 

### Python

IsoSpecPy package is standalone and does not need the C/C++ library to be installed separately.

Requirements:
* Python (v2 and v3 are 
* C++11 compiler: clang++ (v3.3 or later) or g++ (v4.7 or later)
* setuptools

Follow <a href="https://wiki.python.org/moin/WindowsCompilers">these instructions</a> to install a C++ compiler on Windows.

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

If you use either Linux of Mac OSX, then simply:

1. Download the package <a href="https://github.com/MatteoLacki/IsoSpec">from github</a>.
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

# How to run IsoSpec?

## Getting a coverage of the isotopic distribution

To obtain a given coverage, say *P = 99.99%*, of the isotopic distribution, follow these steps.

### Python

```{python}
import IsoSpecPy as iso

# P = 99.99% = 0.9999 = .9999
sp = iso.IsoTotalProb(formula="C254H377N65O75S6", prob_to_cover=.9999)
```
This performs the calculations and dumps the results to RAM.
You can easily access the masses and probabilities by simple iteration:
```{python}
for mass, prob in sp:
    print(mass, prob)
```
or access the cffi arrays directly:
```{python}
print(sp.masses[0])
print(sp.probs[0])
```
The first 0-th isotopologue is always one of the mode(s) of the isotopic distribution.
If you have installed *matplotlib* (separately), you can also easily visualize the outcome
```{python}
sp.plot()
```


Also, if you want to obtain only peaks above a certain heigh, simply run
```{python}
sp = iso.IsoThreshold(formula="C254H377N65O75S6", threshold=.0001)
```
This way, *sp* contains only peaks with individual probabilities heigher than *0.01%*.

The newly added nice feature of IsoSpec2.1 is that you don't need to store these isotopologues in RAM.
You can also iterate over them using the *generator mode*:
```{python}
for m, p in iso.IsoThresholdGenerator(formula="C254H377N65O75S6", threshold=.0001):
    print(m, p)
```

### R

We do not like **R**. But you might. We accept that, as much as we tolerate vegetarians.
**R** does not let us do as much as Python does.

To get here *99.99%* of the isotopic distribution, simply run
```R
library(IsoSpecR)

X = IsoSpecify(molecule=c(C=254,H=377,N=65,O=75,S=6), 
	       stopCondition=.9999, 
               showCounts=True)
print(X)
```
The above result contains the isotopic counts.
If you don't need them, simply run
```{R}
X = IsoSpecify(molecule=c(C=254,H=377,N=65,O=75,S=6), 
	       stopCondition=.9999)
print(X)
```

And to get all isotopologues with probabilities higher than *0.01%*, use:

```{R}
X = IsoSpecify(molecule=c(C=254,H=377,N=65,O=75,S=6),
	       stopCondition=.0001,
	       algo=2)
print(X)
```
You can also pass in custom isotopic frequencies.
Simply, pass in a data-frame as the `isotopes` parameter.
Have a look at `data(isotopicData)`, which is a list of examples of a proper input.

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
* mist@mimuw.edu.pl

and we will help you in intergrate our software in your project!

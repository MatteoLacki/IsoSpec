#include <iostream>
#include "/Users/matteo/Documents/IsoSpec/IsoSpec/IsoSpec++/unity-build.cpp"

using std::cout;
using std::endl;

int main(void){
    int isotopeNumbers[] = {2, 3};
    int config_size = isotopeNumbers[0] + isotopeNumbers[1];
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* generator = setupIsoThresholdGenerator(
        2,
        isotopeNumbers,
        atomCounts,
        isotopeMasses,
        isotopeProbabilities,
        .001,
        true,
        1000,
        1000);

    Tabulator tabulator(generator, true, true, true, true);

    double *masses = tabulator.masses();

    cout << masses[0] << endl;

    return 0;
}

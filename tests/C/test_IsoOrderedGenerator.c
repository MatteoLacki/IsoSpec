#include <iostream>
#include "../../IsoSpec++/unity-build.cpp"

using std::cout;

int main()
{
    int isotopeNumbers[] = {2, 3};
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* p = setupIsoOrderedGenerator(
        2,
        isotopeNumbers,
        atomCounts,
        isotopeMasses,
        isotopeProbabilities,
        1000,
        1000);

    while(advanceToNextConfigurationIsoOrderedGenerator(p))
    {
        cout << "mass="<< massIsoOrderedGenerator(p) << " lprob=" <<
        lprobIsoOrderedGenerator(p) << std::endl;
    }

    deleteIsoOrderedGenerator(p);

    return 0;
}

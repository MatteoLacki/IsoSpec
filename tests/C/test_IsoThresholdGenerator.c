#include <iostream>
#include "../../IsoSpec++/unity-build.cpp"

using std::cout;

int main()
{
    int isotopeNumbers[] = {2, 3};
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* p = setupIsoThresholdGenerator(
        2,
        isotopeNumbers,
        atomCounts,
        isotopeMasses,
        isotopeProbabilities,
        .001,
        true,
        1000,
        1000);

    while(advanceToNextConfigurationIsoThresholdGenerator(p))
    {
        cout << "mass="<< massIsoThresholdGenerator(p) << " lprob=" <<
        lprobIsoThresholdGenerator(p) << std::endl;
    }

    deleteIsoThresholdGenerator(p);

    return 0;
}

#include <iostream>
#include "../../IsoSpec++/unity-build.cpp"

using std::cout;

int main()
{
    int isotopeNumbers[] = {2, 3};
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* iso = setupIso(2, isotopeNumbers, atomCounts, isotopeMasses, isotopeProbabilities);

    void* p = setupIsoThresholdGenerator(
        iso,
        .001,
        true,
        1000,
        1000,
        true);

    while(advanceToNextConfigurationIsoThresholdGenerator(p))
    {
        cout << "mass="<< massIsoThresholdGenerator(p) << " lprob=" <<
        lprobIsoThresholdGenerator(p) << std::endl;
    }

    deleteIsoThresholdGenerator(p);
    deleteIso(iso);

    return 0;
}

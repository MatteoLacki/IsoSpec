#include <iostream>
#include "../../IsoSpec++/unity-build.cpp"

using std::cout;
using std::endl;

using namespace IsoSpec;

int main()
{
    int isotopeNumbers[] = {2, 3};
    int config_size = isotopeNumbers[0] + isotopeNumbers[1];
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* iso = setupIso(
        2,
        isotopeNumbers,
        atomCounts,
        isotopeMasses,
        isotopeProbabilities);

    void* p = setupIsoThresholdGenerator(
        iso,
        .001,
        true,
        1000,
        1000);

    int conf_no(0);
    int *space = new int[config_size];
    while(advanceToNextConfigurationIsoThresholdGenerator(p))
    {
        cout << "mass="<< massIsoThresholdGenerator(p) << " lprob=" <<
        lprobIsoThresholdGenerator(p) << endl;

        get_conf_signatureIsoThresholdGenerator(p, space);
        for(conf_no = 0; conf_no < config_size; conf_no++)
        {
            cout << space[conf_no] << " ";
        }

        cout << endl;

    }

    deleteIsoThresholdGenerator(p);
    delete[] space;
    deleteIso(iso);

    return 0;
}

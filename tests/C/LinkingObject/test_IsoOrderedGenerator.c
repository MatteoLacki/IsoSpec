#include <stdio.h>
#include "../../../IsoSpec++/cwrapper.h"

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
        .5,
        false,
        1000,
        1000);

    while(advanceToNextConfigurationIsoThresholdGenerator(p))
    {
        printf("mass:  %f\n", massIsoThresholdGenerator(p));
        printf("lprob: %f\n", lprobIsoThresholdGenerator(p));
    }

    deleteIsoThresholdGenerator(p);

    return 0;
}

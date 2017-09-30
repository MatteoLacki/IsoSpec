#include <stdio.h>
#include "../../../IsoSpec++/cwrapper.h"

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
        printf("mass:  %f\n", massIsoOrderedGenerator(p));
        printf("lprob: %f\n", lprobIsoOrderedGenerator(p));
    }

    deleteIsoOrderedGenerator(p);

    return 0;
}

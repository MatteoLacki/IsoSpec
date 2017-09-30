#include <stdio.h>
#include "../../../IsoSpec++/cwrapper.h"

int main()
{
    int isotopeNumbers[] = {2, 3};
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* p = setupIsoLayeredGenerator(
        2,
        isotopeNumbers,
        atomCounts,
        isotopeMasses,
        isotopeProbabilities,
        -3.0,
        1000,
        1000);

    while(advanceToNextConfigurationIsoLayeredGenerator(p))
    {
        printf("mass:  %f\n", massIsoLayeredGenerator(p));
        printf("lprob: %f\n", lprobIsoLayeredGenerator(p));
    }

    deleteIsoLayeredGenerator(p);

    return 0;
}

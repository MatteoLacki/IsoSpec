#include <iostream>
#include "/Users/matteo/Documents/IsoSpec/IsoSpec/IsoSpec++/unity-build.cpp"

using std::cout;
using std::endl;

int main()
{
    cout << "Welcome to the test!" << endl;
    int isotopeNumbers[] = {2, 3};
    int atomCounts[] = {10, 10};
    double isotopeMasses[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double isotopeProbabilities[] = {0.5, 0.5, 0.5, 0.3, 0.2};

    void* p = setupIsoOrderedGenerator(2,
                                       isotopeNumbers,
                                       atomCounts,
                                       isotopeMasses,
                                       isotopeProbabilities,
                                       1000,
                                       1000);

    MassSpectrum MS = set_tablesIsoOrderedGenerator(p, 100);

    double *masses = MS.masses, *lprobs = MS.logprobs;
    int confs_no = MS.confs_no;

    for(int i = 0; i < confs_no; i++ )
    {
        cout << masses[i] << " " << lprobs[i] << endl;
    }

    deleteIsoOrderedGenerator(p);

    return 0;
}

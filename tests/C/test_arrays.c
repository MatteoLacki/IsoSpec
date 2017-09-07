#include <iostream>
#include "/Users/matteo/Documents/IsoSpec/IsoSpec/IsoSpec++/unity-build.cpp"

using std::cout;
using std::endl;

int main()
{
    cout << "Hello" << endl;
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

    double *masses, *lprobs;
    int config_no;
    cout << config_no << " "<< endl;
    set_tablesIsoOrderedGenerator(p, &masses, &lprobs, &config_no, 100);

    cout << config_no << endl;
    for(int i = 0; i < config_no; i++ )
    {
        cout << masses[i] << " " << lprobs[i] << endl;
    }

    deleteIsoOrderedGenerator(p);

    return 0;
}

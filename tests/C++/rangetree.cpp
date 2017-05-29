#include "marginalTrek++.h"


const double masses[] = {10000.0, 100.0, 1.0, 1.0};
const double probs[] = {0.5, 0.4, 0.1, 0.001};

const double lCutoff = -1000000.0;
int main()
{
Marginal m(masses, probs, 4, 5);
PrecalculatedMarginal PC(std::move(m), lCutoff, true);
RGTMarginal* RGTM = new RGTMarginal(std::move(m), lCutoff);
unsigned int matches = 0;
unsigned int mismatches = 0;
for(unsigned int ii = 0; ii < PC.get_no_confs(); ii++)
    for(unsigned int jj = 0; jj < PC.get_no_confs(); jj++)
        for(unsigned int kk = 0; kk < PC.get_no_confs(); kk++)
            for(unsigned int ll = 0; ll < PC.get_no_confs(); ll++)
            {
                double pmax = PC.get_lProb(ii);
                double pmin = PC.get_lProb(jj);
                double mmin = PC.get_mass(kk);
                double mmax = PC.get_mass(ll);
                unsigned int simple = 0;
                for(unsigned int mm = 0; mm < PC.get_no_confs(); mm++)
                    if(pmin <= PC.get_lProb(mm) and pmax >= PC.get_lProb(mm) and mmin <= PC.get_mass(mm) and mmax >= PC.get_mass(mm))
                        simple++;
                unsigned int hard = 0;
                RGTM->setup_search(pmin, pmax, mmin, mmax);
                while(RGTM->next())
                    hard++;
                if (simple == hard)
                    matches += 1;
                else
                {
                    mismatches += 1;
                    std::cout << "PARAMS p: " << pmin << " " << pmax << "m: " << mmin << " " << mmax << std::endl;
                    std::cout << "COMPARE " << simple << " " << hard << std::endl;
                    std::cout << "MISMATCH" << std::endl;
                }
            }

/*while(RGTM->next())
{
Conf c = RGTM->current_conf();
std::cout << RGTM->current_lProb() << "\t" << exp(RGTM->current_lProb()) << "\t" << RGTM->current_mass() << '\t' << c[0] << '\t' << c[1] << '\t' << c[2] <<std::endl;
}
*/
std::cout << "matches/mismatches: " << matches << " " << mismatches << std::endl;
delete RGTM;
/*
SyncMarginal* SM = new SyncMarginal(Marginal(masses, probs, 3, 5), -100.0);
unsigned int idx = SM->getNextConfIdxwMass(1000.0, 10000000.0);
while(SM->inRange(idx))
{
    if(SM->get_lProb(idx)>=-3.0)
        std::cout << SM->get_lProb(idx) << "\t" << SM->get_mass(idx) << '\n';
    idx = SM->getNextConfIdxwMass(1000.0, 10000000.0);
}*/
}

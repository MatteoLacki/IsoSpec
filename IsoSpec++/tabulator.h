#ifndef __TABULATOR_H__
#define __TABULATOR_H__

#include "isoSpec++.h"

template <typename T> class Tabulator
{
private:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    int64_t     _confs_no;
public:
    Tabulator(T* generator,
              bool get_masses, bool get_probs,
              bool get_lprobs, bool get_confs);

    ~Tabulator();

    inline double*   masses()   { return _masses; };
    inline double*   lprobs()   { return _lprobs; };
    inline double*   probs()    { return _probs; };
    inline int*      confs()    { return _confs; };
    inline int64_t       confs_no() { return _confs_no; };
};

#endif  // __TABULATOR_H__

#pragma once

#include <stdlib.h>

#include "isoSpec++.h"

#define ISOSPEC_INIT_TABLE_SIZE 1024

namespace IsoSpec
{

//! The class used to store configurations from generators.
template <typename T> class Tabulator
{
private:
    double* _masses;
    double* _lprobs;
    double* _probs;
    int*    _confs;
    size_t  _confs_no;
public:
    //! Constructor
    /*!
        \param generator  An instance of a generator class (or its subclass).
        \param get_masses If true, the tabulator stores the masses of the isotopologues.
        \param get_probs  If true, the tabulator stores the probabilities of the isotopologues.
        \param get_lprobs If true, the tabulator stores the log-probabilities of the isotopologues.
        \param get_confs  If true, the tabulator stores the counts of isotopes corresponding to each isotopologue.
    */
    Tabulator(T* generator,
              bool get_masses, bool get_probs,
              bool get_lprobs, bool get_confs);

    //! Destructor.
    ~Tabulator();

    //! Get the masses of the calculated isotopologues.
    inline double*   masses()   { return _masses; };

    //! Get the log-probabilities of the calculated isotopologues.
    inline double*   lprobs()   { return _lprobs; };

    //! Get the probabilities of the calculated isotopologues.
    inline double*   probs()    { return _probs; };

    //! Get the counts of isotopes that make up the calculated isotopologues.
    inline int*      confs()    { return _confs; };

    //! Get the number of calculated isotopologues.
    inline size_t    confs_no() { return _confs_no; };
};


//! Reallocate the array to a bigger one.
/*!
    \param array An array of arrays to reallocate.
    \param new_size The size of the new array.
*/
void reallocate(double **array, int new_size){
    if( *array != nullptr ){
        *array = (double *) realloc(*array, new_size);
    }
}


template <typename T> Tabulator<T>::Tabulator(T* generator,
                     bool get_masses, bool get_probs,
                     bool get_lprobs, bool get_confs  )
{
    size_t current_size = ISOSPEC_INIT_TABLE_SIZE;
    int confs_tbl_idx = 0;
    _confs_no = 0;

    const int allDimSizeOfInt = sizeof(int)*generator->getAllDim();

    _masses = get_masses ? (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _lprobs = get_lprobs ? (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _probs  = get_probs  ? (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _confs  = get_confs  ? (int *)    malloc(ISOSPEC_INIT_TABLE_SIZE * allDimSizeOfInt): nullptr;


    while(generator->advanceToNextConfiguration()){
        if( _confs_no == current_size )
        {
            current_size *= 2;

	    // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...

            reallocate(&_masses, current_size * sizeof(double));
            reallocate(&_lprobs, current_size * sizeof(double));
            reallocate(&_probs,  current_size * sizeof(double));

            if( _confs != nullptr ){
                _confs = (int *) realloc(_confs, current_size * allDimSizeOfInt);
            }
        }

        if(_masses != nullptr) _masses[_confs_no] = generator->mass();

        if(_lprobs != nullptr) _lprobs[_confs_no] = generator->lprob();

        if(_probs  != nullptr) _probs[_confs_no]  = generator->eprob();

        if(_confs  != nullptr){
            generator->get_conf_signature(_confs + confs_tbl_idx);
            confs_tbl_idx += generator->getAllDim();
        }

        _confs_no++;
    }

    _masses = (double *) realloc(_masses, _confs_no * sizeof(double));
    _lprobs = (double *) realloc(_lprobs, _confs_no * sizeof(double));
    _probs  = (double *) realloc(_probs,  _confs_no * sizeof(double));
    _confs  = (int *)    realloc(_confs,  confs_tbl_idx * sizeof(int));
}

template <typename T> Tabulator<T>::~Tabulator()
{
    if( _masses != nullptr ) free(_masses);
    if( _lprobs != nullptr ) free(_lprobs);
    if( _probs  != nullptr ) free(_probs);
    if( _confs  != nullptr ) free(_confs);
}

} // namespace IsoSpec


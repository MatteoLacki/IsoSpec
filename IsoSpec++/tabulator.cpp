#include "tabulator.h"

#define INIT_TABLE_SIZE 1024

#define PUSH_TO_ARRAY(array, input) \
{ \
    if (confs_no == current_size) \
    { \
        current_size *= 2; \
        array = (double *) realloc(array, current_size * sizeof(double)); \
    } \
    array[confs_no] = input; \
    confs_no++; \
}

#define MEMCPY_TO_ARRAY(array, input) \
{ \
    if (confs_no == current_size) \
    { \
        current_size *= 2; \
        array = (double *) realloc(array, current_size * sizeof(double)); \
    } \
    memcpy(array+confs_tbl_idx, generator->get_conf_signature(), allDimSizeOfInt) \
    confs_tbl_idx += generator->getAllDim(); \
}

void reallocate(double **array, int new_size){
    if( *array != nullptr ){
        *array = (double *) realloc(*array, new_size);
    }
}

Tabulator::Tabulator(IsoThresholdGenerator* generator,
                     bool get_masses, bool get_probs,
                     bool get_lprobs, bool get_confs  )
{
    int current_size = INIT_TABLE_SIZE;
    int current_size_dims = INIT_TABLE_SIZE * generator->getAllDim();
    int confs_tbl_idx = 0;
    int confs_no = 0;

    const int allDimSizeOfInt = sizeof(int)*generator->getAllDim();

    _masses = get_masses ? (double *) malloc(INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _lprobs = get_lprobs ? (double *) malloc(INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _probs  = get_probs  ? (double *) malloc(INIT_TABLE_SIZE * sizeof(double)) : nullptr;
    _confs  = get_confs  ? (int *)    malloc(INIT_TABLE_SIZE * allDimSizeOfInt): nullptr;

    while(advanceToNextConfigurationIsoThresholdGenerator(generator)){
        if( confs_no == current_size )
        {
            current_size *= 2;
            reallocate(&_masses, current_size * sizeof(double));
            reallocate(&_lprobs, current_size * sizeof(double));
            reallocate(&_probs,  current_size * sizeof(double));

            if( _confs != nullptr ){
                _confs = (int *) realloc(_confs, current_size * allDimSizeOfInt);
            }
        }

        if(_masses != nullptr) _masses[confs_no] = generator->mass();

        if(_lprobs != nullptr) _lprobs[confs_no] = generator->lprob();

        if(_probs  != nullptr) _probs[confs_no]  = generator->eprob();

        if(_confs  != nullptr){
            memcpy(_confs + confs_tbl_idx, generator->get_conf_signature(), allDimSizeOfInt);
            confs_tbl_idx += generator->getAllDim();
        }

        confs_no++;
    }

    _masses = (double *) realloc(_masses, confs_no * sizeof(double));
    _lprobs = (double *) realloc(_lprobs, confs_no * sizeof(double));
    _probs  = (double *) realloc(_probs,  confs_no * sizeof(double));
    _confs  = (int *)    realloc(_confs,  confs_tbl_idx * sizeof(int));
}

Tabulator::~Tabulator()
{
    if( _masses != nullptr ) free(_masses);
    if( _lprobs != nullptr ) free(_lprobs);
    if( _probs  != nullptr ) free(_probs);
    if( _confs  != nullptr ) free(_confs);
}

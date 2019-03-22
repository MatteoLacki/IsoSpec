#include "tabulator.h"


namespace IsoSpec
{

Tabulator::Tabulator() : 
_masses(nullptr),
_lprobs(nullptr),
_probs(nullptr),
_confs(nullptr),
_confs_no(0)
{}

Tabulator::~Tabulator()
{
    if( _masses != nullptr ) free(_masses);
    if( _lprobs != nullptr ) free(_lprobs);
    if( _probs  != nullptr ) free(_probs);
    if( _confs  != nullptr ) free(_confs);
}


ThresholdTabulator::ThresholdTabulator(IsoThresholdGenerator* generator, 
                                       bool get_masses, bool get_probs,
                                       bool get_lprobs, bool get_confs) : 
Tabulator()
{
    _confs_no = generator->count_confs();
    const int allDim = generator->getAllDim();

    if(get_masses) _masses = (double *) malloc(_confs_no * sizeof(double));
    if(get_lprobs) _lprobs = (double *) malloc(_confs_no * sizeof(double));
    if(get_probs)  _probs  = (double *) malloc(_confs_no * sizeof(double));
    if(get_confs)  _confs  = (int *)    malloc(_confs_no * allDim*sizeof(int));

    while(generator->advanceToNextConfiguration())
    {
        if(_masses != nullptr) { *_masses = generator->mass();  _masses++; }
        if(_lprobs != nullptr) { *_lprobs = generator->lprob(); _lprobs++; }
        if(_probs  != nullptr) { *_probs  = generator->prob();  _probs++;  }
        if(_confs  != nullptr){
            generator->get_conf_signature(_confs);
            _confs += allDim;
        }
    }

    if(_masses != nullptr) { _masses -= _confs_no; }
    if(_lprobs != nullptr) { _lprobs -= _confs_no; }
    if(_probs  != nullptr) { _probs -= _confs_no; }
    if(_confs  != nullptr) { _confs -= _confs_no*allDim; }

}

ThresholdTabulator::~ThresholdTabulator() {}


LayeredTabulator::LayeredTabulator(IsoLayeredGenerator* ILG,
                     bool get_masses, bool get_probs,
                     bool get_lprobs, bool get_confs,
                     double _target_total_prob, bool _optimize) : 
Tabulator(),
generator(ILG),
target_total_prob(_target_total_prob),
current_size(ISOSPEC_INIT_TABLE_SIZE),
optimize(_optimize),
allDim(generator->getAllDim())
{
    if(_target_total_prob <= 0.0)
        return;

    bool user_wants_probs = get_probs;
    if(optimize)
    // If we want to optimize, we need the probs
        get_probs = true;

    if(get_masses) _masses = (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double));
    if(get_lprobs) _lprobs = (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double));
    if(get_probs)  _probs  = (double *) malloc(ISOSPEC_INIT_TABLE_SIZE * sizeof(double));
    if(get_confs)  _confs  = (int *)    malloc(ISOSPEC_INIT_TABLE_SIZE * generator->getAllDim() * sizeof(int));

    tmasses = _masses;
    tlprobs = _lprobs;
    tprobs = _probs;
    tconfs = _confs;

    size_t last_switch = 0;
    double prob_at_last_switch = 0.0;
    double prob_so_far = 0.0;

    do 
    { // Store confs until we accumulatoe more prob than needed - and, if optimizing, 
      // store also the rest of the layer
        while(generator->advanceToNextConfigurationWithinLayer())
        {
            addConf();
            prob_so_far += generator->prob();
            if(!optimize && prob_so_far >= target_total_prob)
                return;
        }
        if(prob_so_far >= target_total_prob)
            break;
        last_switch = _confs_no;
        prob_at_last_switch = prob_so_far;
    } while(generator->nextLayer(-3.0));

    if(!optimize)
        return;

    throw std::logic_error("hahaha, lol, nope, not yet implemented");

    if(!user_wants_probs)
    {
        free(_probs);
        _probs = nullptr;
    }
}



void LayeredTabulator::addConf()
{
    if( _confs_no == current_size )
    {
        current_size *= 2;

        // FIXME: Handle overflow gracefully here. It definitely could happen for people still stuck on 32 bits...

        if(_masses != nullptr) { _masses = (double*) realloc(_masses, current_size * sizeof(double)); tmasses = _masses + _confs_no; }
        if(_lprobs != nullptr) { _lprobs = (double*) realloc(_lprobs, current_size * sizeof(double)); tlprobs = _lprobs + _confs_no; }
        if(_probs  != nullptr) { _probs  = (double*) realloc(_probs, current_size * sizeof(double));  tprobs  = _probs  + _confs_no; }
        if( _confs != nullptr) { _confs  = (int*)    realloc(_confs, current_size * allDim * sizeof(int)); tconfs = _confs + (allDim * _confs_no); }
    }

    if(_masses != nullptr) { *tmasses = generator->mass();  tmasses++; };
    if(_lprobs != nullptr) { *tlprobs = generator->lprob(); tlprobs++; };
    if(_probs  != nullptr) { *tprobs  = generator->prob();  tprobs++;  };
    if(_confs  != nullptr) { generator->get_conf_signature(tconfs); tconfs += allDim; };

    _confs_no++;
}

LayeredTabulator::~LayeredTabulator() {}

} // namespace IsoSpec

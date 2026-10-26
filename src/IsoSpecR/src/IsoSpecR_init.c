#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP _IsoSpecR_Rinterface(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_RParsePeptideSequence(SEXP, SEXP);
extern SEXP _IsoSpecR_RIsotopicTable();
extern SEXP _IsoSpecR_RBinned(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_RStochastic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_RIsoProperties(SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_REnvelopeConvolve(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_REnvelopeLinearCombination(SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_REnvelopeBin(SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_REnvelopeResample(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_RWasserstein(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _IsoSpecR_RWassersteinMatch(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_IsoSpecR_Rinterface", (DL_FUNC) &_IsoSpecR_Rinterface, 10},
    {"_IsoSpecR_RParsePeptideSequence", (DL_FUNC) &_IsoSpecR_RParsePeptideSequence, 2},
    {"_IsoSpecR_RIsotopicTable", (DL_FUNC) &_IsoSpecR_RIsotopicTable, 0},
    {"_IsoSpecR_RBinned", (DL_FUNC) &_IsoSpecR_RBinned, 6},
    {"_IsoSpecR_RStochastic", (DL_FUNC) &_IsoSpecR_RStochastic, 8},
    {"_IsoSpecR_RIsoProperties", (DL_FUNC) &_IsoSpecR_RIsoProperties, 3},
    {"_IsoSpecR_REnvelopeConvolve", (DL_FUNC) &_IsoSpecR_REnvelopeConvolve, 4},
    {"_IsoSpecR_REnvelopeLinearCombination", (DL_FUNC) &_IsoSpecR_REnvelopeLinearCombination, 3},
    {"_IsoSpecR_REnvelopeBin", (DL_FUNC) &_IsoSpecR_REnvelopeBin, 4},
    {"_IsoSpecR_REnvelopeResample", (DL_FUNC) &_IsoSpecR_REnvelopeResample, 5},
    {"_IsoSpecR_RWasserstein", (DL_FUNC) &_IsoSpecR_RWasserstein, 7},
    {"_IsoSpecR_RWassersteinMatch", (DL_FUNC) &_IsoSpecR_RWassersteinMatch, 6},
    {NULL, NULL, 0}
};

void R_init_IsoSpecR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}


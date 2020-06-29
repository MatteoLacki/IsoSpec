%load_ext autoreload
%autoreload 2
import pandas as pd

from IsoSpecPy.IsoSpecPy import IsoParamsFromDict, IsoThreshold, IsoTotalProb, Iso, ParseFormula
from IsoSpecPy.PeriodicTbl import symbol_to_massNo


formula = 'C100H200'
# peaks that make 99.9% of the isotopic distribution

def iso2df(iso, formula, get_confs=False, **kwds):
    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError(e.msg + "\nThis requires pandas to be installed.")
    element_names = list(ParseFormula(formula))
    symbol2isotopeNames = {el: tuple([str(el)+str(int(cnt)) for cnt in cnts]) for el, cnts in symbol_to_massNo.items() }

    df = pd.DataFrame({'mass': iso.np_masses(),
                       'prob': iso.np_probs()})
    if get_confs:
        confs_df = pd.DataFrame.from_records(sum(r, ()) for r in iso.confs)
        confs_df.columns = sum([symbol2isotopeNames[el] for el in element_names], ())
        df = pd.concat([df, confs_df], axis=1)
    return df

def optimalPsetDF(P, formula, get_confs=False,  **kwds):
    iso = IsoTotalProb(P, formula, get_confs=get_confs, **kwds)
    return iso2df(iso, formula, get_confs=get_confs)
    
def thresholdSetDF(thr, formula, get_confs=False, **kwds):
    iso = IsoThreshold(thr, formula, get_confs=get_confs, **kwds)
    return iso2df(iso, formula, get_confs=get_confs)


optimalPsetDF(.999, 'C100H202')
optimalPsetDF(.999, 'C100H202', get_confs=True)

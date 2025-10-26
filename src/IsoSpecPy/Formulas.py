from .IsoSpecPy import Iso
from collections import Counter



loratadine = "C22H23N2O2Cl1"
desloratadine = "C19H19N2Cl1"
titin = "C169719H270466N45688O52238S911"
bovine_insulin = "C254H377N65O75S6"
ubiquitin = "C378H629N105O118S1"
substance_p = "C63H98N18O13S1"
protein_p16 = "C681H1100N216O208S5"
cholesterol = "C27H46O1"
caffeine = "C8H10N4O2"
water = "H2O1"
glucose = "C6H12O6"
sucrose = "C12H22O11"
horse_myoglobin = "C769H1212N210O218S2"
oxygen = "O2"





averagine_unit = Counter({
'C' : 4.9384,
'H' : 7.7583,
'N' : 1.3577,
'O' : 1.4773,
'S' : 0.0417})

averagine_mass = 111.1254

def fillHs(d, target_avg_mass):
    while True:
        i = Iso(formula = d)
        if i.getTheoreticalAverageMass() > target_avg_mass:
            break
        d['H'] += 1

def averagine_dct(target_avg_mass):
    units = float(target_avg_mass) / averagine_mass
    d = dict((key, int(units*val)) for key, val in averagine_unit.items())
    fillHs(d, target_avg_mass)
    for x in list(d.keys()):
        if d[x] == 0:
            del d[x]
    return d

def averagine(target_avg_mass):
    return ''.join(str(key)+str(val) for key, val in sorted(averagine_dct(target_avg_mass).items(), key= lambda x: x[0]))


acgt_avg_base = Counter({
'C' : 9.75,
'H' : 12.25,
'N' : 3.75,
'O' : 6.0,
'P' : 1.0
})

def acgt_dct(no_bases, add_water = True):
    c = {key : int(val*no_bases+0.5) for key, val in acgt_avg_base.items()}
    if add_water:
        c['H'] += 2
        c['O'] += 1
    return c

def acgt(no_bases, add_water = True):
    c = acgt_dct(no_bases, add_water = add_water)
    return ''.join(str(key)+str(val) for key, val in sorted(c.items(), key= lambda x: x[0]))

'''This file contains a set of functions useful for the analysis of properties of
isotopic envelopes. These are unlikely to be of any interest to end-users, only
for developers interested in testing some properties of the isotopic envelopes.
These are implemented mostly in Python, and so, not as fast as the rest of the
package.'''

def _int_neighbours(subisotopologue_conf):
    conf = list(subisotopologue_conf)
    for i in range(len(conf)):
        for j in range(len(conf)):
            if i != j and conf[i] > 0:
                conf[i] -= 1
                conf[j] += 1
                yield tuple(conf)
                conf[i] += 1
                conf[j] -= 1


def neighbours(conf):
    for i in range(len(conf)):
        for si_neigh in _int_neighbours(conf[i]):
            yield tuple(conf[j] if j != i else si_neigh for j in range(len(conf)))




def conf_mass(iso, conf):
    return sum(sum(mass * cnt for mass, cnt in zip(masses, cnts)) for masses, cnts in zip(iso.isotopeMasses, conf))


from IsoSpecPy import Iso
Iso.conf_mass = conf_mass

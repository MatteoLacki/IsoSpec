from __future__ import print_function
import IsoSpecPy
from IsoSpecPy.Formulas import *
import math


try:
    math.isclose
except AttributeError:

    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

    math.isclose = isclose

glu = IsoSpecPy.IsoThreshold(0.0, formula=glucose)
ca = IsoSpecPy.IsoThreshold(0.0, formula=caffeine)


def test_wasserstein_distance():
    print("Checking Wasserstein distance...", end=" ")
    print(ca.wassersteinDistance(glu), end=" ")
    assert math.isclose(ca.wassersteinDistance(glu), 14.03495145836358)
    print("OK!")


def test_normalization():
    print("Checking normalization... ", end="")

    ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
    print(ubiq.total_prob(), end=" ")
    assert math.isclose(ubiq.total_prob(), 0.9999, rel_tol=0.01)
    ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
    ubiq.scale(0.5)
    assert math.isclose(ubiq.total_prob(), 0.9999 * 0.5, rel_tol=0.01)
    ubiq._recalculate_everything()
    assert math.isclose(ubiq.total_prob(), 0.9999 * 0.5, rel_tol=0.01)
    ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
    ubiq.scale(0.5)
    ubiq.normalize()
    assert math.isclose(ubiq.total_prob(), 1.0)
    ubiq._recalculate_everything()
    assert math.isclose(ubiq.total_prob(), 1.0)
    print("OK!")


def test_addition():
    print("Checking addition...", end=" ")
    wa = IsoSpecPy.IsoThreshold(0.0, formula=water)
    ox = IsoSpecPy.IsoThreshold(0.0, formula=oxygen)
    s = wa + ox
    assert math.isclose(s.total_prob(), 2.0)
    assert len(list(s)) == len(list(wa)) + len(list(ox)) == 15
    print("OK!")


def test_sorting():
    print("Checking sorting...", end=" ")
    ubiq = IsoSpecPy.IsoTotalProb(0.9999, ubiquitin)
    ubiq.sort_by_mass()
    assert list(ubiq.masses) == sorted(ubiq.masses)
    ubiq.sort_by_prob()
    assert list(ubiq.probs) == sorted(ubiq.probs)
    print("OK!")


def test_binning():
    print("Checking binning...", end=" ")
    ubiq = IsoSpecPy.IsoTotalProb(0.999999, ubiquitin)
    # ubiq.plot()
    print(len(ubiq), end=" -> ")
    bu = ubiq.binned()
    ubiq._recalculate_everything()
    print(len(bu), end=" ")
    bu._recalculate_everything()
    assert math.isclose(ubiq.total_prob(), bu.total_prob())
    print("OK!")


def test_convolution():
    print("Checking convolution...", end=" ")
    o = IsoSpecPy.IsoThreshold(0.0, formula="H1")
    sur = IsoSpecPy.IsoThreshold(0.0, formula=sucrose)
    wa = IsoSpecPy.IsoThreshold(0.0, formula=water)
    # (glu*glu).plot()
    # (sur*wa).plot()
    WSD = (sur * wa).wassersteinDistance(glu * glu)
    print(WSD, end=" ")
    assert math.isclose(WSD, 0.0, abs_tol=1e-7)
    print("OK!")


def test_negative_formulas():
    print("Checking negative formulas... ", end="")
    try:
        I = Iso(formula="C-10")
        print("FAIL: exception not thrown")
        raise AssertionError("Exception not thrown for negative formula")
    except Exception as e:
        print(
            """exception successfully obtained, message: "{}" -> OK!""".format((str(e)))
        )


def test_fasta_negative_formulas():
    print("Checking FASTA + negative formulas... ", end="")
    try:
        I = Iso(fasta="C", formula="C-5")
        print("FAIL: exception not thrown")
    except Exception as e:
        print(
            """exception successfully obtained, message: "{}" -> OK!""".format((str(e)))
        )


def test_fasta_modification():
    print("Checking FASTA + modification... ", end="")
    # Selenation of methionine
    I = IsoSpecPy.IsoTotalProb(0.999, formula="C5H9N1O1Se1")
    I2 = IsoSpecPy.IsoTotalProb(0.999, fasta="M", formula="S-1Se1")
    WSD = I.wassersteinDistance(I2)
    print(WSD, end="")
    assert math.isclose(I.wassersteinDistance(I2), 0.0)
    print(" -> OK!")


formulas = "C1 P1 H1 H100 P100 N100 O100 C100H100N100 C100O100".split()


def test_empiric_avg_mass():
    print("Checking empiric avg mass... ", end="")
    isos = [IsoSpecPy.Iso(formula) for formula in formulas]
    dists = [IsoSpecPy.IsoThreshold(0.0, formula) for formula in formulas]
    diffs = [
        abs(i.getTheoreticalAverageMass() - d.empiric_average_mass())
        for i, d in zip(isos, dists)
    ]
    print(max(diffs), end="")
    assert max(diffs) < 1e-6
    print(" -> OK!")


def test_empiric_variance():
    print("Checking empiric variance... ", end="")
    isos = [IsoSpecPy.Iso(formula) for formula in formulas]
    dists = [IsoSpecPy.IsoThreshold(0.0, formula) for formula in formulas]
    diffs = [abs(i.variance() - d.empiric_variance()) for i, d in zip(isos, dists)]
    print(max(diffs), end="")
    assert max(diffs) < 1e-6
    print(" -> OK!")


def test_empiric_stddev():
    print("Checking empiric stddev... ", end="")
    isos = [IsoSpecPy.Iso(formula) for formula in formulas]
    dists = [IsoSpecPy.IsoThreshold(0.0, formula) for formula in formulas]
    diffs = [abs(i.stddev() - d.empiric_stddev()) for i, d in zip(isos, dists)]
    print(max(diffs), end="")
    assert max(diffs) < 1e-6
    print(" -> OK!")

def test_get_monoisotopic_mass():
    print("Checking monoisotopic mass...", end=" ")
    mol = IsoSpecPy.Iso(formula="C100H100N100O100Se100Sn100Pb100U100")
    mono_mass = mol.getMonoisotopicPeakMass()
    print(mono_mass, end=" ")
    assert math.isclose(68885.198515667, mono_mass, rel_tol=1e-9)
    print("OK!")

def test_lightest_peak():
    print("Checking lightest peak...", end=" ")
    formula = "C10B10H10Sn1"
    iso = IsoSpecPy.Iso(formula=formula)
    lightest_mass = iso.getLightestPeakMass()
    lightest_lprob = iso.getLightestPeakLProb()
    lightest_conf = iso.getLightestPeakConf()
    iso_threshold = IsoSpecPy.IsoThreshold(0.0, formula=formula, get_confs=True)
    masses = list(iso_threshold.masses)
    probs = list(iso_threshold.probs)
    confs = list(iso_threshold.confs)
    min_index = masses.index(min(masses))
    print(lightest_conf, lightest_mass, lightest_lprob, end=" ")
    assert lightest_conf == confs[min_index]
    assert lightest_mass == masses[min_index]
    assert math.isclose(lightest_lprob, math.log(probs[min_index]), rel_tol=1e-9)
    print("OK!")

def test_heaviest_peak():
    print("Checking heaviest peak...", end=" ")
    formula = "C10B10H10Sn1"
    iso = IsoSpecPy.Iso(formula=formula)
    heaviest_mass = iso.getHeaviestPeakMass()
    heaviest_lprob = iso.getHeaviestPeakLProb()
    heaviest_conf = iso.getHeaviestPeakConf()
    iso_threshold = IsoSpecPy.IsoThreshold(0.0, formula=formula, get_confs=True)
    masses = list(iso_threshold.masses)
    probs = list(iso_threshold.probs)
    confs = list(iso_threshold.confs)
    max_index = masses.index(max(masses))
    print(heaviest_conf, heaviest_mass, heaviest_lprob, end=" ")
    assert heaviest_conf == confs[max_index]
    assert heaviest_mass == masses[max_index]
    assert math.isclose(heaviest_lprob, math.log(probs[max_index]), rel_tol=1e-9)
    print("OK!")

def test_monoisotopic_peak():
    print("Checking monoisotopic peak...", end=" ")
    formula = "C10B1000H10Sn1"
    iso = IsoSpecPy.Iso(formula=formula)
    monoisotopic_mass = iso.getMonoisotopicPeakMass()
    monoisotopic_lprob = iso.getMonoisotopicPeakLProb()
    monoisotopic_conf = iso.getMonoisotopicPeakConf()
    iso_threshold = IsoSpecPy.IsoThreshold(0.0, formula=formula, get_confs=True)
    masses = list(iso_threshold.masses)
    probs = list(iso_threshold.probs)
    confs = list(iso_threshold.confs)
    monoisotopic_conf2 = ((10, 0), (0, 1000), (10, 0), (0,0,0,0,0,0,0,1,0,0))
    monoisotopic_peak_idx = confs.index(monoisotopic_conf)
    print(monoisotopic_conf, monoisotopic_mass, monoisotopic_lprob, end=" ")
    assert monoisotopic_conf == confs[monoisotopic_peak_idx]
    assert monoisotopic_conf == monoisotopic_conf2
    assert monoisotopic_mass == masses[monoisotopic_peak_idx]
    assert math.isclose(monoisotopic_lprob, math.log(probs[monoisotopic_peak_idx]), rel_tol=1e-9)
    print("OK!")

if __name__ == "__main__":
    test_wasserstein_distance()
    test_normalization()
    test_addition()
    test_sorting()
    test_binning()
    test_convolution()
    test_negative_formulas()
    test_fasta_negative_formulas()
    test_fasta_modification()
    test_empiric_avg_mass()
    test_empiric_variance()
    test_empiric_stddev()
    test_get_monoisotopic_mass()
    test_lightest_peak()
    test_heaviest_peak()
    test_monoisotopic_peak()

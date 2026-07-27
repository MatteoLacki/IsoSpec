/* The C ABI, compiled as C.
 *
 * cwrapper.h is the header non-C++ callers use.  Everything else in the suite
 * exercises it from C++ (where a C++-only construct leaking into the header
 * would go unnoticed), so this translation unit is deliberately compiled with a
 * C compiler in C11 mode with -Wall -Wextra -Werror: it fails the build if the
 * header ever stops being valid C, and fails at runtime if the ABI misbehaves.
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cwrapper.h"

static int failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);           \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

#define CHECK_CLOSE(a, b, tol)                                               \
    do {                                                                     \
        const double va = (a), vb = (b);                                     \
        if (!(fabs(va - vb) <= (tol) * (1.0 + fabs(vb)))) {                  \
            printf("FAIL %s:%d: %s == %s (%.17g vs %.17g)\n",                \
                   __FILE__, __LINE__, #a, #b, va, vb);                      \
            ++failures;                                                      \
        }                                                                    \
    } while (0)

/* Water: two hydrogens (2 isotopes) and one oxygen (3 isotopes). */
static const int isotope_numbers[2] = {2, 3};
static const int atom_counts[2] = {2, 1};
static const double isotope_masses[5] = {
    1.00782503207, 2.0141017778,                       /* H  */
    15.99491461956, 16.99913170, 17.9991610            /* O  */
};
static const double isotope_probabilities[5] = {
    0.999885, 0.000115,
    0.9975716, 0.0003804, 0.002048
};

static void* water_iso(void)
{
    return setupIso(2, isotope_numbers, atom_counts,
                    isotope_masses, isotope_probabilities);
}

static void test_iso_accessors(void)
{
    void* iso = water_iso();
    CHECK(iso != NULL);
    if (iso == NULL) return;

    const double monoisotopic = 2 * isotope_masses[0] + isotope_masses[2];
    CHECK_CLOSE(getMonoisotopicPeakMassIso(iso), monoisotopic, 1e-12);
    CHECK_CLOSE(getLightestPeakMassIso(iso), monoisotopic, 1e-12);
    CHECK_CLOSE(getHeaviestPeakMassIso(iso),
                2 * isotope_masses[1] + isotope_masses[4], 1e-12);

    const double average = 2 * (isotope_masses[0] * isotope_probabilities[0] +
                                isotope_masses[1] * isotope_probabilities[1]) +
                           isotope_masses[2] * isotope_probabilities[2] +
                           isotope_masses[3] * isotope_probabilities[3] +
                           isotope_masses[4] * isotope_probabilities[4];
    CHECK_CLOSE(getTheoreticalAverageMassIso(iso), average, 1e-12);
    CHECK(getIsoVariance(iso) > 0.0);
    CHECK_CLOSE(getIsoStddev(iso), sqrt(getIsoVariance(iso)), 1e-12);
    CHECK(getModeLProbIso(iso) <= 0.0);

    int signature[5];
    getMonoisotopicPeakSignature(iso, signature);
    CHECK(signature[0] == 2);
    CHECK(signature[1] == 0);
    CHECK(signature[2] == 1);
    CHECK(signature[3] == 0);
    CHECK(signature[4] == 0);

    double* estimates = getMarginalLogSizeEstimates(iso, 0.99);
    CHECK(estimates != NULL);
    freeReleasedArray(estimates);

    deleteIso(iso);
}

static void test_error_paths(void)
{
    /* A zero-probability isotope must come back as NULL, not as an exception
     * unwinding into this C frame. */
    const double bad_probs[5] = {0.0, 1.0, 1.0, 0.0, 0.0};
    void* bad = setupIso(2, isotope_numbers, atom_counts,
                         isotope_masses, bad_probs);
    CHECK(bad == NULL);
    deleteIso(NULL);  /* must tolerate NULL */

    /* A zero-dimensional Iso builds, but no generator accepts it. */
    void* empty = setupIso(0, NULL, NULL, NULL, NULL);
    CHECK(empty != NULL);
    void* generator = setupIsoThresholdGenerator(empty, 0.01, true, 1000, 1000, true);
    CHECK(generator == NULL);
    deleteIsoThresholdGenerator(NULL);
    deleteIso(empty);
}

static void test_threshold_generator(void)
{
    void* iso = water_iso();
    void* generator = setupIsoThresholdGenerator(iso, 1e-12, true, 1000, 1000, true);
    CHECK(generator != NULL);
    if (generator == NULL) { deleteIso(iso); return; }

    size_t peaks = 0;
    double total = 0.0;
    int signature[5];
    while (advanceToNextConfigurationIsoThresholdGenerator(generator)) {
        const double prob = probIsoThresholdGenerator(generator);
        const double mass = massIsoThresholdGenerator(generator);
        CHECK(prob >= 1e-12);
        CHECK(mass > 0.0);
        CHECK_CLOSE(prob, exp(lprobIsoThresholdGenerator(generator)), 1e-9);
        get_conf_signatureIsoThresholdGenerator(generator, signature);
        CHECK(signature[0] + signature[1] == 2);
        CHECK(signature[2] + signature[3] + signature[4] == 1);
        total += prob;
        ++peaks;
    }
    CHECK(peaks == 9);  /* 3 hydrogen configurations x 3 oxygen isotopes */
    CHECK_CLOSE(total, 1.0, 1e-9);

    deleteIsoThresholdGenerator(generator);
    deleteIso(iso);
}

static void test_other_generators(void)
{
    void* iso = water_iso();
    void* layered = setupIsoLayeredGenerator(iso, 1000, 1000, true, 0.99);
    CHECK(layered != NULL);
    double total = 0.0;
    while (advanceToNextConfigurationIsoLayeredGenerator(layered))
        total += probIsoLayeredGenerator(layered);
    CHECK_CLOSE(total, 1.0, 1e-9);
    deleteIsoLayeredGenerator(layered);
    deleteIso(iso);

    iso = water_iso();
    void* ordered = setupIsoOrderedGenerator(iso, 1000, 1000);
    CHECK(ordered != NULL);
    double previous = 2.0;
    total = 0.0;
    while (advanceToNextConfigurationIsoOrderedGenerator(ordered)) {
        const double prob = probIsoOrderedGenerator(ordered);
        CHECK(prob <= previous + 1e-15);
        previous = prob;
        total += prob;
    }
    CHECK_CLOSE(total, 1.0, 1e-9);
    deleteIsoOrderedGenerator(ordered);
    deleteIso(iso);

    iso = water_iso();
    void* stochastic = setupIsoStochasticGenerator(iso, 10000, 0.9999, 5.0);
    CHECK(stochastic != NULL);
    total = 0.0;
    while (advanceToNextConfigurationIsoStochasticGenerator(stochastic))
        total += probIsoStochasticGenerator(stochastic);
    CHECK_CLOSE(total, 10000.0, 1e-9);
    deleteIsoStochasticGenerator(stochastic);
    deleteIso(iso);
}

static void test_envelopes(void)
{
    void* iso = water_iso();
    void* envelope = setupThresholdFixedEnvelope(iso, 1e-12, true, false);
    CHECK(envelope != NULL);
    if (envelope == NULL) { deleteIso(iso); return; }

    const size_t size = confs_noFixedEnvelope(envelope);
    CHECK(size == 9);
    CHECK_CLOSE(getTotalProbOfEnvelope(envelope), 1.0, 1e-9);
    CHECK(empiricAverageMass(envelope) > 0.0);
    CHECK_CLOSE(empiricStddev(envelope), sqrt(empiricVariance(envelope)), 1e-9);

    sortEnvelopeByMass(envelope);
    /* massesFixedEnvelope releases the array to the caller. */
    const double* masses = massesFixedEnvelope(envelope);
    const double* probs = probsFixedEnvelope(envelope);
    CHECK(masses != NULL);
    CHECK(probs != NULL);
    if (masses != NULL && probs != NULL) {
        double sum = 0.0;
        for (size_t i = 0; i < size; ++i) {
            if (i > 0) CHECK(masses[i - 1] <= masses[i]);
            sum += probs[i];
        }
        CHECK_CLOSE(sum, 1.0, 1e-9);
    }
    freeReleasedArray((void*)masses);
    freeReleasedArray((void*)probs);
    deleteFixedEnvelope(envelope, false);
    deleteIso(iso);

    /* Envelopes over caller-owned arrays, and the arithmetic on them. */
    double m1[2] = {1.0, 2.0};
    double p1[2] = {0.5, 0.5};
    double m2[2] = {10.0, 20.0};
    double p2[2] = {0.25, 0.75};
    void* a = setupFixedEnvelope(m1, p1, 2, false, false, NAN);
    void* b = setupFixedEnvelope(m2, p2, 2, false, false, NAN);
    CHECK(a != NULL);
    CHECK(b != NULL);

    void* sum = addEnvelopes(a, b);
    CHECK(confs_noFixedEnvelope(sum) == 4);
    CHECK_CLOSE(getTotalProbOfEnvelope(sum), 2.0, 1e-12);

    void* product = convolveEnvelopes(a, b);
    CHECK(confs_noFixedEnvelope(product) == 4);
    CHECK_CLOSE(getTotalProbOfEnvelope(product), 1.0, 1e-12);
    CHECK_CLOSE(empiricAverageMass(product),
                empiricAverageMass(a) + empiricAverageMass(b), 1e-12);

    void* binned = binnedEnvelope(product, 1.0, 0.0);
    CHECK_CLOSE(getTotalProbOfEnvelope(binned), 1.0, 1e-12);

    void* const spectra[2] = {a, b};
    const double intensities[2] = {2.0, 3.0};
    void* combination = linearCombination(spectra, intensities, 2);
    CHECK(confs_noFixedEnvelope(combination) == 4);
    CHECK_CLOSE(getTotalProbOfEnvelope(combination), 5.0, 1e-12);

    void* duplicate = copyFixedEnvelope(product);
    CHECK(confs_noFixedEnvelope(duplicate) == confs_noFixedEnvelope(product));

    /* Distances. */
    CHECK_CLOSE(wassersteinDistance(a, a), 0.0, 1e-12);
    CHECK_CLOSE(abyssalWassersteinDistance(a, a, 1.0, 1.0), 0.0, 1e-12);
    CHECK_CLOSE(orientedWassersteinDistance(a, b),
                -orientedWassersteinDistance(b, a), 1e-9);
    struct ws_match_res match = wassersteinMatch(a, a, 0.5, 1.0);
    CHECK_CLOSE(match.flow, 1.0, 1e-9);
    CHECK_CLOSE(match.res1, 0.0, 1e-9);
    CHECK_CLOSE(match.res2, 0.0, 1e-9);

    /* An unnormalized pair must give NaN, not an exception. */
    double m3[1] = {0.0};
    double p3[1] = {17.0};
    void* unnormalized = setupFixedEnvelope(m3, p3, 1, false, false, NAN);
    CHECK(isnan(wassersteinDistance(a, unnormalized)));

    deleteFixedEnvelope(unnormalized, true);
    deleteFixedEnvelope(duplicate, false);
    deleteFixedEnvelope(combination, false);
    deleteFixedEnvelope(binned, false);
    deleteFixedEnvelope(product, false);
    deleteFixedEnvelope(sum, false);
    deleteFixedEnvelope(b, true);
    deleteFixedEnvelope(a, true);
}

static void test_fasta_and_arrays(void)
{
    void* iso = isoFromFasta("PEPTIDE", false, true);
    CHECK(iso != NULL);
    CHECK(getMonoisotopicPeakMassIso(iso) > 799.0);
    CHECK(getMonoisotopicPeakMassIso(iso) < 800.0);
    deleteIso(iso);

    int counts[6];
    parse_fasta_c("AAA", counts);
    CHECK(counts[0] == 9);   /* C  */
    CHECK(counts[1] == 15);  /* H  */
    CHECK(counts[2] == 3);   /* N  */
    CHECK(counts[3] == 3);   /* O  */
    CHECK(counts[4] == 0);   /* S  */
    CHECK(counts[5] == 0);   /* Se */

    double values[4] = {1.0, 2.0, 3.0, 4.0};
    array_add(values, 4, 10.0);
    CHECK_CLOSE(values[0], 11.0, 1e-12);
    array_mul(values, 4, 2.0);
    CHECK_CLOSE(values[0], 22.0, 1e-12);
    array_fma(values, 4, 0.5, 1.0);
    CHECK_CLOSE(values[0], 12.0, 1e-12);
}

int main(void)
{
    test_iso_accessors();
    test_error_paths();
    test_threshold_generator();
    test_other_generators();
    test_envelopes();
    test_fasta_and_arrays();

    if (failures == 0) {
        printf("C API tests: all checks passed\n");
        return 0;
    }
    printf("C API tests: %d check(s) FAILED\n", failures);
    return 1;
}

// Regression tests for bugs identified during code review, converted to
// doctest.  Each case pins a specific contract so a fix cannot silently
// regress.  Several are only conclusive under ASan/UBSan (built as a separate
// configuration in the Makefile) — they exercise latent UB that a plain build
// rolls over.

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "doctest.h"
#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "pod_vector.h"
#include "fasta.h"

using namespace IsoSpec;

// Bug 1: unsafe_pod_vector<T>::swap must take unsafe_pod_vector<T>&, not
// pod_vector<T>&.
TEST_CASE("unsafe_pod_vector::swap has the correct signature and swaps") {
    using expected_t =
        void (unsafe_pod_vector<int>::*)(unsafe_pod_vector<int>&) noexcept;
    static_assert(std::is_same_v<decltype(&unsafe_pod_vector<int>::swap), expected_t>,
                  "unsafe_pod_vector<T>::swap must take unsafe_pod_vector<T>&");

    unsafe_pod_vector<int> a, b;
    a.init(4);
    b.init(4);
    a.push_back(1);
    a.push_back(2);
    b.push_back(99);

    a.swap(b);

    CHECK(a.size() == 1);
    CHECK(a[0] == 99);
    CHECK(b.size() == 2);
    CHECK(b[0] == 1);
    CHECK(b[1] == 2);

    a.clear();
    b.clear();
}

// Bug 2: FixedEnvelope::Binned() must not wrap a negative bin index to a huge
// size_t.  Conclusive under UBSan/ASan.
TEST_CASE("Binned() with large bin_middle does not wrap the index") {
    Iso iso("H2O1", false);
    const double bin_width = 0.01;
    const double bin_middle = 1.0e6;  // >> H2O monoisotopic mass -> negative offset

    FixedEnvelope env;
    CHECK_NOTHROW(env = FixedEnvelope::Binned(std::move(iso), 0.99, bin_width, bin_middle));
    CHECK(env.confs_no() < 1000000);  // absurd size would indicate a wrap
}

// Bug 3: fasta.h must be self-contained (uses memset/size_t).  The compile-only
// isolation check is fasta_header_isolation.cpp; here we confirm parse_fasta
// actually counts correctly.
TEST_CASE("parse_fasta counts elements correctly") {
    int counts[6] = {-1, -1, -1, -1, -1, -1};
    parse_fasta("AAA", counts);
    CHECK(counts[0] == 9);   // C
    CHECK(counts[1] == 15);  // H
    CHECK(counts[2] == 3);   // N
    CHECK(counts[3] == 3);   // O
    CHECK(counts[4] == 0);   // S
    CHECK(counts[5] == 0);   // Se
}

// Bug 4: constructing a generator from an empty (default-ctor) Iso must not
// underflow `new double[dimNumber-1]` into a SIZE_MAX allocation.  Acceptable
// outcomes: clean invalid_argument, or success with zero configurations.
namespace {
template <typename Generator, typename... Args>
void check_empty_iso_rejected(Args&&... args) {
    Iso iso;  // dimNumber == 0
    try {
        Generator gen(std::move(iso), std::forward<Args>(args)...);
        CHECK_FALSE(gen.advanceToNextConfiguration());
    } catch (const std::invalid_argument&) {
        // Preferred: clean rejection.
    } catch (const std::bad_alloc&) {
        FAIL("empty Iso threw bad_alloc (SIZE_MAX allocation from dimNumber-1)");
    } catch (const std::exception& e) {
        FAIL("empty Iso threw unexpected exception: " << e.what());
    }
}
}  // namespace

TEST_CASE("empty Iso is rejected cleanly by every generator") {
    check_empty_iso_rejected<IsoThresholdGenerator>(0.01, false);
    check_empty_iso_rejected<IsoLayeredGenerator>();
    check_empty_iso_rejected<IsoOrderedGenerator>();
}

// Bug 5: extern "C" wrappers must not let a C++ exception cross the ABI.
extern "C" void* setupIso(int dimNumber, const int* isotopeNumbers,
                          const int* atomCounts, const double* isotopeMasses,
                          const double* isotopeProbabilities);
extern "C" void deleteIso(void* iso);

TEST_CASE("setupIso rejects invalid input with nullptr, not an exception") {
    const int isoNumbers[1] = {2};
    const int atomCounts[1] = {10};
    const double masses[2] = {1.0, 2.0};
    const double probs[2] = {0.0, 1.0};  // 0.0 rejected by Marginal

    void* p = reinterpret_cast<void*>(-1);
    CHECK_NOTHROW(p = setupIso(1, isoNumbers, atomCounts, masses, probs));
    CHECK(p == nullptr);
    deleteIso(p);  // must tolerate nullptr

    const double good_probs[2] = {0.5, 0.5};
    void* q = setupIso(1, isoNumbers, atomCounts, masses, good_probs);
    CHECK(q != nullptr);
    deleteIso(q);
}

// Bug 6: pod_vector with initial capacity 0 must accept push_back / resize.
TEST_CASE("zero-capacity pod_vector accepts push_back and resize") {
    pod_vector<int> v(0);
    CHECK_NOTHROW(v.push_back(7));
    CHECK(v.size() == 1);
    CHECK(v[0] == 7);

    pod_vector<int> v2(0);
    CHECK_NOTHROW(v2.resize(16));
    CHECK(v2.size() == 16);
}

// Bug 7: the const-Iso& overloads of FixedEnvelope's named constructors used to
// take a *shallow* copy (Iso(iso, false)), which shares the marginals with the
// caller's Iso.  The generator built from that copy move-constructs its own
// marginals out of them and frees the underlying arrays when it finishes,
// leaving the caller holding dangling pointers: every subsequent accessor was a
// use-after-free.  Conclusive under ASan.
TEST_CASE("named constructors taking const Iso& leave the Iso usable") {
    Iso iso("C10H20O2");
    const double avg = iso.getTheoreticalAverageMass();
    const double var = iso.variance();
    const double mode = iso.getModeLProb();

    {
        FixedEnvelope a = FixedEnvelope::FromThreshold(iso, 1e-6, true, false);
        CHECK(a.confs_no() > 0);
    }
    CHECK(iso.getTheoreticalAverageMass() == avg);

    {
        FixedEnvelope b = FixedEnvelope::FromTotalProb(iso, 0.99, true, false);
        CHECK(b.confs_no() > 0);
    }
    CHECK(iso.variance() == var);

    {
        FixedEnvelope c = FixedEnvelope::FromStochastic(iso, 1000);
        CHECK(c.confs_no() > 0);
    }
    CHECK(iso.getModeLProb() == mode);

    {
        FixedEnvelope d = FixedEnvelope::Binned(iso, 0.99, 1.0);
        CHECK(d.confs_no() > 0);
    }
    CHECK(iso.getTheoreticalAverageMass() == avg);
    CHECK(iso.variance() == var);

    // ... and the Iso can still be consumed by a generator afterwards.
    IsoThresholdGenerator gen(std::move(iso), 1e-6, true);
    std::size_t n = 0;
    while (gen.advanceToNextConfiguration()) ++n;
    CHECK(n > 0);
}

// Bug 8: OrientedWassersteinDistance accumulated the signed CDF difference, but
// its "this"-tail loop *subtracted* the remaining peaks (a copy-paste from the
// unsigned version, which folds the sign into abs() first).  The result was
// wrong — and not even antisymmetric — whenever the other spectrum ran out
// first.
TEST_CASE("OrientedWassersteinDistance is antisymmetric") {
    double m1[2] = {0.0, 1.0}, p1[2] = {0.5, 0.5};
    double m2[2] = {2.0, 3.0}, p2[2] = {0.5, 0.5};
    FixedEnvelope a(m1, p1, 2);
    FixedEnvelope b(m2, p2, 2);

    const double ab = a.OrientedWassersteinDistance(b);
    const double ba = b.OrientedWassersteinDistance(a);
    CHECK(ab == doctest::Approx(-ba));
    // Everything moves one way, so the magnitude is the plain distance.
    CHECK(std::fabs(ab) == doctest::Approx(a.WassersteinDistance(b)));

    a.release_everything();
    b.release_everything();
}

// Bug 9: AbyssalWassersteinDistance updated a *copy* of the carried peak's
// remaining mass when an incoming peak was only partially matched, so the
// matched mass stayed on the carried list and was counted twice — once as
// transport, once as condemned mass.  Identical spectra came out at
// total_prob * abyss_depth / 2 instead of 0.
TEST_CASE("AbyssalWassersteinDistance of a spectrum to itself is zero") {
    double m1[2] = {0.0, 1.0}, p1[2] = {0.5, 0.5};
    double m2[2] = {0.0, 1.0}, p2[2] = {0.5, 0.5};
    FixedEnvelope a(m1, p1, 2);
    FixedEnvelope b(m2, p2, 2);

    for (double depth : {0.1, 1.0, 10.0}) {
        CAPTURE(depth);
        CHECK(a.AbyssalWassersteinDistance(b, depth) == doctest::Approx(0.0).epsilon(1e-12));
    }

    a.release_everything();
    b.release_everything();
}

// Bug 10 (documentation only): the two-step init lists in Iso(int,...) and
// Marginal::Marginal(...) leak the first-allocated member if the second
// allocation throws.  Catchable only by injecting an allocation failure at the
// right point (custom operator new keyed on allocation count) under
// -fsanitize=leak.  Recorded as a message rather than a fake passing test.
TEST_CASE("DOC: ctor-leak in Iso(int,...) / Marginal on second-alloc failure") {
    MESSAGE("Iso(int,...) and Marginal ctors leak the first array if the second "
            "allocation throws; test by gating the Nth operator new[] to throw "
            "and running under -fsanitize=leak. No automated check yet.");
}

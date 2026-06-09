// Regression tests for bugs identified during code review.
//
// Build (plain):
//   clang++ -std=c++17 -O0 -g -I../../src/IsoSpec++ \
//       ../../src/IsoSpec++/unity-build.cpp \
//       bug_regression_tests.cpp -o bug_regression_tests
//
// Build (with sanitisers — recommended; some bugs only surface here):
//   clang++ -std=c++17 -O0 -g -fsanitize=address,undefined \
//       -I../../src/IsoSpec++ \
//       ../../src/IsoSpec++/unity-build.cpp \
//       bug_regression_tests.cpp -o bug_regression_tests_asan
//
// On a clean tree these tests are expected to FAIL (or produce sanitiser
// reports).  After each underlying bug is fixed, the corresponding test
// must pass.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "isoSpec++.h"
#include "fixedEnvelopes.h"
#include "pod_vector.h"
#include "fasta.h"

using namespace IsoSpec;

static int g_failures = 0;

#define REQUIRE(cond, msg) do { \
    if (!(cond)) { \
        std::fprintf(stderr, "FAIL [%s]: %s  (%s:%d)\n", \
                     __func__, msg, __FILE__, __LINE__); \
        ++g_failures; \
    } \
} while (0)

#define EXPECT_NO_THROW(stmt, msg) do { \
    try { stmt; } \
    catch (const std::exception& e) { \
        std::fprintf(stderr, "FAIL [%s]: %s threw '%s'  (%s:%d)\n", \
                     __func__, msg, e.what(), __FILE__, __LINE__); \
        ++g_failures; \
    } \
} while (0)


// ===========================================================================
// Bug 1: pod_vector.h:367 — unsafe_pod_vector<T>::swap once took
// pod_vector<T>& instead of unsafe_pod_vector<T>&.  Latent bug — no caller
// existed — but the signature was wrong.
//
// Static check: the corrected signature must be a well-formed
// member-function pointer.  Run-time check: two unsafe_pod_vectors can
// actually be swapped and their contents end up swapped.
// ===========================================================================
static void test_unsafe_pod_vector_swap_signature()
{
    using expected_t = void (unsafe_pod_vector<int>::*)(unsafe_pod_vector<int>&) noexcept;
    static_assert(std::is_same_v<decltype(&unsafe_pod_vector<int>::swap), expected_t>,
                  "unsafe_pod_vector<T>::swap must take unsafe_pod_vector<T>&");

    unsafe_pod_vector<int> a, b;
    a.init(4);
    b.init(4);
    a.push_back(1);
    a.push_back(2);
    b.push_back(99);

    a.swap(b);

    REQUIRE(a.size() == 1 && a[0] == 99,
            "after swap, a should hold b's prior contents");
    REQUIRE(b.size() == 2 && b[0] == 1 && b[1] == 2,
            "after swap, b should hold a's prior contents");

    a.clear();
    b.clear();
}


// ===========================================================================
// Bug 2: FixedEnvelope::Binned() casts a possibly-negative double to
// size_t at fixedEnvelopes.cpp:941 and again at :966.
//
// When (min_mass + half_width - bin_middle) is negative the cast wraps
// to a huge size_t, the subsequent acc -= idx_min underflows the pointer,
// and every subsequent acc[...] indexing is UB.
//
// This test is only conclusive under UBSan + ASan, which turn the latent
// UB into hard errors.  Under plain build the call usually returns a
// nonsense envelope without crashing.
// ===========================================================================
static void test_binned_negative_offset_does_not_wrap()
{
    Iso iso("H2O1", false);
    // Pick bin_middle much larger than the H2O monoisotopic mass so
    // min_mass + half_width - bin_middle is strongly negative.
    const double bin_width  = 0.01;
    const double bin_middle = 1.0e6;

    EXPECT_NO_THROW(
        {
            FixedEnvelope env = FixedEnvelope::Binned(
                std::move(iso), /*target_total_prob=*/0.99,
                bin_width, bin_middle);
            REQUIRE(env.confs_no() < 1000000,
                    "Binned() produced an absurdly large envelope, indicating index wrap");
        },
        "Binned() with bin_middle much larger than min_mass");
}


// ===========================================================================
// Bug 3 (header isolation): fasta.h uses memset() and size_t in its inline
// parse_fasta() but does NOT include <cstring> or <cstddef>.  The
// compile-only test lives in fasta_header_isolation.cpp in this directory;
// here we just confirm parse_fasta is callable from a TU that *does* pull
// the missing headers in transitively.  See that file for the real check.
// ===========================================================================
static void test_fasta_header_is_callable()
{
    int counts[6] = {-1, -1, -1, -1, -1, -1};
    parse_fasta("AAA", counts);
    REQUIRE(counts[0] == 9,  "parse_fasta did not count carbon correctly");
    REQUIRE(counts[1] == 15, "parse_fasta did not count hydrogen correctly");
    REQUIRE(counts[2] == 3,  "parse_fasta did not count nitrogen correctly");
    REQUIRE(counts[3] == 3,  "parse_fasta did not count oxygen correctly");
    REQUIRE(counts[4] == 0,  "parse_fasta did not count sulfur correctly");
    REQUIRE(counts[5] == 0,  "parse_fasta did not count selenium correctly");
}


// ===========================================================================
// Bug 4: Iso() default ctor leaves dimNumber == 0, and constructing a
// generator from such an Iso allocates `new double[dimNumber-1]`, which
// underflows to `new double[SIZE_MAX]` and throws bad_alloc.
//
// Constructing a generator from an empty Iso should either succeed and
// produce no configurations, or throw a *meaningful* exception
// (invalid_argument) — not bad_alloc from a SIZE_MAX-sized allocation.
// ===========================================================================
template<typename Generator, typename... Args>
static void check_empty_iso_rejected(const char* name, Args&&... args)
{
    Iso iso;  // dimNumber == 0
    try {
        Generator gen(std::move(iso), std::forward<Args>(args)...);
        // If construction succeeds, the generator must report no work.
        REQUIRE(!gen.advanceToNextConfiguration(),
                "empty Iso somehow produced configurations");
    } catch (const std::invalid_argument&) {
        // Preferred outcome: clean rejection.
    } catch (const std::bad_alloc&) {
        std::fprintf(stderr,
            "FAIL [empty Iso -> %s]: threw bad_alloc; should reject empty "
            "input with std::invalid_argument before allocating "
            "maxConfsLPSum/dimNumber-1.\n", name);
        ++g_failures;
    } catch (const std::exception& e) {
        std::fprintf(stderr,
            "FAIL [empty Iso -> %s]: threw unexpected exception '%s' "
            "(expected std::invalid_argument or success).\n",
            name, e.what());
        ++g_failures;
    }
}

static void test_empty_iso_generator_rejected_cleanly()
{
    check_empty_iso_rejected<IsoThresholdGenerator>(
        "IsoThresholdGenerator", /*threshold=*/0.01, /*absolute=*/false);
    check_empty_iso_rejected<IsoLayeredGenerator>(
        "IsoLayeredGenerator");
    check_empty_iso_rejected<IsoOrderedGenerator>(
        "IsoOrderedGenerator");
}


// ===========================================================================
// Bug 5 (fixed): extern "C" wrappers in cwrapper.cpp must not let a C++
// exception escape the ABI — that is UB for a C / cffi caller.  setupIso()
// with invalid input (e.g. probability 0, rejected by Marginal) must now
// return nullptr instead of throwing.
// ===========================================================================
extern "C" void* setupIso(int dimNumber,
                          const int* isotopeNumbers,
                          const int* atomCounts,
                          const double* isotopeMasses,
                          const double* isotopeProbabilities);
extern "C" void deleteIso(void* iso);

static void test_cwrapper_rejects_invalid_input_with_nullptr()
{
    const int isoNumbers[1] = {2};
    const int atomCounts[1] = {10};
    const double masses[2]  = {1.0, 2.0};
    const double probs[2]   = {0.0, 1.0};  // 0.0 is rejected by Marginal

    bool escaped = false;
    void* p = reinterpret_cast<void*>(-1);  // sentinel: must be overwritten
    try {
        p = setupIso(1, isoNumbers, atomCounts, masses, probs);
    } catch (...) {
        escaped = true;
    }

    REQUIRE(!escaped,
            "setupIso() let a C++ exception cross the extern \"C\" boundary");
    REQUIRE(p == nullptr,
            "setupIso() with invalid input should return nullptr");

    deleteIso(p);  // must tolerate nullptr

    // A valid call must still succeed and return a usable handle.
    const double good_probs[2] = {0.5, 0.5};
    void* q = setupIso(1, isoNumbers, atomCounts, masses, good_probs);
    REQUIRE(q != nullptr, "setupIso() with valid input should return a handle");
    deleteIso(q);
}


// ===========================================================================
// Bug 6 (regression guard): pod_vector with initial_size == 0 must accept
// push_back / resize.  Commit 7e99ddf fixed this; this test pins the
// contract so it cannot silently regress.
// ===========================================================================
static void test_pod_vector_zero_capacity_push_back()
{
    pod_vector<int> v(0);
    EXPECT_NO_THROW(v.push_back(7), "push_back on zero-capacity pod_vector");
    REQUIRE(v.size() == 1, "size mismatch after push to zero-capacity pod_vector");
    REQUIRE(v[0] == 7, "value mismatch after push to zero-capacity pod_vector");

    pod_vector<int> v2(0);
    EXPECT_NO_THROW(v2.resize(16), "resize on zero-capacity pod_vector");
    REQUIRE(v2.size() == 16, "size mismatch after resize from zero capacity");
}


// ===========================================================================
// Bug 7 (documentation): the two-step init lists in Iso(int,...)
// (isoSpec++.cpp:61) and Marginal::Marginal(...) (marginalTrek++.cpp:165)
// leak the first-allocated member if the second allocation throws.
//
// This is only catchable by injecting an allocation failure at the right
// point in the init list — i.e. a custom operator new keyed on allocation
// size, or std::set_new_handler reaching a counter threshold.  We leave a
// stub here to keep the bug visible, but do not run a fake "test" that
// silently passes.
// ===========================================================================
static void test_iso_and_marginal_ctor_leak_documentation()
{
    std::fprintf(stderr,
        "DOC [%s]: ctor-leak bugs in isoSpec++.cpp:61 (Iso(int,...)) and "
        "marginalTrek++.cpp:165 (Marginal) leak the first-allocated array "
        "if the second allocation throws.  To test: install a custom "
        "operator new (or std::set_new_handler) that fails the Nth "
        "allocation and run under -fsanitize=leak.  Suggested harness:\n"
        "    static int gate = 0;\n"
        "    void* operator new[](std::size_t n) {\n"
        "        if (--gate == 0) throw std::bad_alloc();\n"
        "        return ::operator new(n);\n"
        "    }\n"
        "Then loop the gate value over [1..K] for each ctor.\n",
        __func__);
}


// ===========================================================================
// Main: run every test.
// ===========================================================================
int main()
{
    test_unsafe_pod_vector_swap_signature();
    test_binned_negative_offset_does_not_wrap();
    test_fasta_header_is_callable();
    test_empty_iso_generator_rejected_cleanly();
    test_cwrapper_rejects_invalid_input_with_nullptr();
    test_pod_vector_zero_capacity_push_back();
    test_iso_and_marginal_ctor_leak_documentation();

    if (g_failures == 0) {
        std::printf("All regression tests passed.\n");
        return 0;
    }
    std::printf("%d regression test(s) failed.\n", g_failures);
    return 1;
}

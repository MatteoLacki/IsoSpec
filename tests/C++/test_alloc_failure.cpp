// Allocation failures inside constructors: every Iso constructor, and
// addElement, is run with its Nth allocation made to throw std::bad_alloc, for
// N = 0, 1, 2, ... until one completes, and nothing allocated before the throw
// may outlive it. (A constructor that throws runs no destructor, so whatever it
// had already allocated into raw members used to leak.)
//
// This file replaces the global operator new/delete for the whole test binary.
// Disarmed, the replacement is malloc/free behind one thread_local check.
// Armed, it fails the Nth allocation and remembers every allocation made
// meanwhile until it is freed; aligned new and delete keep their defaults, and
// the library allocates nothing through them. Not under ThreadSanitizer, whose
// runtime defines operator new/delete itself (and which this is no test of).

#include <cstddef>
#include <cstdlib>
#include <new>
#include <string>
#include <vector>

#include "doctest.h"
#include "isoSpec++.h"

#if defined(__SANITIZE_THREAD__)
#define ISOSPEC_TEST_UNDER_TSAN 1
#elif defined(__has_feature)
#if __has_feature(thread_sanitizer)
#define ISOSPEC_TEST_UNDER_TSAN 1
#endif
#endif

#ifndef ISOSPEC_TEST_UNDER_TSAN

namespace {

constexpr std::size_t kMaxTracked = 4096;

struct AllocTracker {
    bool armed = false;
    long countdown = -1;       // allocations left before the failing one
    std::size_t live = 0;      // tracked allocations not yet freed
    void* tracked[kMaxTracked];
    bool overflowed = false;   // more live allocations than the table holds
};

thread_local AllocTracker t_alloc;

void* tracked_alloc(std::size_t size) {
    AllocTracker& t = t_alloc;
    if (t.armed) {
        if (t.countdown == 0)
            throw std::bad_alloc();
        --t.countdown;
    }
    void* p = std::malloc(size ? size : 1);
    if (p == nullptr)
        throw std::bad_alloc();
    if (t.armed) {
        if (t.live < kMaxTracked)
            t.tracked[t.live++] = p;
        else
            t.overflowed = true;
    }
    return p;
}

void tracked_free(void* p) noexcept {
    if (p == nullptr)
        return;
    AllocTracker& t = t_alloc;
    for (std::size_t ii = 0; ii < t.live; ++ii)
        if (t.tracked[ii] == p) {
            t.tracked[ii] = t.tracked[--t.live];
            break;
        }
    std::free(p);
}

// Fail the nth allocation from now (0 = the next one), and track until disarm().
void arm(long n) {
    t_alloc.armed = true;
    t_alloc.countdown = n;
    t_alloc.live = 0;
    t_alloc.overflowed = false;
}

void disarm() {
    t_alloc.armed = false;
    t_alloc.countdown = -1;
}

// Runs make() with its nth allocation failing, for n = 0, 1, ... until it gets
// through, checking each time that the allocations made meanwhile are all freed
// once the exception has propagated -- or, when make() succeeds, once what it
// built (destroyed inside make()) is gone. Returns how many allocations a
// successful run made.
template <typename Make>
long check_no_leak_on_any_failure(Make make) {
    for (long n = 0; n < 100000; ++n) {
        bool threw = false;
        arm(n);
        try {
            make();
        } catch (const std::bad_alloc&) {
            threw = true;
        }
        const long left = t_alloc.countdown;
        disarm();
        INFO("failing allocation #" << n);
        REQUIRE_FALSE(t_alloc.overflowed);
        CHECK(t_alloc.live == 0);
        if (!threw)
            return n - left;
    }
    FAIL("never got through");
    return -1;
}

}  // namespace

void* operator new(std::size_t size) { return tracked_alloc(size); }
void* operator new[](std::size_t size) { return tracked_alloc(size); }
void operator delete(void* p) noexcept { tracked_free(p); }
void operator delete[](void* p) noexcept { tracked_free(p); }
void operator delete(void* p, std::size_t) noexcept { tracked_free(p); }
void operator delete[](void* p, std::size_t) noexcept { tracked_free(p); }

using namespace IsoSpec;

TEST_CASE("Iso constructors free everything when any of their allocations fails") {
    const int isotope_numbers[] = {2, 3};
    const int atom_counts[] = {5, 7};
    const double c_masses[] = {12.0, 13.0033548352};
    const double c_probs[] = {0.9893, 0.0107};
    const double o_masses[] = {15.9949146196, 16.9991317566, 17.9991596128};
    const double o_probs[] = {0.99757, 0.00038, 0.00205};
    const double* masses[] = {c_masses, o_masses};
    const double* probs[] = {c_probs, o_probs};
    const double flat_masses[] = {12.0, 13.0033548352, 15.9949146196, 16.9991317566, 17.9991596128};
    const double flat_probs[] = {0.9893, 0.0107, 0.99757, 0.00038, 0.00205};

    SUBCASE("raw arrays, one per element") {
        CHECK(check_no_leak_on_any_failure([&] { Iso iso(2, isotope_numbers, atom_counts, masses, probs); }) > 2);
    }
    SUBCASE("raw arrays, flat") {
        CHECK(check_no_leak_on_any_failure([&] { Iso iso(2, isotope_numbers, atom_counts, flat_masses, flat_probs); }) > 2);
    }
    SUBCASE("formula") {
        CHECK(check_no_leak_on_any_failure([] { Iso iso("C100H202O7S2"); }) > 2);
    }
    SUBCASE("peptide sequence") {
        CHECK(check_no_leak_on_any_failure([] { Iso iso = Iso::FromFASTA("PEPTIDECM"); }) > 2);
    }
    SUBCASE("full copy") {
        const Iso original("C100H202O7S2");
        CHECK(check_no_leak_on_any_failure([&] { Iso copy(original, true); }) > 2);
    }
    SUBCASE("addElement") {
        // The Iso is built before arming; what addElement allocates must be
        // freed by the time the Iso is, whether the call threw or not.
        CHECK(check_no_leak_on_any_failure([&] {
            const long countdown = t_alloc.countdown;
            disarm();
            Iso* iso = new Iso("C10H20");
            t_alloc.armed = true;
            t_alloc.countdown = countdown;
            try {
                iso->addElement(3, 3, o_masses, o_probs);
            } catch (...) {
                delete iso;
                throw;
            }
            delete iso;
        }) > 0);
    }
}

#endif  // ISOSPEC_TEST_UNDER_TSAN

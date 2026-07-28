// aligned_unique_ptr: the small aligned_alloc-backed buffer, its transition
// to a VM-backed (mmap/mremap/mach_vm_remap/VirtualAlloc) region once it
// crosses a page-size threshold, and the two ownership-transfer paths
// (release(), which is always free()-compatible, and release_with_deleter(),
// which never copies but hands back a matching deleter instead).

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

#include "doctest.h"
#include "aligned_ptr.h"
#include "platform.h"  // DOUBLE_SIMD_ALIGNMENT

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC && ISOSPEC_WE_ARE_ON_UNIX_YAY && !(defined(__APPLE__) && defined(__MACH__))
#include <sys/mman.h>
#ifdef MAP_FIXED_NOREPLACE
#define ISOSPEC_TEST_HAVE_MAP_FIXED_NOREPLACE 1
#else
#define ISOSPEC_TEST_HAVE_MAP_FIXED_NOREPLACE 0
#endif
#else
#define ISOSPEC_TEST_HAVE_MAP_FIXED_NOREPLACE 0
#endif

namespace {

bool is_aligned(const void* p, std::size_t alignment) {
    return (reinterpret_cast<std::uintptr_t>(p) % alignment) == 0;
}

// Same pattern as test_generators.cpp's Tally: doctest's assertion machinery
// costs far more than a single comparison, so a loop over many elements
// tallies pass/fail and asserts once instead of per element.
struct Tally {
    std::size_t checked = 0;
    std::size_t failed = 0;
    void expect(bool ok) { ++checked; if (!ok) ++failed; }
    bool clean() const { return failed == 0 && checked > 0; }
};

// The runtime page-size threshold aligned_unique_ptr::realloc() switches
// backend at; only meaningful (and only compiled in) where the VM backend
// exists at all.
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
std::size_t page_size() { return aligned_ptr_detail::os_page_size(); }
#endif

}  // namespace

TEST_CASE("aligned_unique_ptr: construction, indexing and alignment") {
    SUBCASE("alignment 16") {
        aligned_unique_ptr<double, 16> p(10);
        CHECK(static_cast<bool>(p));
        CHECK(is_aligned(p.get(), 16));
        for (int i = 0; i < 10; i++) p[i] = i * 1.5;
        for (int i = 0; i < 10; i++) CHECK(p[i] == doctest::Approx(i * 1.5));
    }
    SUBCASE("alignment 64") {
        aligned_unique_ptr<double, 64> p(3);
        CHECK(is_aligned(p.get(), 64));
    }
    SUBCASE("the actual production instantiation") {
        aligned_unique_ptr<double, DOUBLE_SIMD_ALIGNMENT> p(5);
        CHECK(is_aligned(p.get(), DOUBLE_SIMD_ALIGNMENT));
        p[4] = 42.0;
        CHECK(p[4] == 42.0);
    }
    SUBCASE("zero-length construction is empty") {
        aligned_unique_ptr<double, 32> p(0);
        CHECK_FALSE(static_cast<bool>(p));
    }
}

TEST_CASE("aligned_unique_ptr: reset does not preserve contents") {
    aligned_unique_ptr<double, 32> p(4);
    for (int i = 0; i < 4; i++) p[i] = 99.0;

    p.reset(8);
    CHECK(static_cast<bool>(p));
    // Contents are unspecified after reset(), but the buffer must be usable.
    for (int i = 0; i < 8; i++) p[i] = i + 1.0;
    for (int i = 0; i < 8; i++) CHECK(p[i] == i + 1.0);

    p.reset(0);
    CHECK_FALSE(static_cast<bool>(p));
}

TEST_CASE("aligned_unique_ptr: move construction and move assignment") {
    aligned_unique_ptr<double, 32> a(4);
    for (int i = 0; i < 4; i++) a[i] = i + 10.0;

    aligned_unique_ptr<double, 32> b(std::move(a));
    CHECK_FALSE(static_cast<bool>(a));  // NOLINT(bugprone-use-after-move) -- deliberately checking moved-from state
    CHECK(static_cast<bool>(b));
    for (int i = 0; i < 4; i++) CHECK(b[i] == i + 10.0);

    aligned_unique_ptr<double, 32> c(2);
    c[0] = -1.0;
    c = std::move(b);
    CHECK_FALSE(static_cast<bool>(b));
    for (int i = 0; i < 4; i++) CHECK(c[i] == i + 10.0);

    // Self-move-assignment must not destroy the object (the guard is
    // `this != &other`, easy to regress by refactoring). Routed through a
    // pointer so the compiler can't flag it as a literal self-move.
    aligned_unique_ptr<double, 32>* self = &c;
    c = std::move(*self);
    for (int i = 0; i < 4; i++) CHECK(c[i] == i + 10.0);
}

TEST_CASE("aligned_unique_ptr: realloc grows and shrinks within the small backend") {
    aligned_unique_ptr<double, 32> p(4);
    for (int i = 0; i < 4; i++) p[i] = i + 1.0;

    p.realloc(8);
    for (int i = 0; i < 4; i++) CHECK(p[i] == i + 1.0);
    for (int i = 4; i < 8; i++) p[i] = i + 1.0;

    p.realloc(6);  // shrink, still small
    for (int i = 0; i < 6; i++) CHECK(p[i] == i + 1.0);

    p.realloc(2);
    CHECK(p[0] == 1.0);
    CHECK(p[1] == 2.0);
}

TEST_CASE("aligned_unique_ptr: realloc(0) empties the buffer and it is reusable") {
    aligned_unique_ptr<double, 32> p(4);
    p[0] = 7.0;
    p.realloc(0);
    CHECK_FALSE(static_cast<bool>(p));

    p.realloc(3);
    CHECK(static_cast<bool>(p));
    p[0] = 8.0;
    CHECK(p[0] == 8.0);
}

TEST_CASE("aligned_unique_ptr: n*sizeof(T) overflow is rejected, not wrapped") {
    // A request this large must throw std::bad_alloc rather than silently
    // wrapping size_t and handing back a too-small buffer for the caller to
    // overrun (classic integer-overflow-to-heap-overflow pattern).
    using D32 = aligned_unique_ptr<double, 32>;
    constexpr std::size_t huge = (std::numeric_limits<std::size_t>::max)() / sizeof(double) + 1;
    CHECK_THROWS_AS(D32{huge}, std::bad_alloc);

    D32 p(4);
    CHECK_THROWS_AS(p.realloc(huge), std::bad_alloc);
    // Must still be usable / unharmed after the rejected realloc.
    CHECK(static_cast<bool>(p));
}

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC

TEST_CASE("aligned_unique_ptr: realloc crossing the VM threshold preserves data") {
    const std::size_t page_elems = page_size() / sizeof(double);

    aligned_unique_ptr<double, 32> p(8);
    for (int i = 0; i < 8; i++) p[i] = i + 1.0;

    // Still below threshold: one element short of a full page.
    p.realloc(page_elems - 1);
    for (int i = 0; i < 8; i++) CHECK(p[i] == i + 1.0);
    for (std::size_t i = 8; i < page_elems - 1; i++) p[i] = static_cast<double>(i);

    // Crosses the threshold now.
    const std::size_t past = page_elems + 16;
    p.realloc(past);
    for (int i = 0; i < 8; i++) CHECK(p[i] == i + 1.0);
    for (std::size_t i = 8; i < page_elems - 1; i++) CHECK(p[i] == static_cast<double>(i));
    for (std::size_t i = page_elems - 1; i < past; i++) p[i] = 3.5 * static_cast<double>(i);

    // Grow again while already mapped -- exercises mremap/mach_vm_remap/
    // VirtualAlloc-commit rather than the crossing-transition copy above.
    const std::size_t bigger = past * 4;
    p.realloc(bigger);
    for (int i = 0; i < 8; i++) CHECK(p[i] == i + 1.0);
    for (std::size_t i = 8; i < page_elems - 1; i++) CHECK(p[i] == static_cast<double>(i));
    for (std::size_t i = page_elems - 1; i < past; i++) CHECK(p[i] == doctest::Approx(3.5 * static_cast<double>(i)));
    for (std::size_t i = past; i < bigger; i++) p[i] = static_cast<double>(i) - 1000.0;

    // Shrink while mapped.
    const std::size_t smaller = past + 4;
    p.realloc(smaller);
    for (std::size_t i = past; i < smaller; i++) CHECK(p[i] == doctest::Approx(static_cast<double>(i) - 1000.0));

    // Back down to zero and re-grow from empty.
    p.realloc(0);
    CHECK_FALSE(static_cast<bool>(p));
    p.realloc(4);
    p[0] = 123.0;
    CHECK(p[0] == 123.0);
}

TEST_CASE("aligned_unique_ptr: repeated grow/shrink cycles around the threshold keep data intact (stress)") {
    const std::size_t page_elems = page_size() / sizeof(double);
    std::mt19937 rng(12345);

    aligned_unique_ptr<double, 32> p(1);
    p[0] = 0.0;
    std::vector<double> shadow(1, 0.0);

    Tally t;
    std::uniform_int_distribution<std::size_t> size_dist(1, page_elems * 3);
    for (int round = 0; round < 40; round++) {
        std::size_t new_n = size_dist(rng);
        p.realloc(new_n);

        std::size_t common = std::min(shadow.size(), new_n);
        shadow.resize(new_n, 0.0);  // logical shadow only, for bookkeeping below

        for (std::size_t i = 0; i < common; i++)
            t.expect(p[i] == doctest::Approx(shadow[i]));

        for (std::size_t i = common; i < new_n; i++) {
            double v = static_cast<double>(round) + static_cast<double>(i) * 0.001;
            p[i] = v;
            shadow[i] = v;
        }
    }
    CHECK(t.clean());
}

#endif  // ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC

TEST_CASE("aligned_unique_ptr: release() is always free()-compatible") {
    SUBCASE("small backend: no copy needed, still free()-able") {
        aligned_unique_ptr<double, 32> p(4);
        for (int i = 0; i < 4; i++) p[i] = i + 1.0;
        double* raw = p.release();
        CHECK_FALSE(static_cast<bool>(p));
        REQUIRE(raw != nullptr);
        for (int i = 0; i < 4; i++) CHECK(raw[i] == i + 1.0);
        std::free(raw);  // must not crash / must not be flagged by ASan
    }

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
    SUBCASE("VM-backed: release() materialises a free()-compatible copy") {
        const std::size_t n = page_size() / sizeof(double) + 64;
        aligned_unique_ptr<double, 32> p(1);
        p.realloc(n);  // crosses the threshold -> VM-backed
        for (std::size_t i = 0; i < n; i++) p[i] = static_cast<double>(i);

        double* raw = p.release();
        CHECK_FALSE(static_cast<bool>(p));
        REQUIRE(raw != nullptr);
        for (std::size_t i = 0; i < n; i++) CHECK(raw[i] == static_cast<double>(i));
        std::free(raw);  // would crash/corrupt under ASan if this were still an mmap'd pointer

        // The instance must remain usable after release().
        p.realloc(2);
        p[0] = 5.0;
        CHECK(p[0] == 5.0);
    }
#endif

    SUBCASE("empty instance") {
        aligned_unique_ptr<double, 32> p;
        CHECK(p.release() == nullptr);
    }
}

TEST_CASE("aligned_unique_ptr: release_with_deleter() avoids copying and still frees correctly") {
    SUBCASE("small backend") {
        aligned_unique_ptr<double, 32> p(4);
        for (int i = 0; i < 4; i++) p[i] = i + 1.0;

        auto rr = p.release_with_deleter();
        CHECK_FALSE(static_cast<bool>(p));
        REQUIRE(rr.ptr != nullptr);
        for (int i = 0; i < 4; i++) CHECK(rr.ptr[i] == i + 1.0);
        rr.deleter(rr.ptr, rr.size);  // must not crash under ASan
    }

    SUBCASE("empty instance: deleter is a safe no-op") {
        aligned_unique_ptr<double, 32> p;
        auto rr = p.release_with_deleter();
        CHECK(rr.ptr == nullptr);
        rr.deleter(rr.ptr, rr.size);  // std::free(nullptr) is well-defined
    }

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
    SUBCASE("VM-backed: no copy, and the deleter frees the *actual* backing size") {
        const std::size_t n = page_size() / sizeof(double) + 64;
        aligned_unique_ptr<double, 32> p(1);
        p.realloc(n);
        for (std::size_t i = 0; i < n; i++) p[i] = static_cast<double>(i) * 2.0;

        void* base = static_cast<void*>(p.get());
        auto rr = p.release_with_deleter();
        CHECK(rr.ptr == base);
        for (std::size_t i = 0; i < n; i++) CHECK(rr.ptr[i] == static_cast<double>(i) * 2.0);

        rr.deleter(rr.ptr, rr.size);

#if ISOSPEC_TEST_HAVE_MAP_FIXED_NOREPLACE
        // The deleter must munmap the *whole* underlying region (rr.size is
        // the actual, page-rounded mapping size, not just the logical
        // request) -- if it under-reported the size, part of the old
        // mapping would still occupy this address range and claiming it
        // again with MAP_FIXED_NOREPLACE would fail with EEXIST.
        void* reclaimed = ::mmap(base, rr.size, PROT_READ | PROT_WRITE,
                                  MAP_PRIVATE | MAP_ANONYMOUS | MAP_FIXED_NOREPLACE, -1, 0);
        REQUIRE(reclaimed != MAP_FAILED);
        CHECK(reclaimed == base);
        ::munmap(reclaimed, rr.size);
#endif
    }
#endif
}

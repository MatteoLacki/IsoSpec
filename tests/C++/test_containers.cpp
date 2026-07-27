// The internal containers and allocators: pod_vector, unsafe_pod_vector,
// Allocator<T>, DirtyAllocator, and the summators.
//
// These are small, hot, hand-rolled data structures with no bounds checking,
// so their contracts are worth pinning down explicitly.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <numeric>
#include <utility>
#include <vector>

#include "doctest.h"
#include "allocator.h"
#include "dirtyAllocator.h"
#include "misc.h"
#include "pod_vector.h"
#include "summator.h"

using namespace IsoSpec;

TEST_CASE("pod_vector: growth, indexing and iteration") {
    pod_vector<int> v;
    CHECK(v.size() == 0);
    CHECK(v.empty());
    CHECK(v.capacity() == 16);  // documented default

    for (int i = 0; i < 1000; ++i) v.push_back(i);

    CHECK(v.size() == 1000);
    CHECK_FALSE(v.empty());
    CHECK(v.capacity() >= 1000);
    CHECK(v.front() == 0);
    CHECK(v.back() == 999);
    CHECK(v[500] == 500);
    CHECK(v.data()[999] == 999);

    // Iterators cover exactly the stored range.
    CHECK(v.end() - v.begin() == 1000);
    CHECK(std::accumulate(v.cbegin(), v.cend(), 0) == 999 * 1000 / 2);

    v.pop_back();
    CHECK(v.size() == 999);
    CHECK(v.back() == 998);

    v.clear();
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 0);
    // A cleared vector is still usable.
    v.push_back(42);
    CHECK(v.size() == 1);
    CHECK(v[0] == 42);
}

TEST_CASE("pod_vector: reserve and resize") {
    pod_vector<double> v(4);
    CHECK(v.capacity() == 4);

    v.reserve(100);
    CHECK(v.capacity() >= 100);
    CHECK(v.size() == 0);
    // reserve() below the current capacity is a no-op, not a shrink.
    const std::size_t cap = v.capacity();
    v.reserve(2);
    CHECK(v.capacity() == cap);

    v.resize(50);
    CHECK(v.size() == 50);
    for (std::size_t i = 0; i < 50; ++i) v[i] = static_cast<double>(i);
    CHECK(v[49] == 49.0);

    // Growing past the capacity reallocates but preserves the contents.
    v.resize(5000);
    CHECK(v.size() == 5000);
    CHECK(v.capacity() >= 5000);
    CHECK(v[49] == 49.0);

    // resize_and_wipe zeroes the new tail only.
    pod_vector<int> w(2);
    w.push_back(7);
    w.resize_and_wipe(100);
    CHECK(w.size() == 100);
    CHECK(w[0] == 7);
    for (std::size_t i = 1; i < 100; ++i) CHECK(w[i] == 0);
}

TEST_CASE("pod_vector: move semantics and swap") {
    pod_vector<int> a(8);
    for (int i = 0; i < 8; ++i) a.push_back(i);
    const int* buf = a.data();

    pod_vector<int> b(std::move(a));
    CHECK(b.size() == 8);
    CHECK(b.data() == buf);   // buffer stolen, not copied
    CHECK(a.size() == 0);     // moved-from vector is empty and safe to destroy

    pod_vector<int> c(4);
    c.push_back(99);
    c = std::move(b);
    CHECK(c.size() == 8);
    CHECK(c.data() == buf);

    pod_vector<int> d(2);
    d.push_back(1);
    d.push_back(2);
    c.swap(d);
    CHECK(c.size() == 2);
    CHECK(d.size() == 8);
    CHECK(d[7] == 7);
}

TEST_CASE("pod_vector: nocheck_push_back within reserved capacity") {
    pod_vector<int> v(4);
    v.reserve(64);
    for (int i = 0; i < 64; ++i) v.nocheck_push_back(i);
    CHECK(v.size() == 64);
    CHECK(v[63] == 63);
}

TEST_CASE("unsafe_pod_vector: the same contract, minus the constructor") {
    unsafe_pod_vector<double> v;
    v.init(4);
    CHECK(v.size() == 0);
    CHECK(v.capacity() == 4);
    CHECK(v.empty());

    for (int i = 0; i < 100; ++i) v.push_back(i * 0.5);
    CHECK(v.size() == 100);
    CHECK(v.front() == 0.0);
    CHECK(v.back() == 49.5);
    CHECK(v[10] == 5.0);
    CHECK(std::is_sorted(v.begin(), v.end()));

    v.reserve(1000);
    CHECK(v.capacity() >= 1000);
    CHECK(v.size() == 100);
    CHECK(v[99] == 49.5);

    v.resize(200);
    CHECK(v.size() == 200);
    v.pop_back();
    CHECK(v.size() == 199);

    v.clear();
    CHECK(v.size() == 0);

    // init() with no argument zeroes the whole object: an empty, unallocated
    // vector that push_back can still grow.
    unsafe_pod_vector<int> w;
    w.init();
    CHECK(w.size() == 0);
    CHECK(w.capacity() == 0);
    w.push_back(5);
    CHECK(w.size() == 1);
    CHECK(w[0] == 5);
    w.clear();
}

TEST_CASE("pod_vector can be built from an unsafe_pod_vector") {
    unsafe_pod_vector<int> u;
    u.init(4);
    for (int i = 0; i < 10; ++i) u.push_back(i);
    const int* buf = u.data();

    pod_vector<int> v(std::move(u));
    CHECK(v.size() == 10);
    CHECK(v.data() == buf);
    CHECK(v[9] == 9);
    CHECK(u.size() == 0);
}

TEST_CASE("Allocator hands out distinct, correctly-sized cells") {
    const int dim = 5;
    const int tabSize = 8;  // small, to force several table shifts
    Allocator<int> alloc(dim, tabSize);

    std::vector<int*> cells;
    for (int i = 0; i < 100; ++i) {
        int* c = alloc.newConf();
        REQUIRE(c != nullptr);
        for (int j = 0; j < dim; ++j) c[j] = i * dim + j;
        cells.push_back(c);
    }

    // Every cell must still hold what was written into it: no overlap.
    for (int i = 0; i < 100; ++i)
        for (int j = 0; j < dim; ++j)
            REQUIRE(cells[i][j] == i * dim + j);

    // makeCopy duplicates a configuration into fresh storage.
    int source[dim] = {1, 2, 3, 4, 5};
    int* copy = alloc.makeCopy(source);
    REQUIRE(copy != nullptr);
    CHECK(copy != source);
    CHECK(memcmp(copy, source, sizeof(source)) == 0);
    source[0] = 99;
    CHECK(copy[0] == 1);
}

TEST_CASE("DirtyAllocator hands out distinct cells") {
    const int dim = 4;
    DirtyAllocator alloc(dim, 8);

    std::vector<void*> cells;
    for (int i = 0; i < 100; ++i) {
        void* c = alloc.newConf();
        REQUIRE(c != nullptr);
        // The cell layout is [double lprob][int conf[dim]] — see misc.h.
        *reinterpret_cast<double*>(c) = i;
        int* conf = getConf(c);
        for (int j = 0; j < dim; ++j) conf[j] = i * dim + j;
        cells.push_back(c);
    }

    for (int i = 0; i < 100; ++i) {
        REQUIRE(getLProb(cells[i]) == static_cast<double>(i));
        const int* conf = getConf(cells[i]);
        for (int j = 0; j < dim; ++j) REQUIRE(conf[j] == i * dim + j);
    }
}

TEST_CASE("summators: Kahan and Shewchuk beat naive summation") {
    // A classic ill-conditioned sum: 1 + N*eps, where each eps is lost
    // individually in naive double addition.
    const double eps = 1e-17;
    const int n = 1000000;

    TSummator naive;
    Summator kahan;
    SSummator shewchuk;

    naive.add(1.0);
    kahan.add(1.0);
    shewchuk.add(1.0);
    for (int i = 0; i < n; ++i) {
        naive.add(eps);
        kahan.add(eps);
        shewchuk.add(eps);
    }

    const double expected = 1.0 + n * eps;
    CHECK(naive.get() == 1.0);  // every increment vanished
    CHECK(kahan.get() == doctest::Approx(expected).epsilon(1e-15));
    CHECK(shewchuk.get() == doctest::Approx(expected).epsilon(1e-15));
}

TEST_CASE("summators agree on well-conditioned input") {
    TSummator naive;
    Summator kahan;
    SSummator shewchuk;
    double reference = 0.0;
    for (int i = 1; i <= 1000; ++i) {
        const double x = 1.0 / i;
        naive.add(x);
        kahan.add(x);
        shewchuk.add(x);
        reference += x;
    }
    CHECK(naive.get() == doctest::Approx(reference));
    CHECK(kahan.get() == doctest::Approx(reference));
    CHECK(shewchuk.get() == doctest::Approx(reference));

    // Empty sums are zero.
    CHECK(TSummator().get() == 0.0);
    CHECK(Summator().get() == 0.0);
    CHECK(SSummator().get() == 0.0);
}

TEST_CASE("SSummator handles cancellation exactly") {
    // Shewchuk's algorithm is exact for sums that cancel completely.
    SSummator s;
    s.add(1e100);
    s.add(1.0);
    s.add(-1e100);
    CHECK(s.get() == doctest::Approx(1.0));

    // The copy constructor carries the partial sums along.
    SSummator copy(s);
    CHECK(copy.get() == s.get());
}

TEST_CASE("misc: get_order and impose_order sort in tandem") {
    double keys[5] = {3.0, 1.0, 4.0, 1.5, 2.0};
    int payload[5] = {30, 10, 40, 15, 20};
    double keys_copy[5];
    memcpy(keys_copy, keys, sizeof(keys));

    size_t* order = get_order(keys_copy, 5);
    // order[i] is the index of the i-th smallest key.
    CHECK(keys[order[0]] == 1.0);
    CHECK(keys[order[4]] == 4.0);
    for (int i = 1; i < 5; ++i) CHECK(keys[order[i - 1]] <= keys[order[i]]);

    impose_order(order, 5, keys_copy, payload);
    CHECK(std::is_sorted(keys_copy, keys_copy + 5));
    // The payload followed its key.
    for (int i = 0; i < 5; ++i) CHECK(payload[i] == static_cast<int>(keys_copy[i] * 10));
    delete[] order;

    // get_inverse_order sorts descending.
    double keys2[4] = {1.0, 4.0, 2.0, 3.0};
    size_t* desc = get_inverse_order(keys2, 4);
    for (int i = 1; i < 4; ++i) CHECK(keys2[desc[i - 1]] >= keys2[desc[i]]);
    delete[] desc;
}

TEST_CASE("misc: array_copy variants") {
    const int source[4] = {1, 2, 3, 4};

    int* copy = array_copy<int>(source, 4);
    CHECK(memcmp(copy, source, sizeof(source)) == 0);
    delete[] copy;

    int* mcopy = array_copy_malloc<int>(source, 4);
    CHECK(memcmp(mcopy, source, sizeof(source)) == 0);
    free(mcopy);

    // The _nptr variants pass a null pointer straight through.
    CHECK(array_copy_nptr<int>(nullptr, 4) == nullptr);
    CHECK(array_copy_nptr_malloc<int>(nullptr, 4) == nullptr);

    int* nptr_copy = array_copy_nptr<int>(source, 4);
    CHECK(nptr_copy != nullptr);
    CHECK(nptr_copy[3] == 4);
    delete[] nptr_copy;
}

TEST_CASE("misc: realloc_append grows an array by one") {
    int* array = new int[3]{1, 2, 3};
    realloc_append<int>(&array, 4, 3);
    CHECK(array[0] == 1);
    CHECK(array[3] == 4);
    realloc_append<int>(&array, 5, 4);
    CHECK(array[4] == 5);
    delete[] array;
}

TEST_CASE("misc: calc_mass and unnormalized_logProb") {
    const int conf[3] = {2, 1, 0};
    const double masses[3] = {10.0, 20.0, 30.0};
    CHECK(calc_mass(conf, masses, 3) == doctest::Approx(40.0));

    // unnormalized_logProb = sum(-log(k_i!) + k_i * log p_i), i.e. the
    // multinomial term without the n! factor.
    const double lprobs[3] = {std::log(0.5), std::log(0.3), std::log(0.2)};
    const double expected = -std::lgamma(3.0) + 2 * lprobs[0] + 1 * lprobs[1];
    CHECK(unnormalized_logProb(conf, lprobs, 3) == doctest::Approx(expected));
}

TEST_CASE("misc: combinedSum indexes per-dimension value tables") {
    std::vector<double> a = {1.0, 2.0, 3.0};
    std::vector<double> b = {10.0, 20.0};
    const std::vector<double>* tables[2] = {&a, &b};
    const int conf[2] = {2, 1};
    CHECK(combinedSum(conf, tables, 2) == doctest::Approx(23.0));

    pod_vector<double> pa(4), pb(4);
    for (double x : a) pa.push_back(x);
    for (double x : b) pb.push_back(x);
    const pod_vector<double>* ptables[2] = {&pa, &pb};
    CHECK(combinedSum(conf, ptables, 2) == doctest::Approx(23.0));
}

TEST_CASE("misc: quickselect finds the n-th smallest by log-probability") {
    // Cells are [double lprob][int conf...]; quickselect orders by the double.
    const int dim = 1;
    DirtyAllocator alloc(dim, 64);
    const double lprobs[9] = {5.0, 1.0, 9.0, 3.0, 7.0, 2.0, 8.0, 4.0, 6.0};
    void* cells[9];
    for (int i = 0; i < 9; ++i) {
        cells[i] = alloc.newConf();
        *reinterpret_cast<double*>(cells[i]) = lprobs[i];
    }

    for (std::size_t n = 0; n < 9; ++n) {
        void* cells_copy[9];
        memcpy(cells_copy, cells, sizeof(cells));
        void* got = quickselect(cells_copy, n, 0, 9);
        CHECK(getLProb(got) == doctest::Approx(static_cast<double>(n + 1)));
        // Everything before position n is no larger than the selected value.
        for (std::size_t i = 0; i < n; ++i)
            CHECK(getLProb(cells_copy[i]) <= getLProb(got));
    }

    // Degenerate range: start == end returns that element.
    void* single[1] = {cells[0]};
    CHECK(quickselect(single, 0, 0, 0) == cells[0]);
}

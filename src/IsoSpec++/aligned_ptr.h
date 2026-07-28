#pragma once

#include <cstddef>      // std::size_t, std::nullptr_t
#include <cstdint>      // std::uintptr_t (only if doing runtime alignment checks)
#include <cstdlib>      // std::aligned_alloc, std::malloc, std::free
#include <cstring>      // std::memcpy
#include <algorithm>    // std::min
#include <limits>       // std::numeric_limits
#include <memory>       // std::assume_aligned (C++20)
#include <new>          // std::bad_alloc
#include <utility>      // std::exchange
#include <type_traits>  // std::is_trivial_v

#include "platform.h"   // ISOSPEC_WE_ARE_ON_WINDOWS / ISOSPEC_WE_ARE_ON_UNIX_YAY

// aligned_unique_ptr's realloc() grows/shrinks without copying once the
// allocation is big enough to be worth mapping directly from the OS. Below
// that size it stays on the plain aligned_alloc/copy/free path. Deliberately
// bypasses this project's own mman.h shim/ISOSPEC_GOT_MMAN machinery (used
// elsewhere for plain huge-page mmap): that layer has no remap concept at
// all on Windows and no equivalent of mremap on Darwin, so each platform's
// real growth primitive is used directly instead:
//   - Linux / other real-mman UNIX: mmap + mremap(MREMAP_MAYMOVE) + munmap.
//   - Apple: mmap/munmap for create/destroy (Darwin has no mremap); growth
//     re-wires pages into a fresh mapping via mach_vm_remap instead of
//     copying, shrink just munmaps the tail in place.
//   - Windows: VirtualAlloc(MEM_RESERVE) a generously oversized address
//     range once, then VirtualAlloc(MEM_COMMIT)/VirtualFree(MEM_DECOMMIT)
//     to grow/shrink within it -- the base address never moves, so there is
//     no copy at all in the common case (only if growth ever outgrows the
//     initial reservation does it fall back to allocating a fresh, larger
//     reservation and copying, same as a real mremap() that has to move).
// On any other platform (no mmap-family primitive available at all) the
// VM-backed path is simply never engaged; realloc() keeps using the
// aligned_alloc/copy/free path regardless of size.
#if ISOSPEC_WE_ARE_ON_WINDOWS || ISOSPEC_WE_ARE_ON_UNIX_YAY
#define ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC 1
#else
#define ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC 0
#endif

// The small-allocation backend needs an aligned allocator that can also be
// deallocated with a plain std::free() (release() promises exactly that).
// glibc/libc++'s std::aligned_alloc satisfies both at once. MSVC's CRT has
// no C11 aligned_alloc at all, only _aligned_malloc/_aligned_free -- and
// memory from those two is NOT free()-compatible (the CRT stores its own
// offset bookkeeping ahead of the returned pointer, which plain free() does
// not know how to unwind). So on MSVC the small backend uses
// _aligned_malloc/_aligned_free for its own storage (still genuinely
// Alignment-aligned, as marginalTrek++.h's use of this class relies on),
// and release() falls back to materialising a fresh plain-malloc()'d copy
// to keep its free()-compatibility promise -- the same thing it already
// does to escape a VM-backed region.
#if defined(_MSC_VER)
#define ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE 0
#else
#define ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE 1
#endif

#if ISOSPEC_WE_ARE_ON_WINDOWS
    #include <windows.h>
    #if defined(_MSC_VER)
        #include <malloc.h>   // _aligned_malloc / _aligned_free
    #endif
#elif ISOSPEC_WE_ARE_ON_UNIX_YAY
    #include <sys/mman.h>   // real system mmap/munmap (+ mremap on Linux)
    #include <unistd.h>     // sysconf
    #if defined(__APPLE__) && defined(__MACH__)
        #include <mach/mach.h>
        #include <mach/mach_vm.h>
    #endif
#endif


#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
namespace aligned_ptr_detail {

inline std::size_t os_page_size() noexcept
{
#if ISOSPEC_WE_ARE_ON_WINDOWS
    static const std::size_t sz = [] {
        SYSTEM_INFO si;
        GetSystemInfo(&si);
        return static_cast<std::size_t>(si.dwPageSize);
    }();
#else
    static const std::size_t sz = static_cast<std::size_t>(sysconf(_SC_PAGESIZE));
#endif
    return sz;
}

inline std::size_t round_up_to_page(std::size_t bytes) noexcept
{
    const std::size_t ps = os_page_size();
    return (bytes + ps - 1) / ps * ps;
}

// Releases a VM-backed block given just (ptr, bytes), without the rest of a
// vm_region's bookkeeping -- used by aligned_unique_ptr::release_with_deleter()
// once ownership has been handed off to a caller who only has those two
// values. bytes must be the actual backing size (region_.size / .reserved
// above), not a logical/unrounded request size; Windows ignores it (a
// placeholder-free VirtualFree(MEM_RELEASE) always frees the whole
// reservation from its base).
inline void vm_release_raw(void* ptr, std::size_t bytes) noexcept
{
#if ISOSPEC_WE_ARE_ON_WINDOWS
    (void)bytes;
    ::VirtualFree(ptr, 0, MEM_RELEASE);
#else
    ::munmap(ptr, bytes);
#endif
}


#if ISOSPEC_WE_ARE_ON_WINDOWS
// Reservation only consumes address space, not RAM or page tables, so we can
// afford to be wasteful up front in order to (almost) never have to move
// the block: growth just commits further pages of the same reservation.
struct vm_region {
    void* base = nullptr;
    std::size_t reserved = 0;
    std::size_t committed = 0;
};

inline std::size_t reservation_size_for(std::size_t bytes) noexcept
{
    constexpr std::size_t min_reservation = std::size_t{1} << 26;  // 64 MiB
    constexpr std::size_t growth_factor = 16;
    std::size_t reserved = round_up_to_page(bytes) * growth_factor;
    if (reserved < min_reservation)
        reserved = min_reservation;
    return reserved;
}

inline vm_region vm_create(std::size_t bytes)
{
    vm_region r;
    r.committed = round_up_to_page(bytes);
    r.reserved = reservation_size_for(bytes);

    r.base = ::VirtualAlloc(nullptr, r.reserved, MEM_RESERVE, PAGE_NOACCESS);
    if (!r.base)
        throw std::bad_alloc();

    if (r.committed > 0 &&
        !::VirtualAlloc(r.base, r.committed, MEM_COMMIT, PAGE_READWRITE)) {
        ::VirtualFree(r.base, 0, MEM_RELEASE);
        throw std::bad_alloc();
    }
    return r;
}

inline void vm_destroy(vm_region& r) noexcept
{
    if (r.base)
        ::VirtualFree(r.base, 0, MEM_RELEASE);
    r = vm_region{};
}

inline void vm_resize(vm_region& r, std::size_t new_bytes)
{
    const std::size_t new_committed = round_up_to_page(new_bytes);
    if (new_committed == r.committed)
        return;

    if (new_committed > r.reserved) {
        // Outgrew the initial reservation (rare): there is no way to extend
        // a VirtualAlloc reservation in place, so relocate -- same as a real
        // mremap() that cannot grow in place and has to move.
        vm_region moved = vm_create(new_bytes);
        std::memcpy(moved.base, r.base, r.committed);
        vm_destroy(r);
        r = moved;
        return;
    }

    if (new_committed > r.committed) {
        if (!::VirtualAlloc(static_cast<char*>(r.base) + r.committed,
                             new_committed - r.committed, MEM_COMMIT, PAGE_READWRITE))
            throw std::bad_alloc();
    } else {
        ::VirtualFree(static_cast<char*>(r.base) + new_committed,
                       r.committed - new_committed, MEM_DECOMMIT);
    }
    r.committed = new_committed;
}

#else  // ISOSPEC_WE_ARE_ON_UNIX_YAY

struct vm_region {
    void* base = nullptr;
    std::size_t size = 0;
};

inline vm_region vm_create(std::size_t bytes)
{
    vm_region r;
    r.size = round_up_to_page(bytes);
    r.base = ::mmap(nullptr, r.size, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (r.base == MAP_FAILED)
        throw std::bad_alloc();
    return r;
}

inline void vm_destroy(vm_region& r) noexcept
{
    if (r.base)
        ::munmap(r.base, r.size);
    r = vm_region{};
}

#if defined(__APPLE__) && defined(__MACH__)

inline void vm_resize(vm_region& r, std::size_t new_bytes)
{
    const std::size_t new_size = round_up_to_page(new_bytes);
    if (new_size == r.size)
        return;

    if (new_size < r.size) {
        // Shrink in place: base address is untouched, just drop the tail.
        ::munmap(static_cast<char*>(r.base) + new_size, r.size - new_size);
        r.size = new_size;
        return;
    }

    // Grow: Darwin has no mremap(). Allocate a fresh region and use
    // mach_vm_remap to re-wire the old pages (not copy them) into its head,
    // then release the old mapping.
    mach_vm_address_t new_addr = 0;
    if (mach_vm_allocate(mach_task_self(), &new_addr, new_size, VM_FLAGS_ANYWHERE) != KERN_SUCCESS)
        throw std::bad_alloc();

    mach_vm_address_t target = new_addr;
    vm_prot_t cur_prot, max_prot;
    const kern_return_t kr = mach_vm_remap(
        mach_task_self(), &target, r.size, 0, VM_FLAGS_OVERWRITE,
        mach_task_self(), reinterpret_cast<mach_vm_address_t>(r.base),
        FALSE /* share the same physical pages, don't copy */,
        &cur_prot, &max_prot, VM_INHERIT_COPY);
    if (kr != KERN_SUCCESS) {
        mach_vm_deallocate(mach_task_self(), new_addr, new_size);
        throw std::bad_alloc();
    }

    ::munmap(r.base, r.size);
    r.base = reinterpret_cast<void*>(new_addr);
    r.size = new_size;
}

#else  // Linux / other real-mman UNIX with an actual mremap()

inline void vm_resize(vm_region& r, std::size_t new_bytes)
{
    const std::size_t new_size = round_up_to_page(new_bytes);
    if (new_size == r.size)
        return;

    void* p = ::mremap(r.base, r.size, new_size, MREMAP_MAYMOVE);
    if (p == MAP_FAILED)
        throw std::bad_alloc();
    r.base = p;
    r.size = new_size;
}

#endif  // __APPLE__
#endif  // ISOSPEC_WE_ARE_ON_WINDOWS

}  // namespace aligned_ptr_detail
#endif  // ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC


template<class T, std::size_t Alignment>
class aligned_unique_ptr {
    static_assert((Alignment & (Alignment - 1)) == 0); // meaning: Alignment is a power of two
    static_assert(Alignment >= alignof(T));
    static_assert(Alignment <= 4096, "Alignment must not exceed the smallest page size assumed across supported platforms");
    // Storage is acquired/released as raw bytes (aligned_alloc/free, or
    // mmap/VirtualAlloc) with no constructor/destructor ever run over it --
    // matches how this class has always behaved (it never invoked T's ctor
    // even before this rework, since it called the operator new[] *function*
    // directly rather than a `new T[n]` *expression*), just made explicit.
    static_assert(std::is_trivial_v<T>, "aligned_unique_ptr never runs T's constructor/destructor; T must be trivial");

    T* ptr_ = nullptr;
    std::size_t capacity_bytes_ = 0;

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
    bool is_mapped_ = false;
    aligned_ptr_detail::vm_region region_{};
#endif

    // C11 aligned_alloc requires size to be a multiple of alignment.
    static std::size_t round_up_to_alignment(std::size_t bytes) noexcept
    {
        return (bytes + Alignment - 1) / Alignment * Alignment;
    }

    // n * sizeof(T), then rounding up to Alignment, must not silently wrap:
    // that would turn a request for a huge buffer into a tiny allocation,
    // which the caller (having asked for n elements) would then write past.
    static std::size_t bytes_for(std::size_t n)
    {
        constexpr std::size_t max_size = (std::numeric_limits<std::size_t>::max)();
        if (n > max_size / sizeof(T))
            throw std::bad_alloc();
        const std::size_t bytes = n * sizeof(T);
        if (bytes > max_size - (Alignment - 1))
            throw std::bad_alloc();
        return bytes;
    }

    static T* allocate_small(std::size_t n)
    {
        if (n == 0)
            return nullptr;
        const std::size_t bytes = bytes_for(n);
#if ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE
        void* p = std::aligned_alloc(Alignment, round_up_to_alignment(bytes));
#else
        void* p = ::_aligned_malloc(bytes, Alignment);
#endif
        if (!p)
            throw std::bad_alloc();
        return static_cast<T*>(p);
    }

    static void free_small(void* p) noexcept
    {
#if ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE
        std::free(p);
#else
        ::_aligned_free(p);
#endif
    }

    void free_current() noexcept
    {
        if (!ptr_)
            return;
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        if (is_mapped_) {
            aligned_ptr_detail::vm_destroy(region_);
            return;
        }
#endif
        free_small(ptr_);
    }

    static void deleter_free(void* p, std::size_t) noexcept
    {
        std::free(p);
    }

    static void deleter_small(void* p, std::size_t) noexcept
    {
        free_small(p);
    }

public:
    static constexpr std::size_t alignment = Alignment;

    // Returned by release_with_deleter(): the raw pointer plus everything
    // needed to free it correctly, whichever backend actually produced it.
    // Unlike release(), never copies -- deleter(ptr, size) must be called
    // exactly once to avoid leaking it.
    struct release_result {
        T* ptr;
        std::size_t size;
        void (*deleter)(void*, std::size_t) noexcept;
    };

    T* allocate(std::size_t n)
    {
        return allocate_small(n);
    }

    explicit aligned_unique_ptr(std::size_t n)
        : ptr_(allocate(n))
        , capacity_bytes_(n * sizeof(T))
    {}

    aligned_unique_ptr() noexcept
        : ptr_(nullptr)
    {}

    ~aligned_unique_ptr()
    {
        free_current();
    }

    aligned_unique_ptr(const aligned_unique_ptr&) = delete;
    aligned_unique_ptr& operator=(const aligned_unique_ptr&) = delete;

    aligned_unique_ptr(aligned_unique_ptr&& other) noexcept
        : ptr_(std::exchange(other.ptr_, nullptr))
        , capacity_bytes_(std::exchange(other.capacity_bytes_, 0))
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        , is_mapped_(std::exchange(other.is_mapped_, false))
        , region_(std::exchange(other.region_, aligned_ptr_detail::vm_region{}))
#endif
    {}

    aligned_unique_ptr& operator=(aligned_unique_ptr&& other) noexcept
    {
        if (this != &other) {
            free_current();
            ptr_ = std::exchange(other.ptr_, nullptr);
            capacity_bytes_ = std::exchange(other.capacity_bytes_, 0);
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
            is_mapped_ = std::exchange(other.is_mapped_, false);
            region_ = std::exchange(other.region_, aligned_ptr_detail::vm_region{});
#endif
        }
        return *this;
    }

    T* get() noexcept
    {
        return std::assume_aligned<Alignment>(ptr_);
    }

    const T* get() const noexcept
    {
        return std::assume_aligned<Alignment>(ptr_);
    }

    T& operator[](std::size_t i) noexcept
    {
        return get()[i];
    }

    const T& operator[](std::size_t i) const noexcept
    {
        return get()[i];
    }

    explicit operator bool() const noexcept
    {
        return ptr_ != nullptr;
    }

    void reset(std::size_t n) noexcept
    {
        free_current();
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        is_mapped_ = false;
        region_ = aligned_ptr_detail::vm_region{};
#endif
        if (n > 0) {
            ptr_ = allocate_small(n);
            capacity_bytes_ = n * sizeof(T);
        } else {
            ptr_ = nullptr;
            capacity_bytes_ = 0;
        }
    }

    // Grows or shrinks the buffer to hold n elements, preserving the first
    // min(old capacity, new capacity) bytes. Below the page-size threshold
    // this is the same aligned_alloc/copy/free pattern reset() uses; once a
    // given instance has crossed that threshold it switches permanently to
    // a VM-backed region that the OS resizes in place (see the backends
    // above), so further growth/shrinkage stops needing to copy at all.
    void realloc(std::size_t n)
    {
        const std::size_t new_bytes = n == 0 ? 0 : bytes_for(n);

        if (new_bytes == 0) {
            free_current();
            ptr_ = nullptr;
            capacity_bytes_ = 0;
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
            is_mapped_ = false;
            region_ = aligned_ptr_detail::vm_region{};
#endif
            return;
        }

#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        if (is_mapped_) {
            aligned_ptr_detail::vm_resize(region_, new_bytes);
            ptr_ = static_cast<T*>(region_.base);
            capacity_bytes_ = new_bytes;
            return;
        }

        if (new_bytes >= aligned_ptr_detail::os_page_size()) {
            // Crossing the threshold: switch backend. This one transition
            // still has to copy -- there is nothing yet to remap.
            aligned_ptr_detail::vm_region r = aligned_ptr_detail::vm_create(new_bytes);
            if (ptr_) {
                std::memcpy(r.base, ptr_, std::min(capacity_bytes_, new_bytes));
                free_small(ptr_);
            }
            region_ = r;
            ptr_ = static_cast<T*>(region_.base);
            capacity_bytes_ = new_bytes;
            is_mapped_ = true;
            return;
        }
#endif

        T* new_ptr = allocate_small(n);
        if (ptr_) {
            std::memcpy(new_ptr, ptr_, std::min(capacity_bytes_, new_bytes));
            free_small(ptr_);
        }
        ptr_ = new_ptr;
        capacity_bytes_ = new_bytes;
    }

    // Hands back ownership as a plain T* that the caller can free with a
    // bare std::free() / C free(), regardless of which backend produced it.
    // Whenever the current backend isn't already free()-compatible on its
    // own -- a VM-backed region always needs this, and so does the small
    // backend on MSVC (see ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE
    // above) -- this materialises a fresh plain-malloc()'d copy instead.
    // That copy is the whole cost of calling this in those cases;
    // release_with_deleter() below avoids it when the caller can accept a
    // deleter instead of assuming free().
    T* release()
    {
        if (!ptr_) {
            capacity_bytes_ = 0;
            return nullptr;
        }

        bool needs_materialize = !ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE;
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        needs_materialize = needs_materialize || is_mapped_;
#endif
        if (needs_materialize) {
            void* copy = std::malloc(capacity_bytes_);
            if (!copy)
                throw std::bad_alloc();
            std::memcpy(copy, ptr_, capacity_bytes_);
            free_current();
            ptr_ = nullptr;
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
            is_mapped_ = false;
            region_ = aligned_ptr_detail::vm_region{};
#endif
            capacity_bytes_ = 0;
            return static_cast<T*>(copy);
        }

        capacity_bytes_ = 0;
        return std::assume_aligned<Alignment>(std::exchange(ptr_, nullptr));
    }

    // Fast counterpart to release(): never copies. Returns the raw pointer
    // together with the (size, deleter) pair needed to free it correctly --
    // the caller must invoke result.deleter(result.ptr, result.size) exactly
    // once (or not at all, if it hands the triple further on).
    release_result release_with_deleter() noexcept
    {
        T* p = std::exchange(ptr_, nullptr);
        std::size_t bytes = capacity_bytes_;
        capacity_bytes_ = 0;
#if ISOSPEC_ALIGNED_PTR_HAVE_VM_REALLOC
        if (is_mapped_) {
            is_mapped_ = false;
#if ISOSPEC_WE_ARE_ON_WINDOWS
            const std::size_t region_bytes = region_.reserved;
#else
            const std::size_t region_bytes = region_.size;
#endif
            region_ = aligned_ptr_detail::vm_region{};
            return release_result{p, region_bytes, &aligned_ptr_detail::vm_release_raw};
        }
#endif
#if ISOSPEC_ALIGNED_PTR_SMALL_BACKEND_FREE_COMPATIBLE
        return release_result{p, bytes, &deleter_free};
#else
        return release_result{p, bytes, &deleter_small};
#endif
    }
};

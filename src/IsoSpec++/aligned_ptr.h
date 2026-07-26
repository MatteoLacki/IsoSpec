#pragma once

#include <cstddef>      // std::size_t, std::nullptr_t
#include <cstdint>      // std::uintptr_t (only if doing runtime alignment checks)
#include <memory>       // std::assume_aligned (C++20)
#include <new>          // std::align_val_t, operator new[]
#include <utility>      // std::exchange
#include <type_traits>  // std::is_object_v (if keeping the static_assert)
#include <cassert>      // assert (if adding checks)


template<class T, std::size_t Alignment>
class aligned_unique_ptr {
    static_assert((Alignment & (Alignment - 1)) == 0); // meaning: Alignment is a power of two
    static_assert(Alignment >= alignof(T));

    T* ptr_ = nullptr;

public:
    static constexpr std::size_t alignment = Alignment;

    T* allocate(std::size_t n)
    {
        return static_cast<T*>(
            ::operator new[](n * sizeof(T),
                             std::align_val_t{Alignment}));
    }

    explicit aligned_unique_ptr(std::size_t n)
        : ptr_(allocate(n))
    {}

    aligned_unique_ptr() noexcept
        : ptr_(nullptr)
    {}

    ~aligned_unique_ptr()
    {
        if (ptr_)
        ::operator delete[](ptr_, std::align_val_t{Alignment});
    }

    aligned_unique_ptr(const aligned_unique_ptr&) = delete;
    aligned_unique_ptr& operator=(const aligned_unique_ptr&) = delete;

    aligned_unique_ptr(aligned_unique_ptr&& other) noexcept
        : ptr_(std::exchange(other.ptr_, nullptr))
    {}

    aligned_unique_ptr& operator=(aligned_unique_ptr&& other) noexcept
    {
        if (this != &other) {
            reset();
            ptr_ = std::exchange(other.ptr_, nullptr);
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
        if (ptr_) {
            ::operator delete[](ptr_, std::align_val_t{Alignment});
        }
        if (n > 0) {
            ptr_ = allocate(n);
        } else {
            ptr_ = nullptr;
        }
    }

    T* release() noexcept
    {
        return std::assume_aligned<Alignment>(std::exchange(ptr_, nullptr));
    }
};

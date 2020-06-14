#pragma once

#include <type_traits>
#include <cstdlib>

template<typename T> class pod_vector
{
    static_assert(std::is_trivially_copyable<T>::value, "Cannot use a pod_vector with a non-Plain Old Data type.");

    size_t backend_size;
    size_t first_free;
    T* store;

 public:

    pod_vector(size_t initial_size = 16) : backend_size(initial_size), first_free(0)
    {
        store = reinterpret_cast<T*>(malloc(sizeof(T) * initial_size));
        if(store == NULL)
            throw std::bad_alloc();
    }

    ~pod_vector() { free(store); };

    void fast_reserve(size_t n)
    {
        ISOSPEC_IMPOSSIBLE(n < backend_size);
        T* new_store = reinterpret_cast<T*>(reallocarray(store, n, sizeof(T)));
        if(new_store == NULL)
            throw std::bad_alloc();
        store = new_store;
        backend_size = n;
    }

    void reserve(size_t n)
    {
        if (n > backend_size)
            fast_reserve(n);
    }

    ISOSPEC_FORCE_INLINE void nocheck_push_back(const T& val) noexcept
    {
        ISOSPEC_IMPOSSIBLE(first_free >= backend_size);
        store[first_free] = val;
        first_free++;
    }

    ISOSPEC_FORCE_INLINE void push_back(const T& val)
    {
        if(first_free >= backend_size)
            fast_reserve(backend_size * 2);
        store[first_free] = val;
        first_free++;
    }

    ISOSPEC_FORCE_INLINE T& operator[](size_t n) noexcept
    {
        ISOSPEC_IMPOSSIBLE(n >= first_free);
        return store[n];
    }

    ISOSPEC_FORCE_INLINE const T& operator[](size_t n) const noexcept
    {
        ISOSPEC_IMPOSSIBLE(n >= first_free);
        return store[n];
    }

    ISOSPEC_FORCE_INLINE size_t size() const noexcept
    {
        return first_free;
    }

    ISOSPEC_FORCE_INLINE size_t capacity() const noexcept
    {
        return backend_size;
    }

    ISOSPEC_FORCE_INLINE T* data() noexcept
    {
        return store;
    }

    ISOSPEC_FORCE_INLINE const T* data() const noexcept
    {
        return store;
    }

    ISOSPEC_FORCE_INLINE bool empty() const noexcept
    {
        return first_free == 0;
    }

    ISOSPEC_FORCE_INLINE const T& back() const noexcept
    {
        ISOSPEC_IMPOSSIBLE(first_free > backend_size);
        return store[first_free-1];
    }

    ISOSPEC_FORCE_INLINE void pop_back() noexcept
    {
        // Unlike std::vector we do not ever shrink backend storage unless explicitly requested.
        ISOSPEC_IMPOSSIBLE(first_free == 0);
        first_free--;
    }

    void swap(pod_vector<T>& other) noexcept
    {
        std::swap(backend_size, other.backend_size);
        std::swap(first_free, other.first_free);
        std::swap(store, other.store);
    }
};
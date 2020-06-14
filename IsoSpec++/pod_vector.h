/*!
    Copyright (C) 2015-2020 Mateusz Łącki and Michał Startek.

    This file is part of IsoSpec.

    IsoSpec is free software: you can redistribute it and/or modify
    it under the terms of the Simplified ("2-clause") BSD licence.

    IsoSpec is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

    You should have received a copy of the Simplified BSD Licence
    along with IsoSpec.  If not, see <https://opensource.org/licenses/BSD-2-Clause>.
*/

#pragma once

#include <type_traits>
#include <cstdlib>
#include <utility>

template<typename T> class pod_vector
{
    static_assert(std::is_trivially_copyable<T>::value, "Cannot use a pod_vector with a non-Plain Old Data type.");

    T* backend_past_end;
    T* first_free;
    T* store;

 public:
    explicit pod_vector(size_t initial_size = 16)
    {
        store = reinterpret_cast<T*>(malloc(sizeof(T) * initial_size));
        if(store == NULL)
            throw std::bad_alloc();
        first_free = store;
        backend_past_end = store + initial_size;
    }

    ~pod_vector() { free(store); }

    void fast_reserve(size_t n)
    {
        ISOSPEC_IMPOSSIBLE(n < static_cast<size_t>(backend_past_end - store));
        T* new_store = reinterpret_cast<T*>(reallocarray(store, n, sizeof(T)));
        if(new_store == NULL)
            throw std::bad_alloc();
        first_free = new_store + (first_free - store);
        backend_past_end = new_store + n;
        store = new_store;
    }

    void reserve(size_t n)
    {
        if (n > backend_past_end - store)
            fast_reserve(n);
    }

    ISOSPEC_FORCE_INLINE void nocheck_push_back(const T& val) noexcept
    {
        ISOSPEC_IMPOSSIBLE(first_free >= backend_past_end);
        *first_free = val;
        first_free++;
    }

    ISOSPEC_FORCE_INLINE void push_back(const T& val)
    {
        if(first_free >= backend_past_end)
            fast_reserve((backend_past_end-store) * 2);
        *first_free = val;
        first_free++;
    }

    ISOSPEC_FORCE_INLINE T& operator[](size_t n) noexcept
    {
        ISOSPEC_IMPOSSIBLE(store + n >= first_free);
        return store[n];
    }

    ISOSPEC_FORCE_INLINE const T& operator[](size_t n) const noexcept
    {
        ISOSPEC_IMPOSSIBLE(store + n >= first_free);
        return store[n];
    }

    ISOSPEC_FORCE_INLINE size_t size() const noexcept
    {
        return first_free - store;
    }

    ISOSPEC_FORCE_INLINE size_t capacity() const noexcept
    {
        return backend_past_end - store;
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
        return first_free == store;
    }

    ISOSPEC_FORCE_INLINE const T& back() const noexcept
    {
        ISOSPEC_IMPOSSIBLE(first_free > backend_past_end);
        return *(first_free-1);
    }

    ISOSPEC_FORCE_INLINE void pop_back() noexcept
    {
        // Unlike std::vector we do not ever shrink backend storage unless explicitly requested.
        ISOSPEC_IMPOSSIBLE(first_free == store);
        first_free--;
    }

    void swap(pod_vector<T>& other) noexcept
    {
        std::swap(backend_past_end, other.backend_past_end);
        std::swap(first_free, other.first_free);
        std::swap(store, other.store);
    }

    typedef T* iterator;
    typedef T value_type;
    typedef size_t size_type;
    typedef T& reference;
    typedef const T& const_reference;

    iterator begin()
    {
        return store;
    }

    iterator end()
    {
        return first_free;
    }

    ISOSPEC_FORCE_INLINE const T& front() const noexcept
    {
        ISOSPEC_IMPOSSIBLE(store == first_free);
        return *store;
    }
};

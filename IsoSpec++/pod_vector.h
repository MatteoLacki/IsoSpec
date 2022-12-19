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
#include <cstddef>
#include <utility>
#include <new>
#include <algorithm>
#include "platform.h"
#include "conf.h"



template<typename T> class unsafe_pod_vector;

template<typename T> class pod_vector
{
    T* backend_past_end;
    T* first_free;
    T* store;

 public:
    explicit pod_vector(size_t initial_size = 16)
    {
    #if !ISOSPEC_BUILDING_R
        static_assert(std::is_trivially_copyable<T>::value, "Cannot use a pod_vector with a non-Plain Old Data type.");
    #endif

        store = reinterpret_cast<T*>(malloc(sizeof(T) * initial_size));
        if(store == NULL)
            throw std::bad_alloc();
        first_free = store;
        backend_past_end = store + initial_size;
    }

    pod_vector(const pod_vector<T>& other) = delete;
    pod_vector& operator=(const pod_vector<T>& other) = delete;
    pod_vector& operator=(pod_vector<T>&& other)
    {
        free(store);
        backend_past_end = other.backend_past_end;
        first_free = other.first_free;
        store = other.store;
        other.backend_past_end = other.first_free = other.store = NULL;
        return *this;
    }

    pod_vector(pod_vector<T>&& other)
    {
        backend_past_end = other.backend_past_end;
        first_free = other.first_free;
        store = other.store;
        other.backend_past_end = other.first_free = other.store = NULL;
    }

    ~pod_vector() { free(store); backend_past_end = first_free = store = NULL; }

    explicit pod_vector(unsafe_pod_vector<T>&& other)
    {
        backend_past_end = other.backend_past_end;
        first_free = other.first_free;
        store = other.store;
       other.backend_past_end = other.first_free = other.store = NULL;
    }

    void fast_reserve(size_t n)
    {
        ISOSPEC_IMPOSSIBLE(n < static_cast<size_t>(backend_past_end - store));
        const std::ptrdiff_t store_used_size = first_free - store;
        T* new_store = reinterpret_cast<T*>(realloc(store, n * sizeof(T)));
        if(new_store == NULL)
            throw std::bad_alloc();
        first_free = new_store + store_used_size;
        backend_past_end = new_store + n;
        store = new_store;
    }

    void reserve(size_t n)
    {
        if (n > static_cast<size_t>(backend_past_end - store))
            fast_reserve(n);
    }

    void resize(size_t new_size)
    {
        ISOSPEC_IMPOSSIBLE(static_cast<std::ptrdiff_t>(new_size) < first_free - store);
        size_t cap = capacity();
        if(cap < new_size)
        {
            do {
            cap = cap * 2;
            } while(cap < new_size);
            fast_reserve(cap);
        }
        first_free = store + new_size;
    }

    void resize_and_wipe(size_t new_size)
    {
        size_t old_size = size();
        ISOSPEC_IMPOSSIBLE(new_size <= old_size);
        resize(new_size);
        memset(store+old_size, 0, (new_size-old_size) * sizeof(T));
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
            fast_reserve((std::max<std::ptrdiff_t>)(4, (backend_past_end-store)) * 2);
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
    typedef const T* const_iterator;
    typedef T value_type;
    typedef size_t size_type;
    typedef T& reference;
    typedef const T& const_reference;

    iterator begin() noexcept { return store; };
    const_iterator begin() const noexcept { return store; }
    const_iterator cbegin() const noexcept { return store; }
    iterator end() noexcept { return first_free; }
    const_iterator end() const noexcept { return first_free; }
    const_iterator cend() const noexcept { return first_free; }

    ISOSPEC_FORCE_INLINE const T& front() const noexcept
    {
        ISOSPEC_IMPOSSIBLE(store == first_free);
        return *store;
    }

    void clear()
    {
        free(store);
        first_free = store = backend_past_end = NULL;
    }

    friend class unsafe_pod_vector<T>;
};


template<typename T> class unsafe_pod_vector
{

    T* backend_past_end;
    T* first_free;
    T* store;

 public:
    unsafe_pod_vector() = default;

    void init() {
    #if !ISOSPEC_BUILDING_R
        static_assert(std::is_trivially_copyable<T>::value, "Cannot use a pod_vector with a non-Plain Old Data type.");
        static_assert(std::is_trivially_copyable<unsafe_pod_vector<T> >::value, "Cannot use a pod_vector with a non-Plain Old Data type.");
    #endif
        memset(this, 0, sizeof(*this));
    }

    void init(size_t initial_size)
    {
    #if !ISOSPEC_BUILDING_R
        static_assert(std::is_trivially_copyable<T>::value, "Cannot use a pod_vector with a non-Plain Old Data type.");
        static_assert(std::is_trivially_copyable<unsafe_pod_vector<T> >::value, "Cannot use a pod_vector with a non-Plain Old Data type.");
    #endif
        store = reinterpret_cast<T*>(malloc(sizeof(T) * initial_size));
        if(store == NULL)
            throw std::bad_alloc();
        first_free = store;
        backend_past_end = store + initial_size;
    }

    unsafe_pod_vector(const pod_vector<T>& other) = delete;  // NOLINT(runtime/explicit) - seriously? Deleted constructors have to be marked explicit?
    unsafe_pod_vector& operator=(const pod_vector<T>& other) = delete;
    //unsafe_pod_vector(unsafe_pod_vector<T>&& other) = default;

    ~unsafe_pod_vector() = default;

    void fast_reserve(size_t n)
    {
        ISOSPEC_IMPOSSIBLE(n < static_cast<size_t>(backend_past_end - store));
        T* new_store = reinterpret_cast<T*>(realloc(store, n * sizeof(T)));
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

    void resize(size_t new_size)
    {
        ISOSPEC_IMPOSSIBLE(new_size < first_free - store);
        size_t cap = capacity();
        if(cap < new_size)
        {
            do {
            cap = cap * 2;
            } while(cap < new_size);
            fast_reserve(cap);
        }
        first_free = store + new_size;
    }

    void resize_and_wipe(size_t new_size)
    {
        size_t old_size = size();
        ISOSPEC_IMPOSSIBLE(new_size <= old_size);
        resize(new_size);
        memset(store+old_size, 0, (new_size-old_size) * sizeof(T));
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
            fast_reserve((std::max<std::ptrdiff_t>)(4, (backend_past_end-store)) * 2);
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
    typedef const T* const_iterator;
    typedef T value_type;
    typedef size_t size_type;
    typedef T& reference;
    typedef const T& const_reference;

    iterator begin() noexcept { return store; };
    const_iterator begin() const noexcept { return store; }
    const_iterator cbegin() const noexcept { return store; }
    iterator end() noexcept { return first_free; }
    const_iterator end() const noexcept { return first_free; }
    const_iterator cend() const noexcept { return first_free; }

    ISOSPEC_FORCE_INLINE const T& front() const noexcept
    {
        ISOSPEC_IMPOSSIBLE(store == first_free);
        return *store;
    }

    void clear()
    {
        free(store);
        first_free = store = backend_past_end = NULL;
    }

    friend class pod_vector<T>;
};

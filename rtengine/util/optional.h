/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2026 Daniel Gao <daniel.gao.work@gmail.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <stdexcept>
#include <type_traits>
#include <utility>

namespace rt {

class bad_optional_access : public std::logic_error
{
public:
    explicit bad_optional_access(const std::string& msg) : std::logic_error(msg) {}
};

struct nullopt_t { explicit nullopt_t() = default; };
static constexpr nullopt_t nullopt{};

// Minimal backport of C++17 std::optional
template <class T>
class optional
{
public:
    constexpr optional() : m_has_value(false) {}
    constexpr optional(nullopt_t) noexcept : m_has_value(false) {}
    optional(const T& value) : m_has_value(true) { new(&m_storage) T(value); }
    optional(T&& value) : m_has_value(true) { new(&m_storage) T(std::move(value)); }

    // Copy constructor/assignment
    optional(const optional& other) : m_has_value(other.m_has_value)
    {
        if (m_has_value) {
            new(&m_storage) T(*other);
        }
    }
    optional& operator=(const optional& other)
    {
        if (this != &other) {
            reset();
            if (other) {
                new(&m_storage) T(*other);
                m_has_value = true;
            }
        }
        return *this;
    }

    // Move constructor/assignment
    optional(optional&& other)
        noexcept(std::is_nothrow_move_constructible<T>::value)
        : m_has_value(other.m_has_value)
    {
        if (m_has_value) {
            new(&m_storage) T(std::move(*other));
            other.reset();
        }
    }
    optional& operator=(optional&& other)
        noexcept(std::is_nothrow_move_assignable<T>::value
                 && std::is_nothrow_move_constructible<T>::value)
    {
        if (this != &other) {
            reset();
            if (other) {
                new(&m_storage) T(std::move(*other));
                m_has_value = true;
                other.reset();
            }
        }
        return *this;
    }

    ~optional() { reset(); }

    explicit operator bool() const { return m_has_value; }
    bool has_value() const { return m_has_value; }

    T& value()
    {
        if (!m_has_value) throw bad_optional_access("bad optional access");
        return *ptr();
    }
    const T& value() const
    {
        if (!m_has_value) throw bad_optional_access("bad optional access");
        return *ptr();
    }

    T& operator*() { return *ptr(); }
    const T& operator*() const { return *ptr(); }
    T* operator->() { return  ptr(); }
    const T* operator->() const { return ptr(); }

    template <class U>
    T value_or(U&& fallback) const
    {
        return m_has_value ? *ptr() : static_cast<T>(std::forward<U>(fallback));
    }

    void reset()
    {
        if (m_has_value) {
            ptr()->~T();
            m_has_value = false;
        }
    }

    template <class... Args>
    T& emplace(Args&&... args)
    {
        reset();
        new(&m_storage) T(std::forward<Args>(args)...);
        m_has_value = true;
        return *ptr();
    }

    void swap(optional& other)
        noexcept(std::is_nothrow_move_constructible<T>::value
                 && noexcept(std::swap(std::declval<T&>(), std::declval<T&>())))
    {
        if (*this && other) {
            using std::swap;
            swap(**this, *other);
        } else if (*this && !other) {
            other = std::move(*this);
            reset();
        } else if (!*this && other) {
            *this = std::move(other);
            other.reset();
        }
    }

private:
    T* ptr() { return reinterpret_cast<T*>(&m_storage); }
    const T* ptr() const { return reinterpret_cast<const T*>(&m_storage); }

    alignas(T) unsigned char m_storage[sizeof(T)];
    bool m_has_value;
};

template <class T>
bool operator==(const optional<T>& a, const optional<T>& b)
{
    if (a.has_value() != b.has_value()) return false;
    else if (!a.has_value()) return true;
    else return *a == *b;
}
template <class T>
bool operator!=(const optional<T>& a, const optional<T>& b) { return !(a == b); }

template <typename T>
void swap(optional<T>& a, optional<T>& b) noexcept(noexcept(a.swap(b))) { a.swap(b); }

template <class T>
optional<T> make_optional(T&& value) { return optional<T>(std::forward<T>(value)); }

template <class T, class... Args>
optional<T> make_optional(Args&&... args)
{
    optional<T> o;
    o.emplace(std::forward<Args>(args)...);
    return o;
}

}  // namespace rt

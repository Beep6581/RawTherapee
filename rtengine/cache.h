/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (C) 2016 Flössie <floessie.mail@gmail.com>
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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <algorithm>
#include <list>
#include <map>
#include <memory>
#include <type_traits>
#include <unordered_map>

#include "../rtgui/threadutils.h"

namespace rtengine
{

namespace cache_helper
{

    // See http://stackoverflow.com/a/20790050
    template<typename, typename = void>
    struct has_hash
        : std::false_type
    {
    };

    template<typename T>
    struct has_hash<T, decltype(std::hash<T>()(std::declval<T>()), void())>
        : std::true_type
    {
    };

}

template<class K, class V>
class Cache
{
public:
    class Hook
    {
    public:
        virtual ~Hook()
        {
        }
        virtual void onDiscard(const K& key, const V& value) = 0;
        virtual void onDisplace(const K& key, const V& value) = 0;
        virtual void onRemove(const K& key, const V& value) = 0;
        virtual void onDestroy() = 0;
    };

    Cache(unsigned long _size, Hook* _hook = nullptr) :
        store_size(_size),
        hook(_hook)
    {
    }

    ~Cache()
    {
        if (hook) {
            resize(0);
            hook->onDestroy();
        }
    }

    bool get(const K& key, V& value) const
    {
        mutex.lock();
        const StoreConstIterator store_it = store.find(key);
        const bool present = store_it != store.end();
        if (present) {
            lru_list.splice(
                lru_list.begin(),
                lru_list,
                store_it->second->lru_list_it
            );
            value = store_it->second->value;
        }
        mutex.unlock();

        return present;
    }

    bool set(const K& key, const V& value)
    {
        return set(key, value, Mode::UNCOND);
    }

    bool replace(const K& key, const V& value)
    {
        return set(key, value, Mode::KNOWN);
    }

    bool insert(const K& key, const V& value)
    {
        return set(key, value, Mode::UNKNOWN);
    }

    bool remove(const K& key)
    {
        mutex.lock();
        const StoreIterator store_it = store.find(key);
        const bool present = store_it != store.end();
        if (present) {
            remove(store_it);
        }
        mutex.unlock();

        return present;
    }

    void resize(unsigned long size)
    {
        mutex.lock();
        while (lru_list.size() > size) {
            discard();
        }
        store_size = size;
        mutex.unlock();
    }

    void clear()
    {
        mutex.lock();
        if (hook) {
            for (const auto& entry : store) {
                hook->onRemove(entry.first, entry.second->value);
            }
        }
        lru_list.clear();
        store.clear();
        mutex.unlock();
    }

private:
    struct Value;

    using Store = typename std::conditional<
        cache_helper::has_hash<K>::value,
        std::unordered_map<K, std::unique_ptr<Value>>,
        std::map<K, std::unique_ptr<Value>>
    >::type;
    using StoreIterator = typename Store::iterator;
    using StoreConstIterator = typename Store::const_iterator;

    using LruList = std::list<StoreIterator>;
    using LruListIterator = typename LruList::iterator;

    struct Value {
        V value;
        LruListIterator lru_list_it;
    };

    enum class Mode {
        UNCOND,
        KNOWN,
        UNKNOWN
    };

    void discard()
    {
        const StoreIterator store_it = lru_list.back();
        if (hook) {
            hook->onDiscard(store_it->first, store_it->second->value);
        }
        store.erase(store_it);
        lru_list.pop_back();
    }

    bool set(const K& key, const V& value, Mode mode)
    {
        mutex.lock();
        const StoreIterator store_it = store.find(key);
        const bool is_new_key = store_it == store.end();
        if (is_new_key) {
            if (mode == Mode::UNCOND || mode == Mode::UNKNOWN) {
                if (lru_list.size() >= store_size) {
                    discard();
                }
                lru_list.push_front(store.end());
                std::unique_ptr<Value> v(
                    new Value{
                        value,
                        lru_list.begin()
                    }
                );
                lru_list.front() = store.emplace(key, std::move(v)).first;
            }
        } else {
            if (mode == Mode::UNCOND || mode == Mode::KNOWN) {
                if (hook) {
                    hook->onDisplace(key, store_it->second->value);
                }
                lru_list.splice(
                    lru_list.begin(),
                    lru_list,
                    store_it->second->lru_list_it
                );
                store_it->second->value = value;
            }
        }
        mutex.unlock();

        return is_new_key;
    }

    void remove(const StoreIterator& store_it)
    {
        if (hook) {
            hook->onRemove(store_it->first, store_it->second->value);
        }
        lru_list.erase(store_it->second->lru_list_it);
        store.erase(store_it);
    }

    unsigned long store_size;
    Hook* const hook;
    mutable MyMutex mutex;
    Store store;
    mutable LruList lru_list;
};

}

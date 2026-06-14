/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2025 Daniel Gao <daniel.gao.work@gmail.com>
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

#include "rtengine/curves.h"

#include <glibmm/arrayhandle.h>
#include <glibmm/keyfile.h>
#include <glibmm/ustring.h>

#include <map>
#include <string>
#include <type_traits>
#include <vector>

Glib::ustring expandRelativePath(const Glib::ustring &procparams_fname,
                                 const Glib::ustring &prefix,
                                 Glib::ustring embedded_fname);

Glib::ustring expandRelativePath2(const Glib::ustring &procparams_fname,
                                  const Glib::ustring &procparams_fname2,
                                  const Glib::ustring &prefix,
                                  Glib::ustring embedded_fname);

Glib::ustring relativePathIfInside(const Glib::ustring &procparams_fname,
                                   bool fnameAbsolute,
                                   Glib::ustring embedded_fname);

Glib::ustring relativePathIfInside2(const Glib::ustring &procparams_fname,
                                    const Glib::ustring &procparams_fname2,
                                    bool fnameAbsolute,
                                    Glib::ustring embedded_fname);

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int& value
)
{
    value = keyfile.get_integer(group_name, key);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double& value
)
{
    value = keyfile.get_double(group_name, key);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    float& value
)
{
    value = static_cast<float>(keyfile.get_double(group_name, key));
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool& value
)
{
    value = keyfile.get_boolean(group_name, key);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    Glib::ustring& value
)
{
    value = keyfile.get_string(group_name, key);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<int>& value
)
{
    value = keyfile.get_integer_list(group_name, key);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<double>& value
)
{
    value = keyfile.get_double_list(group_name, key);
    rtengine::sanitizeCurve(value);
}

inline void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<std::string>& value
)
{
    auto tmpval = keyfile.get_string_list(group_name, key);
    value.assign(tmpval.begin(), tmpval.end());
}

template<typename T>
bool assignFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    T& value,
    bool& params_edited_value
)
{
    if (keyfile.has_key(group_name, key)) {
        getFromKeyfile(keyfile, group_name, key, value);

        params_edited_value = true;

        return true;
    }

    return false;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool assignFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::map<std::string, T>& mapping,
    T& value,
    bool& params_edited_value
)
{
    if (keyfile.has_key(group_name, key)) {
        Glib::ustring v;
        getFromKeyfile(keyfile, group_name, key, v);

        const typename std::map<std::string, T>::const_iterator m = mapping.find(v);

        if (m != mapping.end()) {
            value = m->second;
        } else {
            return false;
        }

        params_edited_value = true;

        return true;
    }

    return false;
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_integer(group_name, key, value);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    float value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_double(group_name, key, static_cast<double>(value));
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_double(group_name, key, value);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_boolean(group_name, key, value);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const Glib::ustring& value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_string(group_name, key, value);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<int>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<int> list = value;
    keyfile.set_integer_list(group_name, key, list);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<double>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<double> list = value;
    keyfile.set_double_list(group_name, key, list);
}

inline void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<std::string>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<Glib::ustring> list = value;
    keyfile.set_string_list(group_name, key, list);
}

template<typename T>
bool saveToKeyfile(
    bool save,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    if (save) {
        putToKeyfile(group_name, key, value, keyfile);
        return true;
    }

    return false;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool saveToKeyfile(
    bool save,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::map<T, const char*>& mapping,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    if (save) {
        const typename std::map<T, const char*>::const_iterator m = mapping.find(value);

        if (m != mapping.end()) {
            keyfile.set_string(group_name, key, m->second);
            return true;
        }
    }

    return false;
}

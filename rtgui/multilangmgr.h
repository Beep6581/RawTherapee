/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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

#include <map>
#include <string>
#include <vector>

#include <glibmm/ustring.h>

class TranslationMetadata
{
public:
    TranslationMetadata() = default;
    ~TranslationMetadata() = default;
    TranslationMetadata(const TranslationMetadata &other) = delete;
    TranslationMetadata(TranslationMetadata &&other) = delete;
    explicit TranslationMetadata(std::map<std::string, std::string> &&metadata);

    TranslationMetadata &operator =(const TranslationMetadata &other) = delete;
    TranslationMetadata &operator =(TranslationMetadata &&other) noexcept = default;

    std::string get(const std::string &key, const std::string &default_value) const;
    std::string getLanguageName(const std::string &default_name) const;

private:
    std::map<std::string, std::string> metadata;
};

class MultiLangMgr
{
public:
    MultiLangMgr ();

    void load(const Glib::ustring &language, const std::vector<Glib::ustring> &fnames);
    Glib::ustring getStr(const std::string& key) const;
    const TranslationMetadata *getMetadata(const Glib::ustring &fname) const;
    static bool isOSLanguageDetectSupported();
    static Glib::ustring getOSUserLanguage();

private:
    std::map<std::string, Glib::ustring> translations;
    mutable std::map<Glib::ustring, TranslationMetadata> lang_files_metadata;
};

extern MultiLangMgr langMgr;

inline Glib::ustring M (const std::string& key)
{
    return langMgr.getStr (key);
}

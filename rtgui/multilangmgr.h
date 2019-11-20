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

class MultiLangMgr
{
public:
    MultiLangMgr ();

    void load(const Glib::ustring &language, const std::vector<Glib::ustring> &fnames);
    Glib::ustring getStr(const std::string& key) const;
    static bool isOSLanguageDetectSupported();
    static Glib::ustring getOSUserLanguage();

private:
    std::map<std::string, Glib::ustring> translations;
};

extern MultiLangMgr langMgr;

inline Glib::ustring M (const std::string& key)
{
    return langMgr.getStr (key);
}

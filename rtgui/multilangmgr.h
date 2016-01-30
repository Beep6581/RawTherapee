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
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _MULTILANGMGR_
#define _MULTILANGMGR_

#include <map>
#include <memory>
#include <string>

#include <glibmm/ustring.h>

class MultiLangMgr
{
public:
    MultiLangMgr ();
    MultiLangMgr (const Glib::ustring& fname, MultiLangMgr* fallbackMgr = nullptr);

public:
    bool load (const Glib::ustring& fname, MultiLangMgr* fallbackMgr = nullptr);

public:
    Glib::ustring getStr (const std::string& key) const;

public:
    static bool isOSLanguageDetectSupported ();
    static Glib::ustring getOSUserLanguage ();

private:
    std::map<std::string, Glib::ustring> translations;
    std::unique_ptr<MultiLangMgr> fallbackMgr;

};

extern MultiLangMgr langMgr;

inline Glib::ustring M (const std::string& key)
{
    return langMgr.getStr (key);
}

#endif

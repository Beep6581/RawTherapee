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
#include <string>
#include <glibmm.h>

class MultiLangMgr {

    std::map<std::string, Glib::ustring> transTable;
    MultiLangMgr* fallBack;
    
    Glib::ustring TranslateRFC2Language(Glib::ustring rfcName);

    public: 
    MultiLangMgr () : fallBack (NULL) {}
    MultiLangMgr (Glib::ustring fname) : fallBack (NULL) { load (fname); }
    MultiLangMgr (Glib::ustring fname, MultiLangMgr* fb) : fallBack (NULL) { load (fname, fb); }

    bool load (Glib::ustring fname, MultiLangMgr* fb = NULL);
    bool save (Glib::ustring fname);

    bool isOSLanguageDetectSupported();
    Glib::ustring getOSUserLanguage();
    Glib::ustring getStr (std::string key);
};

extern MultiLangMgr langMgr;

Glib::ustring M (std::string key);

#endif

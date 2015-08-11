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
#ifdef WIN32
// Desired auto detect function is Vista+
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 8) || __GNUC__ > 4
#define WINVER 0x0600 // switching to WINVER for gcc 4.8.1 support on Winx64
#else
#define _WIN32_WINNT 0x0600
#endif
#include <windows.h>
#include <winnls.h>
#if (__GNUC__ == 4 && __GNUC_MINOR__ >= 8) || __GNUC__ > 4
#undef WINVER
#else
#undef _WIN32_WINNT
#endif
#endif
#include <glib/gstdio.h>
#include "multilangmgr.h"
#include <cstring>
#include "../rtengine/safegtk.h"

MultiLangMgr langMgr;

Glib::ustring M (std::string key)
{
    return langMgr.getStr (key);
}

// fb is fallback manager if the first could not be loaded
bool MultiLangMgr::load (Glib::ustring fname, MultiLangMgr* fb)
{
    FILE *f = safe_g_fopen (fname, "rt");

    fallBack = fb;

    if (f == NULL) {
        return false;
    }

    transTable.clear ();

    char* buffer = new char[2048];

    while (fgets (buffer, 2048, f) != 0) {
        // find separator
        int seppos = 0;

        while (buffer[seppos] != 0 && buffer[seppos] != ';') {
            seppos++;
        }

        // no separator found
        if (buffer[seppos] == 0) {
            continue;
        }

        // cut the last \n and \r characters
        int endpos = strlen(buffer) - 1;

        while (buffer[endpos] == '\n' || buffer[endpos] == '\r') {
            endpos--;
        }

        buffer[endpos + 1] = 0;
        // replace "\n" to '\n'
        int j = 0;

        for (int i = 0; i < endpos + 1; i++)
            if (i < endpos && buffer[i] == '\\' && buffer[i + 1] == 'n') {
                buffer[j++] = '\n';
                i++;
            } else {
                buffer[j++] = buffer[i];
            }

        buffer[j] = 0;
        // cut to two parts
        buffer[seppos] = 0;
        transTable[buffer] = buffer + seppos + 1;
    }

    fclose (f);
    delete [] buffer;
    return true;
}

bool MultiLangMgr::save (Glib::ustring fname)
{

    FILE *f = safe_g_fopen (fname, "wt");

    if (f == NULL) {
        return false;
    }

    std::map<std::string, Glib::ustring>::iterator r;

    for (r = transTable.begin (); r != transTable.end(); r++) {
        fprintf (f, "%s;%s\n", r->first.c_str(), safe_locale_from_utf8(r->second).c_str());
    }

    fclose (f);
    return true;
}


bool MultiLangMgr::isOSLanguageDetectSupported()
{
#ifdef WIN32
#ifdef __MINGW64_VERSION_MAJOR
    // Only on Vista or above
    return LOBYTE(LOWORD(GetVersion())) >= 6;
#else
    return false;
#endif
#elif defined(__linux__) || defined(__APPLE__)
    return true;
#else
    return false;
#endif
}


// returns Language name mapped from the currently selected OS language
Glib::ustring MultiLangMgr::getOSUserLanguage()
{
    Glib::ustring langName = Glib::ustring("default");

    if (isOSLanguageDetectSupported()) {

#ifdef WIN32
// When using old versions of MINGW this is not defined
#ifdef __MINGW64_VERSION_MAJOR
        WCHAR langRFCU[64] = {0};

        if (GetUserDefaultLocaleName(langRFCU, 64) != 0 && lstrlenW(langRFCU) >= 2) {
            // convert UNICODE16 to GTK
            char langRFC[64];
            WideCharToMultiByte(CP_UTF8, 0, langRFCU, -1, langRFC, 64, 0, 0);
            Glib::ustring localRFC = Glib::ustring(langRFC);

            langName = TranslateRFC2Language(localRFC);
        }

#endif
#elif defined(__linux__) || defined(__APPLE__)
        char *tmplocale;
        tmplocale = setlocale(LC_CTYPE, "");
        setlocale(LC_NUMERIC, "C"); // to set decimal point to "."

        if(tmplocale) {
            langName = TranslateRFC2Language(tmplocale);
        }

#endif
    }

    return langName;
}

// Translates RFC standard language code to file name, e.g. "de-DE" to "Deutsch"
Glib::ustring MultiLangMgr::TranslateRFC2Language(Glib::ustring rfcName)
{
    if (rfcName.length() < 2) {
        return Glib::ustring("default");
    }

    Glib::ustring major = rfcName.substr(0, 2).lowercase();
    Glib::ustring minor;

    if (rfcName.length() >= 5) {
        minor = rfcName.substr(3, 2).uppercase();
    }

    //printf("Lang: %s - %s\n",major.c_str(),minor.c_str());

    if (major == "ca") {
        return "Catala";
    }

    if (major == "zh") {
        return (minor == "CN" || minor == "SG") ? "Chinese (Simplified)" : "Chinese (Traditional)";
    }

    if (major == "cs") {
        return "Czech";
    }

    if (major == "da") {
        return "Dansk";
    }

    if (major == "de") {
        return "Deutsch";
    }

    if (major == "es") {
        return "Espanol";
    }

    if (major == "eu") {
        return "Euskara";
    }

    if (major == "fr") {
        return "Francais";
    }

    if (major == "el") {
        return "Greek";
    }

    if (major == "he") {
        return "Hebrew";
    }

    if (major == "it") {
        return "Italiano";
    }

    if (major == "ja") {
        return "Japanese";
    }

    if (major == "lv") {
        return "Latvian";
    }

    if (major == "hu") {
        return "Magyar";
    }

    if (major == "nl") {
        return "Nederlands";
    }

    if (major == "nn" || major == "nb") {
        return "Norsk BM";
    }

    if (major == "pl") {
        return "Polish";
    }

    if (major == "pt") {
        return "Portugues (Brasil)";
    }

    if (major == "ru") {
        return "Russian";
    }

    if (major == "sr") {
        return "Serbian (Cyrilic Characters)";
    }

    if (major == "sk") {
        return "Slovak";
    }

    if (major == "fi") {
        return "Suomi";
    }

    if (major == "sv") {
        return "Swedish";
    }

    if (major == "tr") {
        return "Turkish";
    }

    // Don't split en-US, en-GB, etc. since only default english is constantly updated
    return "default";
}

Glib::ustring MultiLangMgr::getStr (std::string key)
{

    std::map<std::string, Glib::ustring>::iterator r = transTable.find (key);

    if (r != transTable.end()) {
        return r->second;
    } else if (fallBack) {
        return fallBack->getStr (key);
    } else {
        return key;
    }
}

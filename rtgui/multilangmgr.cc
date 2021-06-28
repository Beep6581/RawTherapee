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
#include "multilangmgr.h"

#include <fstream>
#include <glib.h>
#ifdef WIN32
#include <windows.h>
#include <winnls.h>
#endif
#ifdef __APPLE__
#include <CoreFoundation/CoreFoundation.h>
#endif

namespace
{

// Maps standard locales to languages, e.g. "de-DE" to "Deutsch".
struct LocaleToLang : private std::map<std::pair<Glib::ustring, Glib::ustring>, Glib::ustring>
{
    static const std::pair<Glib::ustring, Glib::ustring> key (const Glib::ustring& major, const Glib::ustring& minor = Glib::ustring ())
    {
        return std::make_pair (major, minor);
    }

    LocaleToLang ()
    {
        emplace (key ("ca", "ES"), "Catala");
        emplace (key ("cs", "CZ"), "Czech");
        emplace (key ("da", "DK"), "Dansk");
        emplace (key ("de", "DE"), "Deutsch");
#ifdef __APPLE__
        emplace (key ("en", "UK"), "English (UK)");
#else
        emplace (key ("en", "GB"), "English (UK)");
#endif
        emplace (key ("en", "US"), "English (US)");
        emplace (key ("es", "ES"), "Espanol");
        emplace (key ("eu", "ES"), "Euskara");
        emplace (key ("fr", "FR"), "Francais");
        emplace (key ("el", "GR"), "Greek");
        emplace (key ("he", "IL"), "Hebrew");
        emplace (key ("it", "IT"), "Italiano");
        emplace (key ("ja", "JP"), "Japanese");
        emplace (key ("lv", "LV"), "Latvian");
        emplace (key ("hu", "HU"), "Magyar");
        emplace (key ("nl", "NL"), "Nederlands");
        emplace (key ("nn", "NO"), "Norsk BM");
        emplace (key ("nb", "NO"), "Norsk BM");
        emplace (key ("pl", "PL"), "Polish");
        emplace (key ("pt", "PT"), "Portugues (Brasil)");
        emplace (key ("ru", "RU"), "Russian");
        emplace (key ("sr", "RS"), "Serbian (Cyrilic Characters)");
        emplace (key ("sk", "SK"), "Slovak");
        emplace (key ("fi", "FI"), "Suomi");
        emplace (key ("sv", "SE"), "Swedish");
        emplace (key ("tr", "TR"), "Turkish");
        emplace (key ("zh", "CN"), "Chinese (Simplified)");
        emplace (key ("zh", "SG"), "Chinese (Traditional)");
    }

    Glib::ustring operator() (const Glib::ustring& locale) const
    {
        Glib::ustring major, minor;

        if (locale.length () >= 2) {
            major = locale.substr (0, 2).lowercase ();
        }

        if (locale.length () >= 5) {
            minor = locale.substr (3, 2).uppercase ();
        }

        // Look for matching language and country.
        auto iterator = find (key (major, minor));

        if (iterator != end ()) {
            return iterator->second;
        }

        // Look for matching language only.
        iterator = find (key (major, major.uppercase()));

        if (iterator != end ()) {
            return iterator->second;
        }

        return "default";
    }

    std::string getLocale(const Glib::ustring &language) const
    {
        for (auto &p : *this) {
            if (p.second == language) {
                std::string ret = p.first.first;
                if (!p.first.second.empty()) {
                    ret += "_" + p.first.second;
                }
                return ret;
            }
        }
        return "C";
    }
};

const LocaleToLang localeToLang;

void setGtkLanguage(const Glib::ustring &language)
{
    if(language != "default") { // nothing to change when using default
        std::string lang = localeToLang.getLocale(language);
#ifdef __APPLE__

        // On MacOS, LANG environment variable is not defined when running app bundle
        // So we should set all locale data
        const Glib::ustring localeUTF8 = lang + ".UTF-8";

        lang = lang + ".UTF-8"; // According to Apple documentation, UTF-8 is a built-in encoding on all platforms on which macOS runs

        g_setenv("LANG", lang.c_str(), true);
        setlocale(LC_ALL, lang.c_str());
        setlocale (LC_NUMERIC, "C"); // Force decimal point to dot.

#else

        const gchar *env_langc = g_getenv("LANG");
        if(env_langc) {
            const std::string env_lang(env_langc);
            if (!env_lang.empty()) {
                const std::string::size_type suffix_pos = env_lang.find_first_of(".");
                if (suffix_pos != std::string::npos) {
                    lang += env_lang.substr(suffix_pos);
                }
            }
        }

        g_setenv("LANG", lang.c_str(), true);
        
#endif
    }
}

}

MultiLangMgr langMgr;

MultiLangMgr::MultiLangMgr ()
{
}

void MultiLangMgr::load(const Glib::ustring &language, const std::vector<Glib::ustring> &fnames)
{
    setGtkLanguage(language);

    translations.clear();

    for (const auto& fname : fnames) {
        if(fname.empty()) {
            continue;
        }

        std::ifstream file(fname.c_str());
        if (!file.is_open()) {
            continue;
        }

        std::string entry;
        auto hint = translations.begin();
        while (std::getline(file, entry)) {

            if (entry.empty() || entry.front() == '#' || entry.front() == '!') {
                continue;
            }

            std::string key, value;

            std::istringstream line(entry);

            if(std::getline(line, key, ';') && translations.find(key) == translations.end() && std::getline(line, value)) {
                size_t pos = 0;
                while((pos = value.find("\\n", pos)) != std::string::npos) {
                     value.replace(pos, 2, "\n");
                     pos++;
                }
                hint = translations.emplace_hint(hint, std::move(key), std::move(value));
            }
        }
    }
}

Glib::ustring MultiLangMgr::getStr (const std::string& key) const
{
    const auto iterator = translations.find(key);

    if (iterator != translations.end()) {
        return iterator->second;
    }

    return key;
}

bool MultiLangMgr::isOSLanguageDetectSupported ()
{
#if defined (WIN32) || defined (__linux__) || defined (__APPLE__)
    return true;
#else
    return false;
#endif
}

Glib::ustring MultiLangMgr::getOSUserLanguage ()
{
    Glib::ustring langName ("default");

#if defined (WIN32)

    const LCID localeID = GetUserDefaultLCID ();
    TCHAR localeName[18];

    const int langLen = GetLocaleInfo (localeID, LOCALE_SISO639LANGNAME, localeName, 9);
    if (langLen <= 0) {
        return langName;
    }

    localeName[langLen - 1] = '-';

    const int countryLen = GetLocaleInfo (localeID, LOCALE_SISO3166CTRYNAME, &localeName[langLen], 9);
    if (countryLen <= 0) {
        return langName;
    }

    langName = localeToLang (localeName);

#elif defined (__linux__)

    // Query the current locale and force decimal point to dot.
    const char *locale = getenv("LANG");
    if (locale || (locale = setlocale (LC_CTYPE, ""))) {
        langName = localeToLang (locale);
    }

    setlocale (LC_NUMERIC, "C");

#elif defined (__APPLE__)

    // "LANG" environment variable is not defined. Retrieving it from CoreFundation API
    // Get used Mac string encoding
    CFStringEncoding strEncoding = CFStringGetSystemEncoding();
    // Get user locale data
    CFLocaleRef cfLocale = CFLocaleCopyCurrent();
    // Get locale language code
    CFStringRef langCodeStr = (CFStringRef)CFLocaleGetValue(cfLocale, kCFLocaleLanguageCode);
    Glib::ustring langCode("");

    if (langCodeStr != NULL) {
        const auto langCodeStrLength = CFStringGetLength(langCodeStr) + 1;
        char langCodeBuffer[langCodeStrLength];
        CFStringGetCString(langCodeStr, langCodeBuffer, langCodeStrLength, strEncoding);
        langCode = Glib::ustring(langCodeBuffer);
    }

    // Get locale country code
    CFStringRef countryCodeStr = (CFStringRef)CFLocaleGetValue(cfLocale, kCFLocaleCountryCode);
    Glib::ustring countryCode("");

    if (countryCodeStr != NULL) {
        const auto countryCodeStrLength = CFStringGetLength(countryCodeStr) + 1;
        char countryCodeBuffer[countryCodeStrLength];
        CFStringGetCString(countryCodeStr, countryCodeBuffer, countryCodeStrLength, strEncoding);
        countryCode = Glib::ustring(countryCodeBuffer);
    }

    // Concatenate locale data
    Glib::ustring locale = langCode + "_" + countryCode;

    // Release user locale data
    CFRelease(cfLocale);
    CFRelease(langCodeStr);
    CFRelease(countryCodeStr);

    // Set locale environment data
    locale = locale + ".UTF-8"; // According to Apple documentation, UTF-8 is a built-in encoding on all platforms on which macOS runs
    g_setenv("LANG", locale.c_str(), true);
    setlocale(LC_ALL, locale.c_str());
    setlocale (LC_NUMERIC, "C"); // Force decimal point to dot.

    langName = localeToLang(locale);
    
#endif

    return langName;
}

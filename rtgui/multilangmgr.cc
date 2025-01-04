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
#include <iostream>
#include <utility>
#ifdef _WIN32
#include <windows.h>
#include <winnls.h>
#endif
#ifdef __APPLE__
#include <CoreFoundation/CoreFoundation.h>
#endif

#include "rtengine/settings.h"

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
        emplace (key ("de", ""  ), "Deutsch");
#ifdef __APPLE__
        emplace (key ("en", "UK"), "English (UK)");
#else
        emplace (key ("en", "GB"), "English (UK)");
#endif
        emplace (key ("en", "US"), "English (US)");
        emplace (key ("es", ""  ), "Espanol (Latin America)");
        emplace (key ("es", "ES"), "Espanol (Castellano)");
        emplace (key ("eu", "ES"), "Euskara");
        emplace (key ("fr", ""  ), "Francais");
        emplace (key ("el", "GR"), "Greek");
        emplace (key ("he", "IL"), "Hebrew");
        emplace (key ("it", ""  ), "Italiano");
        emplace (key ("ja", "JP"), "Japanese");
        emplace (key ("lv", ""  ), "Latvian");
        emplace (key ("hu", ""  ), "Magyar");
        emplace (key ("nl", ""  ), "Nederlands");
        emplace (key ("nn", "NO"), "Norsk BM");
        emplace (key ("nb", "NO"), "Norsk BM");
        emplace (key ("pl", ""  ), "Polish");
        emplace (key ("pt", ""  ), "Portugues (Brasil)");
        emplace (key ("ru", ""  ), "Russian");
        emplace (key ("sr", "RS"), "Serbian (Cyrilic Characters)");
        emplace (key ("sk", ""  ), "Slovak");
        emplace (key ("fi", ""  ), "Suomi");
        emplace (key ("sv", "SE"), "Swedish");
        emplace (key ("tr", ""  ), "Turkish");
        emplace (key ("uk", ""  ), "Ukrainian");
        emplace (key ("zh", "CN"), "Chinese (Simplified)");
        emplace (key ("zh", "SG"), "Chinese (Traditional)");
    }

    Glib::ustring operator() (const Glib::ustring& locale) const
    {
        Glib::ustring major, minor;

        // TODO: Support 3 character language code when needed.
        if (locale.length () >= 2) {
            major = locale.substr (0, 2).lowercase ();
        }

        if (locale.length () >= 5) {
            const Glib::ustring::size_type length =
                locale.length() > 5 && g_unichar_isalnum(locale[5]) ? 3 : 2;
            minor = locale.substr (3, length).uppercase ();
        }

        // Look for matching language and country.
        auto iterator = find (key (major, minor));

        if (iterator != end ()) {
            return iterator->second;
        }

        // Look for matching language only.
        iterator = find (key (major, ""));

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

TranslationMetadata::TranslationMetadata(std::map<std::string, std::string> &&metadata) :
    metadata(std::move(metadata))
{
}

std::string TranslationMetadata::get(const std::string &key, const std::string &default_value) const
{
    const auto found_entry = metadata.find(key);
    if (found_entry == metadata.end()) {
        return default_value;
    }
    return found_entry->second;
}

std::string TranslationMetadata::getLanguageName(const std::string &default_name) const
{
    return get("LANGUAGE_DISPLAY_NAME", default_name);
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

const TranslationMetadata *MultiLangMgr::getMetadata(const Glib::ustring &fname) const
{
    static const char comment_symbol = '#';
    static const char *space_chars = " \t";
    static const char var_symbol = '@';
    static const char key_value_separator = '=';

    // Look for the metadata in the cache.
    const auto &found_metadata = lang_files_metadata.find(fname);
    if (found_metadata != lang_files_metadata.end()) {
        return &found_metadata->second;
    }

    std::ifstream file(fname.c_str());
    if (!file.is_open()) {
        if (rtengine::settings->verbose) {
            std::cerr << "Unable to open language file " << fname << " to get metadata." << std::endl;
        }
        return nullptr;
    }

    if (rtengine::settings->verbose) {
        std::cout << "Reading metadata from language file " << fname << std::endl;
    }
    std::map<std::string, std::string> raw_metadata;
    const auto read_key_value = [&raw_metadata](const std::string &meta_line) {
        // One metadata key-value pair per line. The format is as follows:
        // #001 @KEY=VALUE
        // The line must begin with the comment symbol (#). After the first
        // sequence of whitespace characters, the metadata variable symbol (@)
        // must appear. It is followed immediately with the key name. The end of
        // the key name is marked with the equal sign (=). All remaining
        // characters until the end of the line make up the metadata value.
        if (meta_line.empty() || meta_line.front() != comment_symbol) {
            return;
        }
        const auto first_space = meta_line.find_first_of(space_chars, 1);
        if (first_space == std::string::npos) {
            return;
        }
        const auto definition_start = meta_line.find_first_not_of(space_chars, first_space + 1);
        if (definition_start == std::string::npos || meta_line[definition_start] != var_symbol) {
            return;
        }
        const auto separator_pos = meta_line.find(key_value_separator, definition_start + 1);
        if (separator_pos == std::string::npos) {
            return;
        }
        std::string key = meta_line.substr(definition_start + 1, separator_pos - definition_start - 1);
        std::string value = meta_line.substr(separator_pos + 1);
        if (rtengine::settings->verbose) {
            std::cout << "Found metadata key " << key << " with value " << value << std::endl;
        }
        raw_metadata.emplace(std::move(key), std::move(value));
    };

    // Read lines in order. Metadata only appear in the first section of each
    // file.
    for (
        std::string line;
        std::getline(file, line) && (line.empty() ||
                                        line.front() == comment_symbol ||
                                        line.find_first_not_of(space_chars) == std::string::npos);) {
        read_key_value(line);
    }

    // Add metadata to cache and return.
    lang_files_metadata[fname] = TranslationMetadata(std::move(raw_metadata));
    return &lang_files_metadata[fname];
}

bool MultiLangMgr::isOSLanguageDetectSupported ()
{
#if defined (_WIN32) || defined (__linux__) || defined (__APPLE__)
    return true;
#else
    return false;
#endif
}

Glib::ustring MultiLangMgr::getOSUserLanguage ()
{
    Glib::ustring langName ("default");

#if defined (_WIN32)

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

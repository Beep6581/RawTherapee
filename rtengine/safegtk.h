#ifndef SAFE_GTK_H_INCLUDED
#define SAFE_GTK_H_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>
#include <giomm.h>

Glib::ustring safe_filename_to_utf8 (const std::string& src);
Glib::ustring safe_locale_to_utf8 (const std::string& src);
std::string safe_locale_from_utf8 (const Glib::ustring& utf8_str);

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode);
bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test);

#endif

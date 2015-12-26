#ifndef SAFE_GTK_H_INCLUDED
#define SAFE_GTK_H_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>
#include <giomm.h>

bool safe_spawn_command_line_async (const Glib::ustring& cmd_utf8);
bool safe_spawn_command_line_sync (const Glib::ustring& cmd_utf8);

Glib::ustring safe_filename_to_utf8 (const std::string& src);
Glib::ustring safe_locale_to_utf8 (const std::string& src); // from rtengine
std::string safe_locale_from_utf8 (const Glib::ustring& utf8_str);
std::string safe_filename_from_utf8 (const Glib::ustring& utf8_str);

FILE * safe_g_fopen_WriteBinLock(const Glib::ustring& fname);
int safe_open_ReadOnly(const char *fname);

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode);
bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test);

#endif

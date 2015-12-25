#ifndef SAFE_GTK_H_INCLUDED
#define SAFE_GTK_H_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>
#include <giomm.h>

Glib::RefPtr<Gio::FileInfo> safe_query_file_info (Glib::RefPtr<Gio::File> &file);
void safe_build_file_list (Glib::RefPtr<Gio::File> &dir, std::vector<Glib::ustring> &names, const Glib::ustring &directory = "", const std::vector<Glib::ustring> *extensions = NULL);

bool safe_spawn_command_line_async (const Glib::ustring& cmd_utf8);
bool safe_spawn_command_line_sync (const Glib::ustring& cmd_utf8);

Glib::ustring safe_filename_to_utf8 (const std::string& src);
Glib::ustring safe_locale_to_utf8 (const std::string& src); // from rtengine
std::string safe_locale_from_utf8 (const Glib::ustring& utf8_str);
std::string safe_filename_from_utf8 (const Glib::ustring& utf8_str);

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode);
FILE * safe_g_fopen_WriteBinLock(const Glib::ustring& fname);
int safe_open_ReadOnly(const char *fname);

bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test);
int safe_g_remove(const Glib::ustring& filename);
int safe_g_rename(const Glib::ustring& oldFilename, const Glib::ustring& newFilename);
int safe_g_mkdir_with_parents(const Glib::ustring& dirName, int mode);

#endif

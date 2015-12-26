#ifndef SAFE_GTK_H_INCLUDED
#define SAFE_GTK_H_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>
#include <giomm.h>

FILE * safe_g_fopen(const Glib::ustring& src, const gchar *mode);
bool safe_file_test (const Glib::ustring& filename, Glib::FileTest test);

#endif

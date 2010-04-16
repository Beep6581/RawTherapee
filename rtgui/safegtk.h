#ifndef SAFE_GTK_H_INCLUDED
#define SAFE_GTK_H_INCLUDED

#include <gtkmm.h>
#include <glibmm.h>

Glib::RefPtr<Gdk::Pixbuf> safe_create_from_file(const std::string& filename);

class FileMTimeInfo {

    public:
        Glib::ustring fname;
        Glib::TimeVal mtime;
        
        FileMTimeInfo (Glib::ustring name, Glib::TimeVal mtime) : fname(name), mtime(mtime) {}
        bool operator<(const FileMTimeInfo& other) const { return mtime<other.mtime; }
};

Glib::RefPtr<Gio::FileInfo> safe_query_file_info (Glib::RefPtr<Gio::File> &file);
void safe_build_file_list (Glib::RefPtr<Gio::File> &dir, std::vector<FileMTimeInfo> &flist);
void safe_build_file_list (Glib::RefPtr<Gio::File> &dir, std::vector<Glib::ustring> &names, const Glib::ustring &directory = "");
void safe_build_subdir_list (Glib::RefPtr<Gio::File> &dir, std::vector<Glib::ustring> &subDirs, bool add_hidden);

bool safe_spawn_command_line_async (const Glib::ustring& cmd_utf8);

Glib::ustring safe_locale_to_utf8 (const std::string& src); // from rtengine
std::string safe_locale_from_utf8 (const Glib::ustring& utf8_str);



#endif

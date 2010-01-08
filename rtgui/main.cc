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
// generated 2004/6/3 19:15:32 CEST by gabor@darkstar.(none)
// using glademm V2.5.0
//
// newer (non customized) versions of this file go to raw.cc_new

// This file is for your program, I won't touch it again!

//#include <config.h>
#include <gtkmm.h>
#include <giomm.h>
#include <iostream>
#include <rtwindow.h>
#include <string.h>
#include <stdlib.h>
#include <options.h>
extern Options options;

//#ifdef WIN32
//#include <windows.h> // included for WinMain
//#endif

Glib::ustring argv0;
Glib::ustring argv1;

int main(int argc, char **argv)
{  

   std::string argv0_, argv1_;
   
#ifndef WIN32
   argv0_ = argv[0];
#else
   char exname[512];
   GetModuleFileName (NULL, exname, 512);
   argv0_ = exname;
#endif
   int i;
   for (i=argv0_.size()-1; (argv0_[i]!='/' && argv0_[i]!='\\') && i>0; i--);
   if (argv0_[i]=='/' || argv0_[i]=='\\')
     argv0_ = argv0_.substr(0,i);

    if (argc>1)
        argv1_ = argv[1];
    else
        argv1_ = "";

   argv0 = Glib::locale_to_utf8 (argv0_);
   argv1 = Glib::locale_to_utf8 (argv1_);

   Glib::thread_init();
   gdk_threads_init();
   Gio::init ();

   Options::load ();

//   Gtk::RC::add_default_file (argv0+"/themes/"+options.theme);
   std::vector<std::string> rcfiles;
   rcfiles.push_back (argv0+"/themes/"+options.theme);
   Gtk::RC::set_default_files (rcfiles);

   Gtk::Main m(&argc, &argv);
//   MainWindow *MainWindow = new class MainWindow();
   RTWindow *rtWindow = new class RTWindow();
   gdk_threads_enter ();
   m.run(*rtWindow);
   gdk_threads_leave ();
   delete rtWindow;
   return 0;
}





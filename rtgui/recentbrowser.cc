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
#include "recentbrowser.h"
#include "multilangmgr.h"
#include "options.h"

using namespace rtengine;

RecentBrowser::RecentBrowser () : listener (NULL) {

    recentDirs = Gtk::manage (new MyComboBoxText ());

    Gtk::Frame* frame = Gtk::manage (new Gtk::Frame (M("MAIN_FRAME_RECENT")));
    frame->add (*recentDirs);
    for(size_t i=0;i<options.recentFolders.size();i++) {
        recentDirs->append_text (options.recentFolders[i]);
    }

    pack_start (*frame, Gtk::PACK_SHRINK, 4);

    conn = recentDirs->signal_changed().connect(sigc::mem_fun(*this, &RecentBrowser::selectionChanged));

    show_all ();
}

void RecentBrowser::selectionChanged () {

    Glib::ustring sel = recentDirs->get_active_text ();
    if (sel!="" && listener)
        listener->selectDir (sel);
}

void RecentBrowser::dirSelected (const Glib::ustring& dirname, const Glib::ustring& openfile) {
    
    size_t numFolders = options.recentFolders.size();
    if(numFolders>0) { // search entry and move to top if it exists
        size_t i;
        for(i=0;i<numFolders;i++) {
            if(options.recentFolders[i] == dirname) {
                break;
            }
        }
        if(i>0) {
            if(i<numFolders) {
               options.recentFolders.erase(options.recentFolders.begin()+i);
            }
            options.recentFolders.insert(options.recentFolders.begin(),dirname);
        }
    } else {
        options.recentFolders.insert(options.recentFolders.begin(),dirname);
    }

    conn.block (true);

    recentDirs->remove_text (dirname);
    recentDirs->prepend_text (dirname);
    recentDirs->set_active_text (dirname);

    conn.block (false);
}


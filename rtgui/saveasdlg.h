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
#ifndef _SAVEASDLG_
#define _SAVEASDLG_

#include <gtkmm.h>
#include <adjuster.h>
#include <saveformatpanel.h>
#include <options.h>

class SaveAsDialog : public Gtk::Dialog, public FormatChangeListener {

  protected:
    Gtk::FileChooserWidget* fchooser;
    SaveFormatPanel* formatOpts;
    Glib::ustring fname;
    Gtk::FileFilter filter_jpg;
    Gtk::FileFilter filter_tif;
    Gtk::FileFilter filter_png;
    Gtk::RadioButton* immediately;
    Gtk::RadioButton* putToQueueHead;
    Gtk::RadioButton* putToQueueTail;

  public:
    SaveAsDialog (Glib::ustring initialDir);
    
    Glib::ustring   getFileName      ();
    Glib::ustring   getDirectory     ();
    SaveFormat      getFormat        ();
    bool            getImmediately   ();
    bool            getToHeadOfQueue ();
    bool            getToTailOfQueue ();    

    void  setInitialFileName (Glib::ustring iname);

    void okPressed ();
    void cancelPressed ();
    void formatChanged (Glib::ustring f);
};


#endif

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
#ifndef _ICMPANEL_
#define _ICMPANEL_

#include <gtkmm.h>
#include <toolpanel.h>

class ICMPanelListener {

    public:
        virtual void saveInputICCReference (Glib::ustring fname) {}
};

class ICMPanel : public Gtk::VBox, public ToolPanel {

    private: 
        Gtk::RadioButton*  iembedded;
        Gtk::RadioButton*  icamera;
        Gtk::RadioButton*  ifromfile;
        Gtk::CheckButton*  igamma;
        Gtk::ComboBoxText* wnames;
        Gtk::ComboBoxText* onames;
        Gtk::RadioButton*  ofromdir;
        Gtk::RadioButton*  ofromfile;
        Gtk::RadioButton*  iunchanged;
        Gtk::FileChooserButton* ipDialog;
        Gtk::FileChooserButton* opDialog;
        Gtk::RadioButton::Group opts;
        Gtk::Button*        saveRef;
        sigc::connection   ipc;
        Glib::ustring      oldip;
        ICMPanelListener*  icmplistener;
        
    public:
        ICMPanel ();
        
        void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
        void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
        void setBatchMode   (bool batchMode);
        
        void wpChanged ();
        void opChanged ();
        void ipChanged ();

        void ipSelectionChanged ();

        void setRaw (bool raw);
        void saveReferencePressed ();

        void setICMPanelListener (ICMPanelListener* ipl) { icmplistener = ipl; }
};

#endif

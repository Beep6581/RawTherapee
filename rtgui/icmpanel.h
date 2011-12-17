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
#include "adjuster.h"
#include "guiutils.h"

#include "toolpanel.h"
#include "../rtengine/imagedata.h"

class ICMPanelListener {

    public:
        virtual void saveInputICCReference (Glib::ustring fname) {}
};

class ICMPanel : public Gtk::VBox, public AdjusterListener, public FoldableToolPanel {

	protected:
	Adjuster* gampos;
	Adjuster* slpos;
	bool lastgamfree;
    sigc::connection  gamcsconn; 
	//bool freegamma;
    private:
        Gtk::CheckButton*  freegamma;
    	Gtk::RadioButton*  inone;
		
        Gtk::RadioButton*  iembedded;
        Gtk::RadioButton*  icamera;
        Gtk::RadioButton*  icameraICC;
        Gtk::RadioButton*  ifromfile;
        Gtk::CheckButton*  ckbBlendCMSMatrix;
        MyComboBoxText*    wnames;
        MyComboBoxText*    wgamma;
		
        MyComboBoxText*    onames;
        Gtk::RadioButton*  ofromdir;
        Gtk::RadioButton*  ofromfile;
        Gtk::RadioButton*  iunchanged;
        MyFileChooserButton* ipDialog;
        Gtk::RadioButton::Group opts;
        Gtk::Button*        saveRef;
        sigc::connection   ipc;
        Glib::ustring      oldip;
        ICMPanelListener*  icmplistener;
        
        static Glib::ustring lastICCWorkDir;
        bool enableLastICCWorkDirChange;

    public:
        ICMPanel ();
        
        void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited=NULL); 
        void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited=NULL);
        void setBatchMode   (bool batchMode);
		void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited=NULL);
        void adjusterChanged (Adjuster* a, double newval);
		void setAdjusterBehavior (bool gammaadd, bool slopeadd);
 
        void wpChanged ();
        void opChanged ();
        void ipChanged ();
        void gpChanged ();
		void GamChanged ();
        void ipSelectionChanged ();
        void iccTogglesChanged();

        void setRawMeta (bool raw, const rtengine::ImageData* pMeta);
        void saveReferencePressed ();

        void setICMPanelListener (ICMPanelListener* ipl) { icmplistener = ipl; }
};

#endif

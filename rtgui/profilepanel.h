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
#ifndef _PROFILEPANEL_
#define _PROFILEPANEL_

#include <gtkmm.h>
#include <vector>
#include "../rtengine/rtengine.h"
#include "pparamschangelistener.h"
#include "profilechangelistener.h"
#include "partialpastedlg.h"
#include "guiutils.h"

class ProfilePanel : public Gtk::VBox, public PParamsChangeListener {

  protected:

    static PartialPasteDlg* partialProfileDlg;
    Gtk::Button* save;
    Gtk::Button* load;
    Gtk::Button* copy;
    Gtk::Button* paste;
    MyComboBoxText* profiles;
    std::vector<Glib::ustring> pparams;
    rtengine::procparams::PartialProfile* custom;
    rtengine::procparams::PartialProfile* lastsaved;
    rtengine::procparams::PartialProfile* lastphoto;
    Glib::ustring old;
    ProfileChangeListener* tpc;
    bool dontupdate;
    sigc::connection changeconn;
    Gtk::FileChooserDialog* savedialog;

    void changeTo (rtengine::procparams::PartialProfile* newpp, Glib::ustring profname);
    void refreshProfileList ();

  public:

    ProfilePanel (bool readOnly=false);
    virtual ~ProfilePanel ();

    void setProfileChangeListener (ProfileChangeListener* ppl) { tpc = ppl; }

    static void init ();
    static void cleanup ();

    void initProfile (const Glib::ustring& profname, rtengine::procparams::ProcParams* lastSaved, rtengine::procparams::ProcParams* lastPhoto);

    // PParamsChangeListener interface
    void procParamsChanged (rtengine::procparams::ProcParams* params, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited=NULL);

    // gui callbacks
    void save_clicked (GdkEventButton* event);
    void load_clicked (GdkEventButton* event);
    void copy_clicked (GdkEventButton* event);
    void paste_clicked (GdkEventButton* event);
    void selection_changed ();
};

#endif

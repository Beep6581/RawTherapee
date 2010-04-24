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
#include <profilepanel.h>
#include <options.h>
#include <profilestore.h>
#include <clipboard.h>
#include <multilangmgr.h>

using namespace rtengine;
using namespace rtengine::procparams;

extern Glib::ustring argv0;

ProfilePanel::ProfilePanel () {

    tpc = NULL;
  
    profiles = Gtk::manage (new Gtk::ComboBoxText ());
    Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
    hbox->show ();
//    pack_start (*profiles, Gtk::PACK_SHRINK, 4);

    pack_start (*hbox, Gtk::PACK_SHRINK, 4);
    
    save = Gtk::manage (new Gtk::Button ());
    save->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
    load = Gtk::manage (new Gtk::Button ());
    load->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));
    copy = Gtk::manage (new Gtk::Button ());
    copy->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-copy"), Gtk::ICON_SIZE_BUTTON)));
    paste = Gtk::manage (new Gtk::Button ());
    paste->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-paste"), Gtk::ICON_SIZE_BUTTON)));

    hbox->pack_start (*profiles);
    hbox->pack_start (*load, Gtk::PACK_SHRINK, 1);
    hbox->pack_start (*save, Gtk::PACK_SHRINK, 1);
    hbox->pack_start (*copy, Gtk::PACK_SHRINK, 1);
    hbox->pack_start (*paste, Gtk::PACK_SHRINK, 1);

    load->signal_clicked().connect( sigc::mem_fun(*this, &ProfilePanel::load_clicked) );
    save->signal_clicked().connect( sigc::mem_fun(*this, &ProfilePanel::save_clicked) );
    copy->signal_clicked().connect( sigc::mem_fun(*this, &ProfilePanel::copy_clicked) );
    paste->signal_clicked().connect( sigc::mem_fun(*this, &ProfilePanel::paste_clicked) );

    custom = NULL;
    lastphoto = NULL;
    lastsaved = NULL;
    dontupdate = false;

    refreshProfileList ();

    profiles->set_active (0);
    old = profiles->get_active_text();
    changeconn = profiles->signal_changed().connect( sigc::mem_fun(*this, &ProfilePanel::selection_changed) );

    save->set_tooltip_text (M("PROFILEPANEL_TOOLTIPSAVE"));
    load->set_tooltip_text (M("PROFILEPANEL_TOOLTIPLOAD"));
    copy->set_tooltip_text (M("PROFILEPANEL_TOOLTIPCOPY"));
    paste->set_tooltip_text (M("PROFILEPANEL_TOOLTIPPASTE"));
    
    show_all_children ();
}

ProfilePanel::~ProfilePanel () {

    delete custom;
    delete lastsaved;
    delete lastphoto;
}

void ProfilePanel::refreshProfileList () {

    Glib::ustring oldsel = profiles->get_active_text ();
    changeconn.block (true);

    // clear items
    profiles->clear_items ();
    pparams.clear ();

    // re-parse profile directories (deletes old ones)
    profileStore.parseProfiles ();
    pparams = profileStore.getProfileNames ();
    for (int i=0; i<pparams.size(); i++)
        profiles->append_text (pparams[i]);

    if (custom)
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
    if (lastsaved)
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")");
    if (lastphoto)
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")");

    profiles->set_active_text (oldsel);
    changeconn.block (false);
}

void ProfilePanel::save_clicked () {

    Gtk::FileChooserDialog dialog(M("PROFILEPANEL_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
    if (options.multiUser)
    dialog.set_current_folder (Options::rtdir + "/" + options.profilePath);
    else
    dialog.set_current_folder (argv0 + "/" + options.profilePath);

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("PROFILEPANEL_FILEDLGFILTERPP"));
    filter_pp.add_pattern("*.pp2");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("PROFILEPANEL_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

//    dialog.set_do_overwrite_confirmation (true);

    savedialog = &dialog;

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {

        std::string fname = dialog.get_filename();

        bool hasext = true;
        int dotpos = fname.find_last_of ('.');
        if (dotpos==Glib::ustring::npos)
            hasext = false;
        int dirpos1 = fname.find_last_of ('/');
        if (dirpos1!=Glib::ustring::npos && dirpos1>dotpos)
            hasext = false;
        int dirpos2 = fname.find_last_of ('\\');
        if (dirpos2!=Glib::ustring::npos && dirpos2>dotpos)
            hasext = false;

        if (!hasext) 
            fname = fname + ".pp2";

        if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS)) {
          Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
          Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
          int response = msgd.run ();
          if (response==Gtk::RESPONSE_NO)
            return;
        }
        
        ProcParams* toSave = NULL;
        if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")") 
            toSave = custom;
        else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")") 
            toSave = lastsaved; 
        else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")") 
            toSave = lastphoto; 
        else
            toSave = profileStore.getProfile (profiles->get_active_text());
            
        if (toSave) {
            toSave->save (fname);
            refreshProfileList ();
        }
    }
}

void ProfilePanel::copy_clicked () {

    ProcParams* toSave = NULL;
    if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")") 
        toSave = custom;
    else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")") 
        toSave = lastsaved; 
    else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")") 
        toSave = lastphoto; 
    else
        toSave = profileStore.getProfile (profiles->get_active_text());
        
    if (toSave)
        clipboard.setProcParams (*toSave);
}

void ProfilePanel::load_clicked () {

    Gtk::FileChooserDialog dialog(M("PROFILEPANEL_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);
    if (options.multiUser)
    dialog.set_current_folder (Options::rtdir + "/" + options.profilePath);
    else
    dialog.set_current_folder (argv0 + "/" + options.profilePath);

    //Add response buttons the the dialog:
    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    //Add filters, so that only certain file types can be selected:

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("PROFILEPANEL_FILEDLGFILTERPP"));
    filter_pp.add_pattern("*.pp2");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("PROFILEPANEL_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {
        if (!custom) {
            custom = new ProcParams ();
            profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
        }
        custom->load (dialog.get_filename());
        profiles->set_active_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
        old = profiles->get_active_text();
        changeTo (custom, M("PROFILEPANEL_PFILE"));
    }
}

void ProfilePanel::paste_clicked () {

    if (!clipboard.hasProcParams())
        return;

    if (!custom) {
        custom = new ProcParams ();
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
    }
    *custom = clipboard.getProcParams ();
    profiles->set_active_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
    old = profiles->get_active_text();
    changeTo (custom, M("HISTORY_FROMCLIPBOARD"));
}

void ProfilePanel::changeTo (ProcParams* newpp, Glib::ustring profname) {

    if (!newpp)
        return;

  // Keep transformation parameters while changing the profile

/*  int cropx = working->crop_x;
  int cropy = working->crop_y;
  int cropw = working->crop_w;
  int croph = working->crop_h;
  bool crope = working->crop_enabled;
  int rotcor = working->rotate_coarse;
  double rotfine  = working->rotate_fine;
  double lenscorr = working->lens_distortion;
  bool hflip      = working->horizontal_flip;
  bool vflip      = working->vertical_flip;

  working->copy (newpp);

  working->crop_x = cropx;
  working->crop_y = cropy;
  working->crop_w = cropw;
  working->crop_h = croph;
  working->crop_enabled    = crope;
  working->rotate_coarse   = rotcor;
  working->rotate_fine     = rotfine;
  working->lens_distortion = lenscorr;
  working->horizontal_flip = hflip;
  working->vertical_flip   = vflip;
*/
    if (tpc)
        tpc->profileChange (newpp, EvProfileChanged, profname);  
}

void ProfilePanel::selection_changed () {

    if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")") {
        if (!dontupdate) 
            changeTo (custom, Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");          
    }             
    else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")") 
            changeTo (lastsaved, Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")");                       
    else if (profiles->get_active_text() == Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")") 
            changeTo (lastphoto, Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")");                       
    else {
        ProcParams* s = profileStore.getProfile (profiles->get_active_text());
        if (s)
            changeTo (s, profiles->get_active_text());
    }
    old = profiles->get_active_text ();
    dontupdate = false;
}

void ProfilePanel::procParamsChanged (rtengine::procparams::ProcParams* p, rtengine::ProcEvent ev, Glib::ustring descr, ParamsEdited* paramsEdited) {

    // to prevent recursion filter out the events caused by the profilepanel
    if (ev==EvProfileChanged || ev==EvPhotoLoaded)
        return;

    if (profiles->get_active_text() != Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")") {
        dontupdate = true;
        if (!custom) {
            custom = new ProcParams ();
            profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
        }
        *custom = *p;
        profiles->set_active_text (Glib::ustring("(") + M("PROFILEPANEL_PCUSTOM") + ")");
        old = profiles->get_active_text();
    }
    else
        *custom = *p;
}

void ProfilePanel::initProfile (const Glib::ustring& profname, ProcParams* lastSaved, ProcParams* lastPhoto) {

    changeconn.block (true);

    profiles->clear_items ();
    pparams.clear ();

    pparams = profileStore.getProfileNames ();
    for (int i=0; i<pparams.size(); i++)
        profiles->append_text (pparams[i]);

    delete custom;
    custom = NULL;
    delete lastsaved;
    lastsaved = lastSaved;
    delete lastphoto;
    lastphoto = lastPhoto;
    
    Glib::ustring defline = profname;
    ProcParams* defprofile = profileStore.getProfile (profname);

    if (lastphoto) 
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PLASTPHOTO") + ")");  

    if (lastsaved) {
        defline = Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")";
        defprofile = lastsaved;
        profiles->append_text (Glib::ustring("(") + M("PROFILEPANEL_PLASTSAVED") + ")");  
    }

    if (tpc) {
        if (lastsaved)
            tpc->setDefaults (lastsaved);
        else
            tpc->setDefaults (profileStore.getProfile (profname));
    }
    if (defprofile) {
        old = defline;
        profiles->set_active_text (defline);
        changeconn.block (false);
        if (tpc)
            tpc->profileChange (defprofile, EvPhotoLoaded, defline);  
    }
    else {
        // select first valid profile
        old = "";
        profiles->set_active (0);
        ProcParams* s = profileStore.getProfile (profiles->get_active_text());
        if (!s)
            s = new ProcParams ();
        changeconn.block (false);
        if (tpc)
            tpc->profileChange (s, EvPhotoLoaded, profiles->get_active_text());  
    }
}



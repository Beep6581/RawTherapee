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
#include "filterpanel.h"
#include "multilangmgr.h"
#include "../rtengine/rtengine.h"
#include "rtimage.h"

using namespace rtengine;

FilterPanel::FilterPanel () : listener (NULL)
{

    set_border_width (4);

    enabled = Gtk::manage (new Gtk::CheckButton (M("EXIFFILTER_METADATAFILTER")));
    pack_start (*enabled, Gtk::PACK_SHRINK, 2);
    pack_start (*Gtk::manage(new Gtk::HSeparator ()), Gtk::PACK_SHRINK, 2);

    enaFNumber = Gtk::manage (new Gtk::CheckButton (M("EXIFFILTER_APERTURE") + ":"));
    Gtk::VBox* fnvb = Gtk::manage(new Gtk::VBox ());
    Gtk::HBox* fnhb = Gtk::manage(new Gtk::HBox ());
    fnvb->pack_start (*enaFNumber, Gtk::PACK_SHRINK, 0);
    fnumberFrom = Gtk::manage(new Gtk::Entry ());
    fnumberTo = Gtk::manage(new Gtk::Entry ());
    fnhb->pack_start (*fnumberFrom, true, true, 2);
    fnhb->pack_start (*Gtk::manage(new Gtk::Label(" - ")), false, false, 4);
    fnhb->pack_start (*fnumberTo, true, true, 2);
    fnvb->pack_start (*fnhb, Gtk::PACK_SHRINK, 0);
    pack_start (*fnvb, Gtk::PACK_SHRINK, 4);

    enaShutter = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_SHUTTER") + ":"));
    Gtk::VBox* svb = Gtk::manage(new Gtk::VBox ());
    Gtk::HBox* shb = Gtk::manage(new Gtk::HBox ());
    svb->pack_start (*enaShutter, Gtk::PACK_SHRINK, 0);
    shutterFrom = Gtk::manage(new Gtk::Entry ());
    shutterTo = Gtk::manage(new Gtk::Entry ());
    shb->pack_start (*shutterFrom, true, true, 2);
    shb->pack_start (*Gtk::manage(new Gtk::Label(" - ")), false, false, 4);
    shb->pack_start (*shutterTo, true, true, 2);
    svb->pack_start (*shb, Gtk::PACK_SHRINK, 0);
    pack_start (*svb, Gtk::PACK_SHRINK, 4);

    enaISO = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_ISO") + ":"));
    Gtk::VBox* ivb = Gtk::manage(new Gtk::VBox ());
    Gtk::HBox* ihb = Gtk::manage(new Gtk::HBox ());
    ivb->pack_start (*enaISO, Gtk::PACK_SHRINK, 0);
    isoFrom = Gtk::manage(new Gtk::Entry ());
    isoTo = Gtk::manage(new Gtk::Entry ());
    ihb->pack_start (*isoFrom, true, true, 2);
    ihb->pack_start (*Gtk::manage(new Gtk::Label(" - ")), false, false, 4);
    ihb->pack_start (*isoTo, true, true, 2);
    ivb->pack_start (*ihb, Gtk::PACK_SHRINK, 0);
    pack_start (*ivb, Gtk::PACK_SHRINK, 4);

    enaFocalLen = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_FOCALLEN") + ":"));
    Gtk::VBox* fvb = Gtk::manage(new Gtk::VBox ());
    Gtk::HBox* fhb = Gtk::manage(new Gtk::HBox ());
    fvb->pack_start (*enaFocalLen, Gtk::PACK_SHRINK, 0);
    focalFrom = Gtk::manage(new Gtk::Entry ());
    focalTo = Gtk::manage(new Gtk::Entry ());
    fhb->pack_start (*focalFrom, true, true, 2);
    fhb->pack_start (*Gtk::manage(new Gtk::Label(" - ")), false, false, 4);
    fhb->pack_start (*focalTo, true, true, 2);
    fvb->pack_start (*fhb, Gtk::PACK_SHRINK, 0);
    pack_start (*fvb, Gtk::PACK_SHRINK, 4);

    enaExpComp = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_EXPOSURECOMPENSATION") + ":"));
    Gtk::VBox* evb = Gtk::manage(new Gtk::VBox ());
    evb->pack_start (*enaExpComp, Gtk::PACK_SHRINK, 0);
    expcomp = Gtk::manage(new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE));
    expcomp->set_headers_visible (false);
    Gtk::ScrolledWindow* sexpcomp = Gtk::manage(new Gtk::ScrolledWindow());
    sexpcomp->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    sexpcomp->set_size_request(-1, 80);
    sexpcomp->add(*expcomp);
    evb->pack_start (*sexpcomp, Gtk::PACK_SHRINK, 0);
    pack_start (*evb, Gtk::PACK_SHRINK, 4);

    enaCamera = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_CAMERA") + ":"));
    Gtk::VBox* cvb = Gtk::manage(new Gtk::VBox ());
    cvb->pack_start (*enaCamera, Gtk::PACK_SHRINK, 0);
    camera = Gtk::manage(new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE));
    camera->set_headers_visible (false);
    Gtk::ScrolledWindow* scamera = Gtk::manage(new Gtk::ScrolledWindow());
    scamera->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    scamera->set_size_request(-1, 80);
    scamera->add(*camera);
    cvb->pack_start (*scamera, Gtk::PACK_SHRINK, 0);
    pack_start (*cvb, Gtk::PACK_SHRINK, 4);

    enaLens = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_LENS") + ":"));
    Gtk::VBox* lvb = Gtk::manage(new Gtk::VBox ());
    lvb->pack_start (*enaLens, Gtk::PACK_SHRINK, 0);
    lens = Gtk::manage(new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE));
    lens->set_headers_visible (false);
    Gtk::ScrolledWindow* slens = Gtk::manage(new Gtk::ScrolledWindow());
    slens->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    slens->set_size_request(-1, 80);
    slens->add(*lens);
    lvb->pack_start (*slens, Gtk::PACK_SHRINK, 0);
    pack_start (*lvb, Gtk::PACK_SHRINK, 4);

    enaFiletype = Gtk::manage(new Gtk::CheckButton(M("EXIFFILTER_FILETYPE") + ":"));
    Gtk::VBox* ftvb = Gtk::manage(new Gtk::VBox ());
    ftvb->pack_start (*enaFiletype, Gtk::PACK_SHRINK, 0);
    filetype = Gtk::manage(new Gtk::ListViewText (1, false, Gtk::SELECTION_MULTIPLE));
    filetype->set_headers_visible (false);
    Gtk::ScrolledWindow* sfiletype = Gtk::manage(new Gtk::ScrolledWindow());
    sfiletype->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_ALWAYS);
    sfiletype->set_size_request(-1, 80);
    sfiletype->add(*filetype);
    ftvb->pack_start (*sfiletype, Gtk::PACK_SHRINK, 0);
    pack_start (*ftvb, Gtk::PACK_SHRINK, 4);

    // add panel ending
    Gtk::VBox* vboxpe = Gtk::manage (new Gtk::VBox ());
    Gtk::HSeparator* hseptpe = Gtk::manage (new Gtk::HSeparator ());
    Gtk::Image* peImg = Gtk::manage (new RTImage("PanelEnding.png"));
    vboxpe->pack_start(*hseptpe, Gtk::PACK_SHRINK, 4);
    vboxpe->pack_start(*peImg);
    pack_start(*vboxpe, Gtk::PACK_SHRINK, 0);

    conns = 0;
    sChange[conns++] = fnumberFrom->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = fnumberTo->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = shutterFrom->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = shutterTo->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = isoFrom->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = isoTo->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = focalFrom->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = focalTo->signal_changed().connect (sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = expcomp->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = filetype->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = camera->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = lens->get_selection()->signal_changed().connect(sigc::mem_fun(*this, &FilterPanel::valueChanged));
    sChange[conns++] = enaFNumber->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaShutter->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaFocalLen->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaISO->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaExpComp->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaCamera->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaLens->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enabled->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );
    sChange[conns++] = enaFiletype->signal_toggled().connect( sigc::mem_fun(*this, &FilterPanel::valueChanged) );

    set_size_request (0, -1);

    show_all ();
}

void FilterPanel::setFilter (ExifFilterSettings& defefs, bool updateLists)
{


    for (int i = 0; i < conns; i++) {
        sChange[i].block (true);
    }

//  enaFNumber->set_active (curefs.filterFNumber);
    fnumberFrom->set_text (ImageMetaData::apertureToString (defefs.fnumberFrom));
    curefs.fnumberFrom = defefs.fnumberFrom;
    fnumberTo->set_text (ImageMetaData::apertureToString (defefs.fnumberTo));
    curefs.fnumberTo = defefs.fnumberTo;

//  enaShutter->set_active (curefs.filterShutter);
    shutterFrom->set_text (ImageMetaData::shutterToString (defefs.shutterFrom));
    curefs.shutterFrom = defefs.shutterFrom;
    shutterTo->set_text (ImageMetaData::shutterToString (defefs.shutterTo));
    curefs.shutterTo = defefs.shutterTo;

//  enaISO->set_active (curefs.filterISO);
    isoFrom->set_text (Glib::ustring::format (defefs.isoFrom));
    curefs.isoFrom = defefs.isoFrom;
    isoTo->set_text (Glib::ustring::format (defefs.isoTo));
    curefs.isoTo = defefs.isoTo;

//  enaFocalLen->set_active (curefs.filterFocalLen);
    focalFrom->set_text (Glib::ustring::format (defefs.focalFrom));
    curefs.focalFrom = defefs.focalFrom;
    focalTo->set_text (Glib::ustring::format (defefs.focalTo));
    curefs.focalTo = defefs.focalTo;

//  enaCompExp->set_active (curefs.filterExpComp);
    Glib::RefPtr<Gtk::TreeSelection> eselection = expcomp->get_selection ();

//  enaFiletype->set_active (curefs.filterFiletype);
    Glib::RefPtr<Gtk::TreeSelection> ftselection = filetype->get_selection ();

//  enaCamera->set_active (curefs.filterCamera);
    Glib::RefPtr<Gtk::TreeSelection> cselection = camera->get_selection ();

//  enaLens->set_active (curefs.filterLens);
    Glib::RefPtr<Gtk::TreeSelection> lselection = lens->get_selection ();

    if( updateLists ) {
        expcomp->clear_items();
        curefs.expcomp.clear();

        for (std::set<std::string>::iterator i = defefs.expcomp.begin(); i != defefs.expcomp.end(); i++) {
            expcomp->append (*i);
            curefs.expcomp.insert(*i);
        }

        eselection->select_all();

        lens->clear_items();
        curefs.lenses.clear();

        for (std::set<std::string>::iterator i = defefs.lenses.begin(); i != defefs.lenses.end(); i++) {
            lens->append (*i);
            curefs.lenses.insert(*i);
        }

        lselection->select_all();

        camera->clear_items();
        curefs.cameras.clear();

        for (std::set<std::string>::iterator i = defefs.cameras.begin(); i != defefs.cameras.end(); i++) {
            camera->append(*i);
            curefs.cameras.insert(*i);
        }

        cselection->select_all();

        filetype->clear_items();
        curefs.filetypes.clear();

        for (std::set<std::string>::iterator i = defefs.filetypes.begin(); i != defefs.filetypes.end(); i++) {
            filetype->append(*i);
            curefs.filetypes.insert(*i);
        }

        ftselection->select_all();
    } else {
        for( Gtk::TreeModel::Children::iterator iter = expcomp->get_model()->children().begin(); iter != expcomp->get_model()->children().end(); iter++) {
            Glib::ustring v;
            iter->get_value(0, v);

            if( defefs.expcomp.find( v ) != defefs.expcomp.end() ) {
                eselection->select( iter );
            } else {
                eselection->unselect( iter );
            }
        }

        for( Gtk::TreeModel::Children::iterator iter = lens->get_model()->children().begin(); iter != lens->get_model()->children().end(); iter++) {
            Glib::ustring v;
            iter->get_value(0, v);

            if( defefs.lenses.find( v ) != defefs.lenses.end() ) {
                lselection->select( iter );
            } else {
                lselection->unselect( iter );
            }
        }

        for( Gtk::TreeModel::Children::iterator iter = camera->get_model()->children().begin(); iter != camera->get_model()->children().end(); iter++) {
            Glib::ustring v;
            iter->get_value(0, v);

            if( defefs.cameras.find( v ) != defefs.cameras.end() ) {
                cselection->select(iter);
            } else {
                cselection->unselect(iter);
            }
        }

        for( Gtk::TreeModel::Children::iterator iter = filetype->get_model()->children().begin(); iter != filetype->get_model()->children().end(); iter++) {
            Glib::ustring v;
            iter->get_value(0, v);

            if( defefs.filetypes.find( v ) != defefs.filetypes.end() ) {
                ftselection->select(iter);
            } else {
                ftselection->unselect(iter);
            }
        }
    }

    curefs = defefs;

    for (int i = 0; i < conns; i++) {
        sChange[i].block (false);
    }
}

bool FilterPanel::isEnabled ()
{

    return enabled->get_active () && is_sensitive();
}

ExifFilterSettings FilterPanel::getFilter ()
{

    ExifFilterSettings efs;
    efs.fnumberFrom = atof (fnumberFrom->get_text().c_str());
    efs.fnumberTo   = atof (fnumberTo->get_text().c_str());
    efs.focalFrom   = atof (focalFrom->get_text().c_str());
    efs.focalTo     = atof (focalTo->get_text().c_str());
    efs.isoFrom     = atoi (isoFrom->get_text().c_str());
    efs.isoTo       = atoi (isoTo->get_text().c_str());
    efs.shutterFrom = ImageMetaData::shutterFromString (shutterFrom->get_text());
    efs.shutterTo   = ImageMetaData::shutterFromString (shutterTo->get_text());

    efs.filterFNumber  = enaFNumber->get_active ();
    efs.filterShutter  = enaShutter->get_active ();
    efs.filterFocalLen = enaFocalLen->get_active ();
    efs.filterISO      = enaISO->get_active ();
    efs.filterExpComp  = enaExpComp->get_active ();
    efs.filterCamera   = enaCamera->get_active ();
    efs.filterLens     = enaLens->get_active ();
    efs.filterFiletype = enaFiletype->get_active ();

    std::vector<int> sel = camera->get_selected ();

    for (size_t i = 0; i < sel.size(); i++) {
        efs.cameras.insert (camera->get_text (sel[i]));
    }

    sel = expcomp->get_selected ();

    for (size_t i = 0; i < sel.size(); i++) {
        efs.expcomp.insert (expcomp->get_text (sel[i]));
    }

    sel = lens->get_selected ();

    for (size_t i = 0; i < sel.size(); i++) {
        efs.lenses.insert (lens->get_text (sel[i]));
    }

    sel = filetype->get_selected ();

    for (size_t i = 0; i < sel.size(); i++) {
        efs.filetypes.insert (filetype->get_text (sel[i]));
    }

    return efs;
}

// Called within GTK UI thread
void FilterPanel::valueChanged ()
{

    if (listener) {
        listener->exifFilterChanged ();
    }
}

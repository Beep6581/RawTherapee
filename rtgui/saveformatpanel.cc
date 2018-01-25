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
#include "saveformatpanel.h"
#include "multilangmgr.h"
#include "guiutils.h"

SaveFormatPanel::SaveFormatPanel () : listener (nullptr)
{


    // ---------------------  FILE FORMAT SELECTOR


    Gtk::Grid* hb1 = Gtk::manage (new Gtk::Grid ());
    hb1->set_column_spacing(5);
    hb1->set_row_spacing(5);
    setExpandAlignProperties(hb1, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    Gtk::Label* flab = Gtk::manage (new Gtk::Label (M("SAVEDLG_FILEFORMAT") + ":"));
    setExpandAlignProperties(flab, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    format = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(format, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    format->signal_changed ().connect (sigc::mem_fun (*this, &SaveFormatPanel::formatChanged));

    format->append ("JPEG (8 bit)");
    format->append ("TIFF (8 bit)");
    format->append ("TIFF (16 bit)");
    format->append ("TIFF (32 bit float)");
    format->append ("PNG (8 bit)");
    format->append ("PNG (16 bit)");

    fstr[0] = "jpg";
    fstr[1] = "tif";
    fstr[2] = "tif";
    fstr[3] = "tif";
    fstr[4] = "png";
    fstr[5] = "png";

    hb1->attach (*flab, 0, 0, 1, 1);
    hb1->attach (*format, 1, 0, 1, 1);
    hb1->show_all();

    // ---------------------  JPEG OPTIONS


    jpegOpts = Gtk::manage (new Gtk::Grid ());
    jpegOpts->set_column_spacing(15);
    jpegOpts->set_row_spacing(5);
    setExpandAlignProperties(jpegOpts, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);

    jpegQual = new Adjuster (M("SAVEDLG_JPEGQUAL"), 0, 100, 1, 100);
    setExpandAlignProperties(jpegQual, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    jpegQual->setAdjusterListener (this);

    jpegSubSampLabel = Gtk::manage (new Gtk::Label (M("SAVEDLG_SUBSAMP") + Glib::ustring(":")) );
    setExpandAlignProperties(jpegSubSampLabel, true, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);

    jpegSubSamp = Gtk::manage (new MyComboBoxText ());
    setExpandAlignProperties(jpegSubSamp, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    jpegSubSamp->append (M("SAVEDLG_SUBSAMP_1"));
    jpegSubSamp->append (M("SAVEDLG_SUBSAMP_2"));
    jpegSubSamp->append (M("SAVEDLG_SUBSAMP_3"));
    jpegSubSamp->set_tooltip_text (M("SAVEDLG_SUBSAMP_TOOLTIP"));
    jpegSubSamp->set_active (2);
    jpegSubSamp->signal_changed().connect( sigc::mem_fun(*this, &SaveFormatPanel::formatChanged) );

    jpegOpts->attach(*jpegQual, 0, 0, 1, 2);
    jpegOpts->attach(*jpegSubSampLabel, 1, 0, 1, 1);
    jpegOpts->attach(*jpegSubSamp, 1, 1, 1, 1);
    jpegOpts->show_all ();

    // ---------------------  TIFF OPTIONS


    tiffUncompressed = new Gtk::CheckButton (M("SAVEDLG_TIFFUNCOMPRESSED"));
    setExpandAlignProperties(tiffUncompressed, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    tiffUncompressed->signal_toggled().connect( sigc::mem_fun(*this, &SaveFormatPanel::formatChanged));
    tiffUncompressed->show_all();


    // ---------------------  MAIN BOX


    savesPP = Gtk::manage (new Gtk::CheckButton (M("SAVEDLG_SAVESPP")));
    savesPP->signal_toggled().connect( sigc::mem_fun(*this, &SaveFormatPanel::formatChanged));
    savesPP->show_all();

    set_column_spacing(5);
    set_row_spacing(5);

    attach (*hb1, 0, 0, 1, 1);
    attach (*jpegOpts, 0, 1, 1, 1);
    attach (*tiffUncompressed, 0, 2, 1, 1);
    attach (*savesPP, 0, 4, 1, 2);
}
SaveFormatPanel::~SaveFormatPanel ()
{
    delete jpegQual;
    delete tiffUncompressed;
}

void SaveFormatPanel::init (SaveFormat &sf)
{

    FormatChangeListener* tmp = listener;
    listener = nullptr;

    if (sf.format == "jpg") {
        format->set_active (0);
    } else if (sf.format == "png" && sf.pngBits == 16) {
        format->set_active (5);
    } else if (sf.format == "png" && sf.pngBits == 8) {
        format->set_active (4);
    } else if (sf.format == "tif" && sf.tiffBits == 32) {
        format->set_active (3);
    } else if (sf.format == "tif" && sf.tiffBits == 16) {
        format->set_active (2);
    } else if (sf.format == "tif" && sf.tiffBits == 8) {
        format->set_active (1);
    }

    jpegSubSamp->set_active (sf.jpegSubSamp - 1);

    jpegQual->setValue (sf.jpegQuality);
    savesPP->set_active (sf.saveParams);
    tiffUncompressed->set_active (sf.tiffUncompressed);
    listener = tmp;
}

SaveFormat SaveFormatPanel::getFormat ()
{

    SaveFormat sf;

    int sel = format->get_active_row_number();
    sf.format = fstr[sel];

    if (sel == 5) {
        sf.pngBits = 16;
    } else {
        sf.pngBits = 8;
    }

    if (sel == 2) {
        sf.tiffBits = 16;
    } else if (sel == 3) {
        sf.tiffBits = 32;
    } else {
        sf.tiffBits = 8;
    }

    sf.jpegQuality      = (int) jpegQual->getValue ();
    sf.jpegSubSamp      = jpegSubSamp->get_active_row_number() + 1;
    sf.tiffUncompressed = tiffUncompressed->get_active();
    sf.saveParams       = savesPP->get_active ();
    return sf;
}

void SaveFormatPanel::formatChanged ()
{

    int act = format->get_active_row_number();

    if (act < 0 || act > 4) {
        return;
    }

    Glib::ustring fr = fstr[act];

    if (fr == "jpg") {
        jpegOpts->show_all();
        tiffUncompressed->hide();
    } else if (fr == "png") {
        jpegOpts->hide();
        tiffUncompressed->hide();
    } else if (fr == "tif") {
        jpegOpts->hide();
        tiffUncompressed->show_all();
    }

    if (listener) {
        listener->formatChanged (fr);
    }
}

void SaveFormatPanel::adjusterChanged (Adjuster* a, double newval)
{

    int act = format->get_active_row_number();

    if (act < 0 || act > 4) {
        return;
    }

    if (listener) {
        listener->formatChanged (fstr[act]);
    }
}

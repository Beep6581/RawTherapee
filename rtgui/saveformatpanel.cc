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
#include <saveformatpanel.h>
#include <multilangmgr.h>
#include <guiutils.h>

SaveFormatPanel::SaveFormatPanel () : listener (NULL) {

    jpegqual = new Adjuster (M("SAVEDLG_JPEGQUAL"), 0, 100, 1, 100);
    jpegqual->setAdjusterListener (this);
    jpegqual->show ();
    pngcompr = new Adjuster (M("SAVEDLG_PNGCOMPR"), 0, 6, 1, 6);
    pngcompr->setAdjusterListener (this);
    pngcompr->show ();

    Gtk::HBox* hb1 = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* flab = Gtk::manage (new Gtk::Label (M("SAVEDLG_FILEFORMAT")+":"));
    hb1->pack_start (*flab, Gtk::PACK_SHRINK,4);
    format = Gtk::manage (new Gtk::ComboBoxText ());
    format->append_text ("JPEG (8 bit)");
    format->append_text ("TIFF (8 bit)");
    format->append_text ("TIFF (16 bit)");
    format->append_text ("PNG (8 bit)");
    format->append_text ("PNG (16 bit)");
    format->set_active (0);
    oformat = 0;
    format->signal_changed().connect( sigc::mem_fun(*this, &SaveFormatPanel::formatChanged) );
    hb1->pack_start (*format);
    pack_start (*hb1, Gtk::PACK_SHRINK, 4);

    formatopts = Gtk::manage (new Gtk::VBox ());
    formatopts->pack_start (*jpegqual, Gtk::PACK_SHRINK, 4);
    pack_start (*formatopts, Gtk::PACK_SHRINK, 4);
    
    savespp = Gtk::manage (new Gtk::CheckButton (M("SAVEDLG_SAVESPP")));
    pack_start (*savespp, Gtk::PACK_SHRINK, 4);

    show_all ();
    set_border_width (4);
    
    fstr[0] = "jpg";
    fstr[1] = "tif";
    fstr[2] = "tif";
    fstr[3] = "png";
    fstr[4] = "png";
}

void SaveFormatPanel::init (SaveFormat &sf) {
  
    FormatChangeListener* tmp = listener;
    listener = NULL;    
    
    if (sf.format=="jpg")
        format->set_active (0);
    else if (sf.format=="png" && sf.pngBits==16)
        format->set_active (4);
    else if (sf.format=="png" && sf.pngBits==8)
        format->set_active (3);
    else if (sf.format=="tif" && sf.tiffBits==16)
        format->set_active (2);
    else if (sf.format=="tif" && sf.tiffBits==8)
        format->set_active (1);
       
    pngcompr->setValue (sf.pngCompression);
    jpegqual->setValue (sf.jpegQuality);
    savespp->set_active (sf.saveParams);
    listener = tmp;
}
 
SaveFormat SaveFormatPanel::getFormat () {

    SaveFormat sf;

    int sel = format->get_active_row_number();
    sf.format = fstr[sel];
    if (sel==4)
        sf.pngBits = 16;
    else
        sf.pngBits = 8;
    if (sel==2)
        sf.tiffBits = 16;
    else
        sf.tiffBits = 8;
    sf.pngCompression   = (int) pngcompr->getValue ();
    sf.jpegQuality      = (int) jpegqual->getValue ();
    sf.saveParams       = savespp->get_active ();
    return sf;
}
        
void SaveFormatPanel::formatChanged () {

    if (oformat==0)
        removeIfThere (formatopts, jpegqual);
    else if (oformat==3 || oformat==4)
        removeIfThere (formatopts, pngcompr);

    int act = format->get_active_row_number();
    if (act<0 || act>4)
        return;
        
    Glib::ustring fr = fstr[act];
    if (fr=="jpg") 
        formatopts->pack_start (*jpegqual, Gtk::PACK_SHRINK,4);
    else if (fr=="png") 
        formatopts->pack_start (*pngcompr, Gtk::PACK_SHRINK,4);

    oformat = act;
  
    if (listener) 
        listener->formatChanged (fr);
}

void SaveFormatPanel::adjusterChanged (Adjuster* a, double newval) {

    formatChanged ();
}

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
#include <curveeditor.h>
#include <fstream>
#include <string>
#include <multilangmgr.h>

CurveEditor::CurveEditor () {

    curve = Gtk::manage (new MyCurve ());
    Gtk::AspectFrame* af = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
    af->add (*curve);
    curve->set_size_request (-1, 200);
    pack_start (*af, Gtk::PACK_EXPAND_WIDGET);
    
    Gtk::HBox* bbox = Gtk::manage (new Gtk::HBox ());
    
    linear = Gtk::manage (new Gtk::Button (M("CURVEEDITOR_LINEAR")));
    save = Gtk::manage (new Gtk::Button ());
    Gtk::Image* saveImg = Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON));
    saveImg->show ();
    save->add (*saveImg);
    load = Gtk::manage (new Gtk::Button ());
    Gtk::Image* loadImg = Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON));
    loadImg->show ();
    load->add (*loadImg);
    
    bbox->pack_start (*linear);
    bbox->pack_end (*save, Gtk::PACK_SHRINK, 4);
    bbox->pack_end (*load, Gtk::PACK_SHRINK, 4);
    
    pack_end (*bbox, Gtk::PACK_SHRINK, 2);
    show_all ();

    linear->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::linearPressed) );
    save->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::savePressed) );
    load->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::loadPressed) );
    
    linear->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLINEAR"));
    save->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
    load->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));
}

void CurveEditor::linearPressed () {

    std::vector<double> lcurve (5);
    lcurve[0] = 1.0;
    lcurve[1] = 0.0;
    lcurve[2] = 0.0;
    lcurve[3] = 1.0;
    lcurve[4] = 1.0;
    curve->setPoints (lcurve);
    curve->queue_draw ();
    curve->notifyListener ();
}

void CurveEditor::savePressed () {

    Gtk::FileChooserDialog dialog(M("CURVEEDITOR_SAVEDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_SAVE);
//    if (options.multiUser)
//        dialog.set_current_folder (Options::rtdir + "/" + options.profilePath);
//    else
//        dialog.set_current_folder (argv0 + "/" + options.profilePath);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-save"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("CURVEEDITOR_FILEDLGFILTERCURVE"));
    filter_pp.add_pattern("*.rtc");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("CURVEEDITOR_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);
 
    dialog.set_do_overwrite_confirmation (true);

    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {

        std::string fname = dialog.get_filename();

        bool hasext = true;
        int dotpos = fname.find_last_of ('.');
        if (dotpos==Glib::ustring::npos)
            hasext = false;
        int dirpos1 = fname.find_last_of ('/');
        if (dirpos1!=Glib::ustring::npos || dirpos1>dotpos)
            hasext = false;
        int dirpos2 = fname.find_last_of ('\\');
        if (dirpos2!=Glib::ustring::npos || dirpos2>dotpos)
            hasext = false;

        if (!hasext) 
            fname = fname + ".rtc";

        if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS)) {
            Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
            Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
            int response = msgd.run ();
            if (response==Gtk::RESPONSE_NO)
                return;
        }
    
        std::ofstream f (fname.c_str());
        std::vector<double> p = curve->getPoints ();
        int ix = 0;
        if (p[ix++]<0)
            f << "Linear\n";
        else
            f << "Spline\n";
        for (int i=0; i<p.size()/2; i++, ix+=2)
            f << p[ix] << ' ' << p[ix+1] << std::endl;
        f.close ();
    }
}

void CurveEditor::loadPressed () {

    Gtk::FileChooserDialog dialog(M("CURVEEDITOR_LOADDLGLABEL"), Gtk::FILE_CHOOSER_ACTION_OPEN);

    dialog.add_button(Gtk::StockID("gtk-cancel"), Gtk::RESPONSE_CANCEL);
    dialog.add_button(Gtk::StockID("gtk-open"), Gtk::RESPONSE_OK);

    Gtk::FileFilter filter_pp;
    filter_pp.set_name(M("CURVEEDITOR_FILEDLGFILTERCURVE"));
    filter_pp.add_pattern("*.rtc");
    dialog.add_filter(filter_pp);

    Gtk::FileFilter filter_any;
    filter_any.set_name(M("CURVEEDITOR_FILEDLGFILTERANY"));
    filter_any.add_pattern("*");
    dialog.add_filter(filter_any);
 
    int result = dialog.run();

    if (result==Gtk::RESPONSE_OK) {
        std::ifstream f (dialog.get_filename().c_str());
        if (f) {
            std::vector<double> p;
            std::string s;
            f >> s;
            if (s=="Linear")
                p.push_back (-1);
            else if (s=="Spline")
                p.push_back (1);
            else return;
            double x;
            while (f) {
                f >> x;
                if (f)
                    p.push_back (x);
            }
            curve->setPoints (p);
            curve->queue_draw ();
            curve->notifyListener ();
        }
    }
}

void CurveEditor::setCurve (const std::vector<double>& c) {

    if (c.size()>4) {
        curve->setPoints (c);
        curve->queue_draw ();
    }
    else
        linearPressed ();
}

std::vector<double> CurveEditor::getCurve () {

    return curve->getPoints ();
}


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
#include <guiutils.h>
#include <multilangmgr.h>

extern Glib::ustring argv0;

CurveEditor::CurveEditor () : cl(NULL), realized(false), curveTypeIx(Linear), activeParamControl(-1) {

    Gtk::HBox* tsbox = Gtk::manage (new Gtk::HBox ());
    Gtk::Label* tslab = Gtk::manage (new Gtk::Label (M("CURVEEDITOR_TYPE")));
    curveType = Gtk::manage (new Gtk::ComboBoxText ());
    curve_reset = Gtk::manage (new Gtk::Button ());
    curve_reset->add (*Gtk::manage (new Gtk::Image (argv0+"/images/undo.png")));
    curve_reset->set_relief (Gtk::RELIEF_NONE);
    curve_reset->set_border_width (0);
    curve_reset->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLINEAR"));

    tsbox->pack_start (*tslab, Gtk::PACK_SHRINK, 8);
    tsbox->pack_start (*curveType);
    tsbox->pack_start (*curve_reset, Gtk::PACK_SHRINK, 0);
    
    pack_start (*tsbox);

    // Order set in the same order than "enum CurveType". Shouldn't change, for compatibility reason
    curveType->append_text (M("CURVEEDITOR_LINEAR"));		// 0 Linear
    curveType->append_text (M("CURVEEDITOR_CUSTOM"));		// 1 Spline
    curveType->append_text (M("CURVEEDITOR_PARAMETRIC"));	// 2 Parametric
    curveType->append_text (M("CURVEEDITOR_NURBS"));		// 3 NURBS
    curveType->set_active (Linear);

    curve_reset->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::curveResetPressed) );

    // custom curve
    customCurveBox = new Gtk::VBox ();
    Gtk::HBox* tmpa = Gtk::manage (new Gtk::HBox ());
    customCurve = Gtk::manage (new MyCurve ());
    //Gtk::AspectFrame* af = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
    //af->add (*customCurve);
    customCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
    customCurve->setType (Spline);
    //customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
    tmpa->pack_start (*customCurve, true, false, 4);
    customCurveBox->pack_start (*tmpa, true, true,4);

    Gtk::HBox* custombbox = Gtk::manage (new Gtk::HBox ());
    saveCustom = Gtk::manage (new Gtk::Button ());
    saveCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
    loadCustom = Gtk::manage (new Gtk::Button ());
    loadCustom->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

    custombbox->pack_end (*saveCustom, Gtk::PACK_EXPAND_WIDGET, 4);
    custombbox->pack_end (*loadCustom, Gtk::PACK_EXPAND_WIDGET, 4);

    customCurveBox->pack_end (*custombbox, Gtk::PACK_SHRINK, 2);
    customCurveBox->show_all ();

    saveCustom->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::savePressed) );
    loadCustom->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::loadPressed) );
    saveCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
    loadCustom->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

    // NURBS curve
    NURBSCurveBox = new Gtk::VBox ();
    Gtk::HBox* tmpb = Gtk::manage (new Gtk::HBox ());
    NURBSCurve = Gtk::manage (new MyCurve ());
    //Gtk::AspectFrame* af = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
    //af->add (*customCurve);
    NURBSCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
    NURBSCurve->setType (NURBS);
    //customCurve->set_tooltip_text (M("CURVEEDITOR_TOOLTIPMOVESPEED"));
    tmpb->pack_start (*NURBSCurve, true, false, 4);
    NURBSCurveBox->pack_start (*tmpb, true, true,4);

    Gtk::HBox* NURBSbbox = Gtk::manage (new Gtk::HBox ());
    saveNURBS = Gtk::manage (new Gtk::Button ());
    saveNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-save"), Gtk::ICON_SIZE_BUTTON)));
    loadNURBS = Gtk::manage (new Gtk::Button ());
    loadNURBS->add (*Gtk::manage (new Gtk::Image (Gtk::StockID("gtk-open"), Gtk::ICON_SIZE_BUTTON)));

    NURBSbbox->pack_end (*saveNURBS, Gtk::PACK_EXPAND_WIDGET, 4);
    NURBSbbox->pack_end (*loadNURBS, Gtk::PACK_EXPAND_WIDGET, 4);

    NURBSCurveBox->pack_end (*NURBSbbox, Gtk::PACK_SHRINK, 2);
    NURBSCurveBox->show_all ();

    saveNURBS->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::savePressed) );
    loadNURBS->signal_clicked().connect( sigc::mem_fun(*this, &CurveEditor::loadPressed) );
    saveNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPSAVE"));
    loadNURBS->set_tooltip_text (M("CURVEEDITOR_TOOLTIPLOAD"));

    // parametric curve
    paramCurveBox = new Gtk::VBox ();
    paramCurve = Gtk::manage (new MyCurve ());
    Gtk::Table* paramctab = Gtk::manage (new Gtk::Table (2,1));
    //Gtk::AspectFrame* afp = Gtk::manage (new Gtk::AspectFrame ("",Gtk::ALIGN_CENTER,Gtk::ALIGN_CENTER,1,false));
    //afp->add (*paramCurve);
    paramCurve->set_size_request (GRAPH_SIZE+2*RADIUS, GRAPH_SIZE+2*RADIUS);
    paramCurve->setType (Parametric);
    shcSelector = Gtk::manage (new SHCSelector ());
    shcSelector->set_size_request (GRAPH_SIZE, 20);

    paramctab->attach (*paramCurve, 0, 1, 0, 1, Gtk::FILL, Gtk::SHRINK, 4, 4);
    paramctab->attach (*shcSelector, 0, 1, 1, 2, Gtk::FILL, Gtk::SHRINK, RADIUS+4, 0);

    Gtk::HBox* tmpc = Gtk::manage (new Gtk::HBox ());
    tmpc->pack_start (*paramctab, true, false);

    paramCurveBox->pack_start (*tmpc, true, true);

    highlights = Gtk::manage (new Adjuster (M("CURVEEDITOR_HIGHLIGHTS"), -100, 100, 1, 0));
    lights     = Gtk::manage (new Adjuster (M("CURVEEDITOR_LIGHTS"), -100, 100, 1, 0));
    darks      = Gtk::manage (new Adjuster (M("CURVEEDITOR_DARKS"), -100, 100, 1, 0));
    shadows    = Gtk::manage (new Adjuster (M("CURVEEDITOR_SHADOWS"), -100, 100, 1, 0));

    Gtk::EventBox* evhighlights = Gtk::manage (new Gtk::EventBox ());
    Gtk::EventBox* evlights = Gtk::manage (new Gtk::EventBox ());
    Gtk::EventBox* evdarks = Gtk::manage (new Gtk::EventBox ());
    Gtk::EventBox* evshadows = Gtk::manage (new Gtk::EventBox ());
    
    evhighlights->add (*highlights);
    evlights->add (*lights);
    evdarks->add (*darks);
    evshadows->add (*shadows);
    
    paramCurveBox->pack_start (*Gtk::manage (new Gtk::HSeparator ()));
    paramCurveBox->pack_start (*evhighlights);
    paramCurveBox->pack_start (*evlights);
    paramCurveBox->pack_start (*evdarks);
    paramCurveBox->pack_start (*evshadows);
    paramCurveBox->show_all ();
    
    customCurveBox->reference ();
    paramCurveBox->reference ();
    
    customCurve->setCurveListener (this); 
    NURBSCurve->setCurveListener (this);
    paramCurve->setCurveListener (this);
    shcSelector->setSHCListener (this);
    
    highlights->setAdjusterListener (this);
    lights->setAdjusterListener (this);
    darks->setAdjusterListener (this);
    shadows->setAdjusterListener (this);

    evhighlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK); 
    evlights->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK); 
    evdarks->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK); 
    evshadows->set_events(Gdk::ENTER_NOTIFY_MASK | Gdk::LEAVE_NOTIFY_MASK); 
    typeconn = curveType->signal_changed().connect (sigc::mem_fun(*this, &CurveEditor::typeSelectionChanged) );
    evhighlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterEntered), 4));
    evlights->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterEntered), 5));
    evdarks->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterEntered), 6));
    evshadows->signal_enter_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterEntered), 7));
    evhighlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterLeft), 4));
    evlights->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterLeft), 5));
    evdarks->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterLeft), 6));
    evshadows->signal_leave_notify_event().connect (sigc::bind(sigc::mem_fun(*this, &CurveEditor::adjusterLeft), 7));

    show_all ();
}

CurveEditor::~CurveEditor () {

    delete customCurveBox;
    delete paramCurveBox;
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

        if (getExtension (fname)!="rtc") 
            fname = fname + ".rtc";

        if (Glib::file_test (fname, Glib::FILE_TEST_EXISTS)) {
            Glib::ustring msg_ = Glib::ustring("<b>") + fname + ": " + M("MAIN_MSG_ALREADYEXISTS") + "\n" + M("MAIN_MSG_QOVERWRITE") + "</b>";
            Gtk::MessageDialog msgd (msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
            int response = msgd.run ();
            if (response==Gtk::RESPONSE_NO)
                return;
        }
    
        std::ofstream f (fname.c_str());
        std::vector<double> p = customCurve->getPoints ();

        switch (curveType->get_active_row_number()) {
        case Spline:		// custom
        	p = customCurve->getPoints ();
        	break;
        case NURBS:		// NURBS
        	p = NURBSCurve->getPoints ();
        	break;
        default:
        	break;
        }

        int ix = 0;
        if (p[ix]==(double)(Linear))
            f << "Linear\n";
        else if (p[ix]==(double)(Spline))
            f << "Spline\n";
        else if (p[ix]==(double)(NURBS))
            f << "NURBS\n";
        ix++;
        for (unsigned int i=0; i<p.size()/2; i++, ix+=2)
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
                p.push_back ((double)(Linear));
            else if (s=="Spline")
                p.push_back ((double)(Spline));
            else if (s=="NURBS")
                p.push_back ((double)(NURBS));
            else return;
            double x;
            while (f) {
                f >> x;
                if (f)
                    p.push_back (x);
            }
            if (p[0] == (double)(Spline)) {
				customCurve->setPoints (p);
				customCurve->queue_draw ();
				customCurve->notifyListener ();
            }
            else if (p[0] == (double)(NURBS)) {
				NURBSCurve->setPoints (p);
				NURBSCurve->queue_draw ();
				NURBSCurve->notifyListener ();
            }
        }
    }
}

void CurveEditor::on_realize () {
    
    Gtk::VBox::on_realize();
    realized = true; 
    setCurve (tmpCurve);
}

void CurveEditor::setCurve (const std::vector<double>& c) {

    tmpCurve = c;
    
    if (realized && curveType->get_active_row_number()<=Unchanged) { // if it is not realized or "unchanged" is selected, just store the curve (prev line) and do not change gui

        typeconn.block(true);
        if (c.size()==0 || c[0]==(double)(Linear)) {
            curveType->set_active (Linear);		// Change the combo selection
            curveTypeIx = Linear;
        }
        else if (c[0]==(double)(Spline)) {
            curveType->set_active (Spline);
            curveTypeIx = Spline;
            customCurve->setPoints (c);
        }
        else if (c[0]==(double)(Parametric)) {
            curveType->set_active (Parametric);
            curveTypeIx = Parametric;
            paramCurve->setPoints (c);
            shcSelector->setPositions (c[1], c[2], c[3]);
            highlights->setValue (c[4]);
            lights->setValue (c[5]);
            darks->setValue (c[6]);
            shadows->setValue (c[7]);
        }
        else if (c[0]==(double)(NURBS)) {
            curveType->set_active (NURBS);
            curveTypeIx = NURBS;
            NURBSCurve->setPoints (c);
        }
        removeIfThere (this, customCurveBox, false);
        removeIfThere (this, paramCurveBox, false);
        removeIfThere (this, NURBSCurveBox, false);

        if (curveType->get_active_row_number()==Spline)
            pack_start (*customCurveBox);
        else if (curveType->get_active_row_number()==Parametric)
            pack_start (*paramCurveBox);   
        else if (curveType->get_active_row_number()==NURBS)
            pack_start (*NURBSCurveBox);

        typeconn.block(false);
    }
}

std::vector<double> CurveEditor::getCurve () {

    if (!realized || curveType->get_active_row_number()==Unchanged)
        return tmpCurve;

    // linear
    if (curveTypeIx<=Linear) {
        std::vector<double> lcurve (1);
        lcurve[0] = (double)(Linear);
        return lcurve;
    }
    // parametric
    else if (curveTypeIx==Parametric) {
        std::vector<double> lcurve (8);
        lcurve[0] = (double)(Parametric);
        shcSelector->getPositions (lcurve[1], lcurve[2], lcurve[3]);
        lcurve[4] = highlights->getValue ();
        lcurve[5] = lights->getValue ();
        lcurve[6] = darks->getValue ();
        lcurve[7] = shadows->getValue ();
        return lcurve;
    }
    // spline (custom)
    else if (curveTypeIx==Spline)
        return customCurve->getPoints ();
    // NURBS (control cage)
    else // Default solution, we return the NURBS curve (curveTypeIx==NURBS)
        return NURBSCurve->getPoints ();
}

void CurveEditor::typeSelectionChanged () {

    removeIfThere (this, customCurveBox, false);
    removeIfThere (this, paramCurveBox, false);
    removeIfThere (this, NURBSCurveBox, false);

    if (curveType->get_active_row_number()==Spline) {
        pack_start (*customCurveBox);
    }
    else if (curveType->get_active_row_number()==Parametric) {
        pack_start (*paramCurveBox);   
    }
    else if (curveType->get_active_row_number()==NURBS) {
        pack_start (*NURBSCurveBox);
    }

    if (curveType->get_active_row_number() < Unchanged)
        curveTypeIx = (CurveType)curveType->get_active_row_number();

    curveChanged ();
}

void CurveEditor::curveChanged () {

    if (cl) {
    	if (cl->isMulti())
    		cl->curveChanged (this);
    	else
    		cl->curveChanged ();
    }
}

void CurveEditor::curveResetPressed () {
	switch (curveTypeIx) {
	case NURBS :	// = Control cage
		NURBSCurve->reset ();
		break;
	case Spline :	// = Custom
		customCurve->reset ();
		break;
	case Parametric :
	    highlights->resetPressed();
	    lights->resetPressed();
	    darks->resetPressed();
	    shadows->resetPressed();
	    shcSelector->reset();
	    paramCurve->reset ();
		break;
	default:
		break;
	}
	curveChanged ();
}

void CurveEditor::shcChanged () {

    paramCurve->setPoints (getCurve());
    if (cl)
        cl->curveChanged ();
}

void CurveEditor::adjusterChanged (Adjuster* a, double newval) {

    paramCurve->setPoints (getCurve());
    if (cl)
        cl->curveChanged ();
}

bool CurveEditor::adjusterEntered (GdkEventCrossing* ev, int ac) {

    if (ev->detail != GDK_NOTIFY_INFERIOR) {    
        activeParamControl = ac;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

bool CurveEditor::adjusterLeft (GdkEventCrossing* ev, int ac) {
    
    if (ev->detail != GDK_NOTIFY_INFERIOR) {    
        activeParamControl = -1;
        paramCurve->setActiveParam (activeParamControl);
    }
    return true;
}

void CurveEditor::setBatchMode (bool batchMode) {

  curveType->append_text (M("GENERAL_UNCHANGED"));
}

bool CurveEditor::isUnChanged () {

    return curveType->get_active_row_number()==Unchanged;
}

void CurveEditor::setUnChanged (bool uc) {

    if (uc) {
        typeconn.block(true);
        removeIfThere (this, customCurveBox, false);
        removeIfThere (this, paramCurveBox, false);
        removeIfThere (this, NURBSCurveBox, false);
        curveType->set_active (Unchanged);
        typeconn.block(false);
    }
    else {
        typeconn.block(true);
        curveType->set_active (-1); // hack: if it remains 3 (unchanged), then setcurve does not switch selection in the combo
        setCurve (getCurve ());
        typeconn.block(false);
    }
}

void CurveEditor::updateBackgroundHistogram (unsigned int* hist) {

    paramCurve->updateBackgroundHistogram (hist);
    customCurve->updateBackgroundHistogram (hist);
    NURBSCurve->updateBackgroundHistogram (hist);
}

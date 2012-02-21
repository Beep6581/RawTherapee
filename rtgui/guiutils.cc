/*
 *  This file is part of RawTherapee.
 *
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
#include "guiutils.h"
#include "options.h"
#include "../rtengine/utils.h"
#include "rtimage.h"

bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference) {

    Glib::ListHandle<Gtk::Widget*> list = cont->get_children ();
    Glib::ListHandle<Gtk::Widget*>::iterator i = list.begin ();
    for (; i!=list.end() && *i!=w; i++);
    if (i!=list.end()) {
        if (increference)
            w->reference ();
        cont->remove (*w);
        return true;
    }
    else
        return false;
}

void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh) {

    if (options.thumbInterp==0)
        rtengine::nearestInterp (src, sw, sh, dst, dw, dh);
    else if (options.thumbInterp==1)
        rtengine::bilinearInterp (src, sw, sh, dst, dw, dh);
}

Glib::ustring removeExtension (const Glib::ustring& filename) {

    Glib::ustring bname = Glib::path_get_basename(filename);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        return filename.substr (0, filename.size()-(bname.size()-lastdot));
    else
        return filename;
}

Glib::ustring getExtension (const Glib::ustring& filename) {

    Glib::ustring bname = Glib::path_get_basename(filename);
    int lastdot = bname.find_last_of ('.');
    if (lastdot!=bname.npos)
        return filename.substr (filename.size()-(bname.size()-lastdot)+1, filename.npos);
    else
        return "";
}

void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams) {

    cr->set_line_width (0.);
    cr->rectangle (imx, imy, imw, imh);
    cr->clip ();

    double c1x = (cparams.x-startx)*scale;
    double c1y = (cparams.y-starty)*scale;
    double c2x = (cparams.x+cparams.w-1-startx)*scale;
    double c2y = (cparams.y+cparams.h-1-starty)*scale;

    cr->set_source_rgba (options.cutOverlayBrush[0], options.cutOverlayBrush[1], options.cutOverlayBrush[2], options.cutOverlayBrush[3]);

    cr->rectangle (imx, imy, imw, c1y);
    cr->rectangle (imx, imy+c2y, imw, imh-c2y);
    cr->rectangle (imx, imy+c1y, c1x, c2y-c1y+1);
    cr->rectangle (imx+c2x, imy+c1y, imw-c2x, c2y-c1y+1);
    cr->fill ();

    // rectangle around the cropped area and guides
    if (cparams.guide!="None") {
        double rectx1 = c1x + imx + 0.5;
        double recty1 = c1y + imy + 0.5;
        double rectx2 = c2x + imx + 0.5;
        double recty2 = c2y + imy + 0.5;
        cr->set_line_width (1.0);
        cr->set_source_rgb (1.0, 1.0, 1.0);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        cr->set_source_rgb (0.0, 0.0, 0.0);
        std::valarray<double> ds (1);
        ds[0] = 4;
        cr->set_dash (ds, 0);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        ds.resize (0);
        cr->set_dash (ds, 0);
        
        if (cparams.guide!="Rule of diagonals") {
            // draw guide lines
            std::vector<double> horiz_ratios;
            std::vector<double> vert_ratios;
            
            if (cparams.guide=="Rule of thirds") {
                horiz_ratios.push_back (1.0/3.0);
                horiz_ratios.push_back (2.0/3.0);
                vert_ratios.push_back (1.0/3.0);
                vert_ratios.push_back (2.0/3.0);
            }
            else if (cparams.guide=="Harmonic means 1") {
                horiz_ratios.push_back (1.0-0.618);
                vert_ratios.push_back (1.0-0.618);
            }
            else if (cparams.guide=="Harmonic means 2") {
                horiz_ratios.push_back (0.618);
                vert_ratios.push_back (1.0-0.618);
            }
            else if (cparams.guide=="Harmonic means 3") {
                horiz_ratios.push_back (1.0-0.618);
                vert_ratios.push_back (0.618);
            }
            else if (cparams.guide=="Harmonic means 4") {
                horiz_ratios.push_back (0.618);
                vert_ratios.push_back (0.618);
            } 
            else if (cparams.guide=="Grid") {
                // To have even distribution, normalize it a bit
                const int longSideNumLines=10;

                int w=rectx2-rectx1, h=recty2-recty1, shortSideNumLines;
                if (w>longSideNumLines && h>longSideNumLines) {
                    if (w>h) {
                        for (int i=1;i<longSideNumLines;i++) vert_ratios.push_back ((double)i/longSideNumLines);
                    
                        shortSideNumLines=(int)round(h*(double)longSideNumLines/w);
                        for (int i=1;i<shortSideNumLines;i++) horiz_ratios.push_back ((double)i/shortSideNumLines);
                    } else {
                        for (int i=1;i<longSideNumLines;i++) horiz_ratios.push_back ((double)i/longSideNumLines);
                    
                        shortSideNumLines=(int)round(w*(double)longSideNumLines/h);
                        for (int i=1;i<shortSideNumLines;i++) vert_ratios.push_back ((double)i/shortSideNumLines);
                    }
                }
            }   
            else if (cparams.guide=="ePassport") {
                /* Official measurements do not specify exact ratios, just min/max measurements within which the eyes and chin-crown distance must lie. I averaged those measurements to produce these guides.
                 * The first horizontal guide is for the crown, the second is rougly for the nostrils, the third is for the chin.
                 * http://www.homeoffice.gov.uk/agencies-public-bodies/ips/passports/information-photographers/
                 * "(...) the measurement of the face from the bottom of the chin to the crown (ie the top of the head, not the top of the hair) is between 29mm and 34mm."
                 */
                horiz_ratios.push_back (7.0/45.0);
                horiz_ratios.push_back (26.0/45.0);
                horiz_ratios.push_back (37.0/45.0);
                vert_ratios.push_back (0.5);
            }

            // Horizontals
            for (int i=0; i<horiz_ratios.size(); i++) {
                cr->set_source_rgb (1.0, 1.0, 1.0);
                cr->move_to (rectx1, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->line_to (rectx2, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->stroke ();
                cr->set_source_rgb (0.0, 0.0, 0.0);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (rectx1, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->line_to (rectx2, recty1 + (recty2-recty1) * horiz_ratios[i]);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }
            // Verticals
            for (int i=0; i<vert_ratios.size(); i++) {
                cr->set_source_rgb (1.0, 1.0, 1.0);
                cr->move_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty1);
                cr->line_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty2);
                cr->stroke ();
                cr->set_source_rgb (0.0, 0.0, 0.0);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty1);
                cr->line_to (rectx1 + (rectx2-rectx1) * vert_ratios[i], recty2);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }           
        }
        else {
            int corners_from[4][2];
            int corners_to[4][2];
            int mindim = MIN (rectx2-rectx1+1, recty2-recty1+1);
            corners_from[0][0] = rectx1;
            corners_from[0][1] = recty1;
            corners_to[0][0]   = rectx1 + mindim;
            corners_to[0][1]   = recty1 + mindim;
            corners_from[1][0] = rectx1;
            corners_from[1][1] = recty2;
            corners_to[1][0]   = rectx1 + mindim;
            corners_to[1][1]   = recty2 - mindim;
            corners_from[2][0] = rectx2;
            corners_from[2][1] = recty1;
            corners_to[2][0]   = rectx2 - mindim;
            corners_to[2][1]   = recty1 + mindim;
            corners_from[3][0] = rectx2;
            corners_from[3][1] = recty2;
            corners_to[3][0]   = rectx2 - mindim;
            corners_to[3][1]   = recty2 - mindim;
            for (int i=0; i<4; i++) {
                cr->set_source_rgb (1.0, 1.0, 1.0);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                cr->set_source_rgb (0.0, 0.0, 0.0);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }
        }
    }
    cr->reset_clip ();
}


/*
 *
 * Derived class of some widgets to properly handle the scroll wheel ;
 * the user has to use the Shift key to be able to change the widget's value,
 * otherwise the mouse wheel will scroll the editor's tabs content.
 *
 */
MyScrolledWindow::MyScrolledWindow () {
	set_size_request(-1,30);
}

bool MyScrolledWindow::on_scroll_event (GdkEventScroll* event) {
	if (!options.hideTPVScrollbar) {
		Gtk::ScrolledWindow::on_scroll_event (event);
		return true;
	}

    Gtk::Adjustment *adjust = get_vadjustment();
    Gtk::VScrollbar *scroll = get_vscrollbar();
    if (adjust && scroll) {
    	double upper = adjust->get_upper();
    	double lower = adjust->get_lower();
    	double value = adjust->get_value();
    	double step  = adjust->get_step_increment();
    	double value2 = 0.;
    	if (event->direction == GDK_SCROLL_DOWN) {
    		value2 = value+step;
    		if (value2 > upper)
    			value2 = upper;
    		if (value2 != value) {
    			scroll->set_value(value2);
    		}
    	}
    	else {
    		value2 = value-step;
    		if (value2 < lower)
    			value2 = lower;
    		if (value2 != value) {
    			scroll->set_value(value2);
    		}
    	}
    }
    return true;
}

MyComboBoxText::MyComboBoxText () {
	set_size_request(40, -1);
}

bool MyComboBoxText::on_scroll_event (GdkEventScroll* event) {

	// If Shift is pressed, the widget is modified
	if (event->state & GDK_SHIFT_MASK) {
		Gtk::ComboBoxText::on_scroll_event(event);
		return true;
	}
	// ... otherwise the scroll event is sent back to an upper level
	return false;
}

MyComboBox::MyComboBox () {
	set_size_request(40, -1);
}

bool MyComboBox::on_scroll_event (GdkEventScroll* event) {

	// If Shift is pressed, the widget is modified
	if (event->state & GDK_SHIFT_MASK) {
		Gtk::ComboBox::on_scroll_event(event);
		return true;
	}
	// ... otherwise the scroll event is sent back to an upper level
	return false;
}

MySpinButton::MySpinButton () {
	Gtk::Border border;
	border.bottom = 0;
	border.top = 0;
	border.left = 3;
	border.right = 3;
	set_inner_border(border);
	set_numeric(true);
	set_wrap(true);
	set_alignment(Gtk::ALIGN_RIGHT);
}

void MySpinButton::updateSize() {
	double vMin, vMax;
	double step, page;
	int maxAbs;
	unsigned int digits, digits2;
	unsigned int maxLen;

	get_range(vMin, vMax);
	get_increments (step, page);

	digits = get_digits();
	maxAbs = (int)(fmax(fabs(vMin), fabs(vMax))+0.000001);
	if (maxAbs==0)
		digits2 = 1;
	else {
		digits2 = (int)(log10(double(maxAbs))+0.000001);
		digits2++;
	}
	maxLen = digits+digits2+(vMin<0?1:0)+(digits>0?1:0);
	set_max_length(maxLen);
	set_width_chars(maxLen);
}

bool MySpinButton::on_key_press_event (GdkEventKey* event) {
	bool rcode = Gtk::Widget::on_key_press_event(event);
	if ( (event->string[0] >= 'a' && event->string[0] <= 'z')
	   ||(event->string[0] >= 'A' && event->string[0] <= 'Z')
	   || event->string[0] == '+'
	)
		return false;
	return rcode;
}

bool MySpinButton::on_scroll_event (GdkEventScroll* event) {
	// If Shift is pressed, the widget is modified
	if (event->state & GDK_SHIFT_MASK) {
		Gtk::SpinButton::on_scroll_event(event);
		return true;
	}
	// ... otherwise the scroll event is sent back to an upper level
	return false;
}

bool MyHScale::on_scroll_event (GdkEventScroll* event) {

	// If Shift is pressed, the widget is modified
	if (event->state & GDK_SHIFT_MASK) {
		Gtk::HScale::on_scroll_event(event);
		return true;
	}
	// ... otherwise the scroll event is sent back to an upper level
	return false;
}

MyFileChooserButton::MyFileChooserButton (const Glib::ustring& title, Gtk::FileChooserAction action) : Gtk::FileChooserButton(title, action) {
	set_size_request(20, -1);
};

// For an unknown reason (a bug ?), it doesn't work when action = FILE_CHOOSER_ACTION_SELECT_FOLDER !
bool MyFileChooserButton::on_scroll_event (GdkEventScroll* event) {

	// If Shift is pressed, the widget is modified
	if (event->state & GDK_SHIFT_MASK) {
		Gtk::FileChooserButton::on_scroll_event(event);
		return true;
	}
	// ... otherwise the scroll event is sent back to an upper level
	return false;
}

TextOrIcon::TextOrIcon (Glib::ustring fname, Glib::ustring labelTx, Glib::ustring tooltipTx, TOITypes type) {

	imgIcon = 0;
	label = 0;
	filename = fname;
	labelText = labelTx;
	tooltipText = tooltipTx;

	switchTo(type);
}

TextOrIcon::~TextOrIcon () {
	if (imgIcon) delete imgIcon;
	if (label) delete label;
}

void TextOrIcon::switchTo(TOITypes type) {
	switch (type) {
	case (TOI_ICON):
		if (!imgIcon) {
			removeIfThere(this, label, false);
			delete label;
			label = 0;
			imgIcon = new RTImage (filename);
			pack_start(*imgIcon, Gtk::PACK_SHRINK, 0);
			set_tooltip_markup ("<span font_size=\"large\" font_weight=\"bold\">" + labelText  + "</span>\n" + tooltipText);
		}
		// do nothing if imgIcon exist, which mean that it is currently being displayed
		break;
	case(TOI_TEXT):
	default:
		if (!label) {
			removeIfThere(this, imgIcon, false);
			delete imgIcon;
			imgIcon = 0;
			label = new Gtk::Label (labelText, Gtk::ALIGN_CENTER);
			pack_start(*label, Gtk::PACK_EXPAND_WIDGET, 0);
			set_tooltip_markup (tooltipText);
		}
		// do nothing if label exist, which mean that it is currently being displayed
		break;
	}
	show_all();
}

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
#include <navigator.h>

Navigator::Navigator () {

	set_label ("Navigator");
	Gtk::VBox* mbox = Gtk::manage (new Gtk::VBox ());
	previewWindow = Gtk::manage (new PreviewWindow ());
	mbox->pack_start (*previewWindow, Gtk::PACK_SHRINK, 2);
	position = Gtk::manage (new Gtk::Label ());
	mbox->pack_start (*position, Gtk::PACK_SHRINK, 2);

	R = Gtk::manage (new Gtk::Label ());
	G = Gtk::manage (new Gtk::Label ());
	B = Gtk::manage (new Gtk::Label ());
	H = Gtk::manage (new Gtk::Label ());
	S = Gtk::manage (new Gtk::Label ());
	V = Gtk::manage (new Gtk::Label ());

    Gtk::Table* table = new Gtk::Table (3, 2);
    table->attach (*R, 0, 1, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*G, 0, 1, 1, 2, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*B, 0, 1, 2, 3, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*H, 1, 2, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*S, 1, 2, 1, 2, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*V, 1, 2, 2, 3, Gtk::EXPAND, Gtk::SHRINK, 0, 0);

	mbox->pack_start (*table, Gtk::PACK_SHRINK, 2);
	
	add (*mbox);
	
	setInvalid ();
	show_all ();
}

void Navigator::setInvalid () {

	position->set_text ("x = n/a, y = n/a");
	R->set_text ("R = n/a");
	G->set_text ("G = n/a");
	B->set_text ("B = n/a");
	H->set_text ("H = n/a");
	S->set_text ("S = n/a");
	V->set_text ("V = n/a");
}

void Navigator::pointerMoved (bool validPos, int x, int y, int r, int g, int b) {

	if (!validPos)
		setInvalid ();
	else {
		position->set_text (Glib::ustring::compose ("x = %1, y = %2", x, y));
		R->set_text (Glib::ustring::compose ("R = %1", r));
		G->set_text (Glib::ustring::compose ("G = %1", g));
		B->set_text (Glib::ustring::compose ("B = %1", b));
		int h, s, v;
		rgb2hsv (r, g, b, h, s, v);
		H->set_text (Glib::ustring::compose ("H = %1", h));
		S->set_text (Glib::ustring::compose ("S = %1", s));
		V->set_text (Glib::ustring::compose ("V = %1", v));
	}
}

void Navigator::rgb2hsv (int r, int g, int b, int &h, int &s, int &v) {

	volatile double var_R = r / 255.0;
	volatile double var_G = g / 255.0;
	volatile double var_B = b / 255.0;

	volatile double var_Min = MIN(MIN(var_R,var_G),var_B);
	volatile double var_Max = MAX(MAX(var_R,var_G),var_B);
	double del_Max = var_Max - var_Min;
    double	V = (var_Max + var_Min) / 2;
    double H, S;
	if (fabs(del_Max)<0.0000001) {
		H = 0;
		S = 0;
	}
	else {
		if (V < 0.5) S = del_Max / (var_Max + var_Min);
		else         S = del_Max / (2 - var_Max - var_Min);

		double del_R = ( ( ( var_Max - var_R ) / 6.0 ) + ( del_Max / 2.0 ) ) / del_Max;
		double del_G = ( ( ( var_Max - var_G ) / 6.0 ) + ( del_Max / 2.0 ) ) / del_Max;
		double del_B = ( ( ( var_Max - var_B ) / 6.0 ) + ( del_Max / 2.0 ) ) / del_Max;
		if      ( var_R == var_Max ) H = del_B - del_G; 
		else if ( var_G == var_Max ) H = (1.0 / 3.0) + del_R - del_B; 
		else if ( var_B == var_Max ) H = (2.0 / 3.0) + del_G - del_R; 

		if ( H < 0 )  H += 1;
		if ( H > 1 )  H -= 1;
	}
    
    h = (int)(H*255.0);
    s = (int)(S*255.0);
    v = (int)(V*255.0);
}

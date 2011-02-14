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
#include <toolpanel.h>

Navigator::Navigator () {

	set_label (M("MAIN_MSG_NAVIGATOR"));
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
	LAB_A = Gtk::manage (new Gtk::Label ());
	LAB_B = Gtk::manage (new Gtk::Label ());
	LAB_L = Gtk::manage (new Gtk::Label ());

    Gtk::Table* table = new Gtk::Table (3, 3);
    table->attach (*R, 0, 1, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*G, 0, 1, 1, 2, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*B, 0, 1, 2, 3, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*H, 1, 2, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*S, 1, 2, 1, 2, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*V, 1, 2, 2, 3, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*LAB_L, 2, 3, 0, 1, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*LAB_A, 2, 3, 1, 2, Gtk::EXPAND, Gtk::SHRINK, 0, 0);
    table->attach (*LAB_B, 2, 3, 2, 3, Gtk::EXPAND, Gtk::SHRINK, 0, 0);

	mbox->pack_start (*table, Gtk::PACK_SHRINK, 2);
	
	add (*mbox);
	
	setInvalid ();
	show_all ();
}

void Navigator::setInvalid () {

        position->set_text (M("NAVIGATOR_XY_NA"));
	R->set_text (M("NAVIGATOR_R_NA"));
	G->set_text (M("NAVIGATOR_G_NA"));
	B->set_text (M("NAVIGATOR_B_NA"));
	H->set_text (M("NAVIGATOR_H_NA"));
	S->set_text (M("NAVIGATOR_S_NA"));
	V->set_text (M("NAVIGATOR_V_NA"));
	LAB_A->set_text (M("NAVIGATOR_LAB_A_NA"));
	LAB_B->set_text (M("NAVIGATOR_LAB_B_NA"));
	LAB_L->set_text (M("NAVIGATOR_LAB_L_NA"));
}

void Navigator::pointerMoved (bool validPos, int x, int y, int r, int g, int b) {

	if (!validPos)
		setInvalid ();
	else {
		position->set_text (Glib::ustring::compose ("x = %1, y = %2", x, y));
		R->set_text (Glib::ustring::compose (M("NAVIGATOR_R_VALUE"), r));
		G->set_text (Glib::ustring::compose (M("NAVIGATOR_G_VALUE"), g));
		B->set_text (Glib::ustring::compose (M("NAVIGATOR_B_VALUE"), b));
		int h, s, v;
		rgb2hsv (r, g, b, h, s, v);
		H->set_text (Glib::ustring::compose (M("NAVIGATOR_H_VALUE"), h));
		S->set_text (Glib::ustring::compose (M("NAVIGATOR_S_VALUE"), s));
		V->set_text (Glib::ustring::compose (M("NAVIGATOR_V_VALUE"), v));
		int LAB_a, LAB_b, LAB_l;
		rgb2lab (r, g, b, LAB_l, LAB_a, LAB_b);
		LAB_A->set_text (Glib::ustring::compose (M("NAVIGATOR_LAB_A_VALUE"), LAB_a));
		LAB_B->set_text (Glib::ustring::compose (M("NAVIGATOR_LAB_B_VALUE"), LAB_b));
		LAB_L->set_text (Glib::ustring::compose (M("NAVIGATOR_LAB_L_VALUE"), LAB_l));

	}
}

void Navigator::rgb2hsv (int r, int g, int b, int &h, int &s, int &v) {

	volatile double var_R = r / 255.0;
	volatile double var_G = g / 255.0;
	volatile double var_B = b / 255.0;
	
	double var_Min = MIN(MIN(var_R,var_G),var_B);
	double var_Max = MAX(MAX(var_R,var_G),var_B);
	double del_Max = var_Max - var_Min;
	double V = var_Max;
	double H, S;
	if (fabs(del_Max)<0.001) {
		H = 0;
		S = 0;
	}
	else {
		S = del_Max/var_Max;
		
		if      ( var_R == var_Max ) H = (var_G - var_B)/del_Max; 
		else if ( var_G == var_Max ) H = 2.0 + (var_B - var_R)/del_Max; 
		else if ( var_B == var_Max ) H = 4.0 + (var_R - var_G)/del_Max; 
		H /= 6.0;
		
		if ( H < 0 )  H += 1;
		if ( H > 1 )  H -= 1;
	}
    
    h = (int)(H*255.0);
    s = (int)(S*255.0);
    v = (int)(V*255.0);
}

void Navigator::rgb2lab (int r, int g, int b, int &LAB_l, int &LAB_a, int &LAB_b) {

	volatile double var_R = r / 255.0;
	volatile double var_G = g / 255.0;
	volatile double var_B = b / 255.0;

	if ( var_R > 0.04045 ) 
		var_R = pow ( ( ( var_R + 0.055 ) / 1.055 ), 2.4);
	else    
		var_R = var_R / 12.92;
	if ( var_G > 0.04045 ) 
		var_G = pow ( ( ( var_G + 0.055 ) / 1.055 ), 2.4);
	else                   
		var_G = var_G / 12.92;
	if ( var_B > 0.04045 ) 
		var_B = pow ( ( ( var_B + 0.055 ) / 1.055 ), 2.4);
	else    
		var_B = var_B / 12.92;

	var_R = var_R * 100;
	var_G = var_G * 100;
	var_B = var_B * 100;

	double var_X = ( var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805 ) / 95.047;
	double var_Y = ( var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722 ) / 100.000;
	double var_Z = ( var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505 ) / 108.883;

	if ( var_X > 0.008856 ) 
		var_X = pow (var_X, ( 1.0/3.0 ));
	else
		var_X = ( 7.787 * var_X ) + ( 16.0 / 116.0 );
	if ( var_Y > 0.008856 ) 
		var_Y = pow (var_Y, ( 1.0/3.0 ));
	else 
		var_Y = ( 7.787 * var_Y ) + ( 16.0 / 116.0 );
	if ( var_Z > 0.008856 ) 
		var_Z = pow (var_Z, ( 1.0/3.0 ));
	else  
		var_Z = ( 7.787 * var_Z ) + ( 16.0 / 116.0 );

	LAB_l = ( 116 * var_Y ) - 16;
	LAB_a = 500 * ( var_X - var_Y );
	LAB_b = 200 * ( var_Y - var_Z );

}

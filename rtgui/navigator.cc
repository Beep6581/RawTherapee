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
#include <iccmatrices.h> // from rtengine

const double (*wprofiles[])[3]  = {xyz_sRGB, xyz_adobe, xyz_prophoto, xyz_widegamut, xyz_bruce, xyz_beta, xyz_best};
const double (*iwprofiles[])[3] = {sRGB_xyz, adobe_xyz, prophoto_xyz, widegamut_xyz, bruce_xyz, beta_xyz, best_xyz};
const char* wpnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB"};
const int numprofiles = 7;

#define D50x 0.96422
#define D50z 0.82521

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

//void Navigator::pointerMoved (bool validPos, int x, int y, int r, int g, int b) {
void Navigator::pointerMoved (bool validPos, Glib::ustring profile, int x, int y, int r, int g, int b) {

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
//		rgb2lab (r, g, b, LAB_l, LAB_a, LAB_b);
		rgb2lab (profile, r, g, b, LAB_l, LAB_a, LAB_b);
		
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

//void Navigator::rgb2lab (int r, int g, int b, int &LAB_l, int &LAB_a, int &LAB_b) {
void Navigator::rgb2lab (Glib::ustring profile, int r, int g, int b, int &LAB_l, int &LAB_a, int &LAB_b) {
	
	double xyz_rgb[3][3];
	double ep=216.0/24389.0;
	double ka=24389.0/27.0;

	volatile double var_R = r / 255.0;
	volatile double var_G = g / 255.0;
	volatile double var_B = b / 255.0;

//	if (profile=="sRGB") {//apply sRGB inverse gamma

// si on souhaite affichage = working space
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
/*	} 
// si on souhaite affichage = espace de sortie
	else 
	if (profile=="ProPhoto") {// apply inverse gamma 1.8
		var_R = pow ( var_R, 1.8);
		var_G = pow ( var_G, 1.8);
		var_B = pow ( var_B, 1.8);
	}
	else {// apply inverse gamma 2.2
		var_R = pow ( var_R, 2.2);
		var_G = pow ( var_G, 2.2);
		var_B = pow ( var_B, 2.2);
	}
*/
	for (int i=0; i<numprofiles; i++) {
		if (profile==wpnames[i]) {
			for (int m=0; m<3; m++) 
				for (int n=0; n<3; n++) {
					xyz_rgb[m][n] = wprofiles[i][m][n];
			}
			break;
		}
	}

	double varxx,varyy,varzz;
	double var_X = ( xyz_rgb[0][0]*var_R + xyz_rgb[0][1]*var_G + xyz_rgb[0][2]*var_B ) / D50x;
	double var_Y = ( xyz_rgb[1][0]*var_R + xyz_rgb[1][1]*var_G + xyz_rgb[1][2]*var_B ) ;
	double var_Z = ( xyz_rgb[2][0]*var_R + xyz_rgb[2][1]*var_G + xyz_rgb[2][2]*var_B ) / D50z;

	varxx = var_X>ep?pow (var_X, 1.0/3.0):( ka * var_X  +  16.0) / 116.0 ;
	varyy = var_Y>ep?pow (var_Y, 1.0/3.0):( ka * var_Y  +  16.0) / 116.0 ;
	varzz = var_Z>ep?pow (var_Z, 1.0/3.0):( ka * var_Z  +  16.0) / 116.0 ;
	LAB_l = ( 116 * varyy ) - 16;
	LAB_a = 500 * ( varxx - varyy );
	LAB_b = 200 * ( varyy - varzz );

}

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
#include <iomanip>
#include "navigator.h"
#include "toolpanel.h"
#include "../rtengine/iccmatrices.h"
#include "../rtengine/iccstore.h"
#include "../rtengine/curves.h"
#include "../rtengine/color.h"
#include "../rtengine/rt_math.h"
#include "options.h"

using namespace rtengine;

Navigator::Navigator ()
{

    set_label (M("MAIN_MSG_NAVIGATOR"));
    Gtk::VBox* mbox = Gtk::manage (new Gtk::VBox ());
    previewWindow = Gtk::manage (new PreviewWindow ());
    mbox->pack_start (*previewWindow, Gtk::PACK_SHRINK, 2);
    position = Gtk::manage (new Gtk::Label ());
    mbox->pack_start (*position, Gtk::PACK_SHRINK, 2);

    //labels
    lR = Gtk::manage (new Gtk::Label (M("NAVIGATOR_R")));
    lG = Gtk::manage (new Gtk::Label (M("NAVIGATOR_G")));
    lB = Gtk::manage (new Gtk::Label (M("NAVIGATOR_B")));
    lH = Gtk::manage (new Gtk::Label (M("NAVIGATOR_H")));
    lS = Gtk::manage (new Gtk::Label (M("NAVIGATOR_S")));
    lV = Gtk::manage (new Gtk::Label (M("NAVIGATOR_V")));
    lLAB_A = Gtk::manage (new Gtk::Label (M("NAVIGATOR_LAB_A")));
    lLAB_B = Gtk::manage (new Gtk::Label (M("NAVIGATOR_LAB_B")));
    lLAB_L = Gtk::manage (new Gtk::Label (M("NAVIGATOR_LAB_L")));

    // left-align labels
    lR->set_alignment(Gtk::ALIGN_START);
    lG->set_alignment(Gtk::ALIGN_START);
    lB->set_alignment(Gtk::ALIGN_START);
    lH->set_alignment(Gtk::ALIGN_START);
    lS->set_alignment(Gtk::ALIGN_START);
    lV->set_alignment(Gtk::ALIGN_START);
    lLAB_A->set_alignment(Gtk::ALIGN_START);
    lLAB_B->set_alignment(Gtk::ALIGN_START);
    lLAB_L->set_alignment(Gtk::ALIGN_START);

    //values
    R = Gtk::manage (new Gtk::Label ());
    G = Gtk::manage (new Gtk::Label ());
    B = Gtk::manage (new Gtk::Label ());
    H = Gtk::manage (new Gtk::Label ());
    S = Gtk::manage (new Gtk::Label ());
    V = Gtk::manage (new Gtk::Label ());
    LAB_A = Gtk::manage (new Gtk::Label ());
    LAB_B = Gtk::manage (new Gtk::Label ());
    LAB_L = Gtk::manage (new Gtk::Label ());

    // right-align values
    R->set_alignment(Gtk::ALIGN_END);
    G->set_alignment(Gtk::ALIGN_END);
    B->set_alignment(Gtk::ALIGN_END);
    H->set_alignment(Gtk::ALIGN_END);
    S->set_alignment(Gtk::ALIGN_END);
    V->set_alignment(Gtk::ALIGN_END);
    LAB_A->set_alignment(Gtk::ALIGN_END);
    LAB_B->set_alignment(Gtk::ALIGN_END);
    LAB_L->set_alignment(Gtk::ALIGN_END);

    // set font family and size
    Glib::ustring fontname;

#ifdef WIN32
    fontname = "Droid Sans Mono Slashed"; // font file is provided in the source tree in rtdata/fonts to be installed by the windows installer
#endif

#ifdef __linux__
    fontname = "Monospace";
#endif

#ifdef __APPLE__
    fontname = "Menlo";
#endif

    if (0) { // (fontname!=""){
        Glib::RefPtr<Gtk::CssProvider> cssProvider = Gtk::CssProvider::create();

        if (cssProvider) {
            cssProvider->load_from_data("Label { font-name: " + fontname + " }");
            R->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            G->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            B->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            H->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            S->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            V->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            LAB_A->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            LAB_B->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            LAB_L->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);

            lR->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lG->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lB->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lH->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lS->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lV->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lLAB_A->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lLAB_B->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
            lLAB_L->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);

            position->get_style_context()->add_provider(cssProvider, GTK_STYLE_PROVIDER_PRIORITY_USER);
        }
    }

    // setup the tables
    Gtk::Table* table0 = Gtk::manage (new Gtk::Table (1, 3)); //rows, cols The main table container
    // let's pack tables1,2-3 into table0


    // RGB
    Gtk::HBox* hbox1 = Gtk::manage (new Gtk::HBox ()); // container
    Gtk::Table* table1 = Gtk::manage (new Gtk::Table (3, 2));

    table1->attach (*lR, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table1->attach (*R,  1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table1->attach (*lG, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table1->attach (*G,  1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table1->attach (*lB, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table1->attach (*B,  1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    hbox1->pack_start (*table1, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox1->pack_start (*Gtk::manage (new  Gtk::VSeparator()), Gtk::PACK_SHRINK, 4);
    table0->attach (*hbox1, 0, 1, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    // HSV
    Gtk::HBox* hbox2 = Gtk::manage (new Gtk::HBox ()); // container
    Gtk::Table* table2 = Gtk::manage (new Gtk::Table (3, 2));

    table2->attach (*lH, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table2->attach (*H,  1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table2->attach (*lS, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table2->attach (*S,  1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table2->attach (*lV, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table2->attach (*V,  1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    hbox2->pack_start (*table2, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox2->pack_start (*Gtk::manage (new  Gtk::VSeparator()), Gtk::PACK_SHRINK, 4);
    table0->attach (*hbox2, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    // LAB
    Gtk::HBox* hbox3 = Gtk::manage (new Gtk::HBox ()); // container
    Gtk::Table* table3 = Gtk::manage (new Gtk::Table (3, 2));

    table3->attach (*lLAB_L, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table3->attach (*LAB_L,  1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table3->attach (*lLAB_A, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table3->attach (*LAB_A,  1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table3->attach (*lLAB_B, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 4, 0);
    table3->attach (*LAB_B,  1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    hbox3->pack_start (*table3, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox3->pack_start (*Gtk::manage (new  Gtk::HBox()), Gtk::PACK_SHRINK, 2);
    table0->attach (*hbox3, 2, 3, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 0, 0);

    table0->set_homogeneous(true); // all cells will be the same size as the largest cell.

    mbox->pack_start (*table0, Gtk::PACK_EXPAND_WIDGET, 2);
    add (*mbox);

    setInvalid ();
    show_all ();
}

void Navigator::setInvalid (int fullWidth, int fullHeight)
{
    if (fullWidth > 0 && fullHeight > 0) {
        position->set_text (Glib::ustring::compose (M("NAVIGATOR_XY_FULL"), fullWidth, fullHeight));
    } else {
        position->set_text (M("NAVIGATOR_XY_NA"));
    }

    R->set_text (M("NAVIGATOR_NA"));
    G->set_text (M("NAVIGATOR_NA"));
    B->set_text (M("NAVIGATOR_NA"));
    H->set_text (M("NAVIGATOR_NA"));
    S->set_text (M("NAVIGATOR_NA"));
    V->set_text (M("NAVIGATOR_NA"));
    LAB_A->set_text (M("NAVIGATOR_NA"));
    LAB_B->set_text (M("NAVIGATOR_NA"));
    LAB_L->set_text (M("NAVIGATOR_NA"));
}

// if !validPos then x/y contain the full image size
void Navigator::pointerMoved (bool validPos, Glib::ustring profile, Glib::ustring profileW, int x, int y, int r, int g, int b)
{

    if (!validPos) {
        setInvalid (x, y);
    } else {
        position->set_text (Glib::ustring::compose ("x: %1, y: %2", x, y));

        R->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), r * 100.f / 255.f) + Glib::ustring("%"));
        G->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), g * 100.f / 255.f) + Glib::ustring("%"));
        B->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), b * 100.f / 255.f) + Glib::ustring("%"));

        float h, s, v;
        Color::rgb2hsv (r * 0xffff / 0xff, g * 0xffff / 0xff, b * 0xffff / 0xff, h, s, v);
        H->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), h * 360.f) + Glib::ustring("\xc2\xb0"));
        S->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), s * 100.f) + Glib::ustring("%"));
        V->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), v * 100.f) + Glib::ustring("%"));

        float LAB_a, LAB_b, LAB_l;
        //rgb2lab (r, g, b, LAB_l, LAB_a, LAB_b);
        rgb2lab (profile, profileW, r, g, b, LAB_l, LAB_a, LAB_b);  // TODO: Really sure this function works?
        LAB_A->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), LAB_a));
        LAB_B->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), LAB_b));
        LAB_L->set_text (Glib::ustring::format(std::fixed, std::setprecision(1), LAB_l));
    }
}


void Navigator::rgb2lab (Glib::ustring profile, Glib::ustring profileW, int r, int g, int b, float &LAB_l, float &LAB_a, float &LAB_b)
{
    double xyz_rgb[3][3];
    const double ep = 216.0 / 24389.0;
    const double ka = 24389.0 / 27.0;

    double var_R = r / 255.0;
    double var_G = g / 255.0;
    double var_B = b / 255.0;

    Glib::ustring profileCalc;
    profileCalc = "sRGB"; //default

    if(options.rtSettings.HistogramWorking) {
        profileCalc = profileW;    //display working
    }

    else {// if you want display = output space
        if (profile == "RT_sRGB" || profile == "RT_sRGB_gBT709" || profile == "RT_sRGB_g10") {
            profileCalc = "sRGB";
        }

        if (profile == "ProPhoto" || profile == "RT_Large_gBT709" || profile == "RT_Large_g10"  || profile == "RT_Large_gsRGB") {
            profileCalc = "ProPhoto";
        }

        if (profile == "AdobeRGB1998" || profile == "RT_Medium_gsRGB") {
            profileCalc = "Adobe RGB";
        }

        if (profile == "WideGamutRGB") {
            profileCalc = "WideGamut";
        }
    }

    if(options.rtSettings.HistogramWorking) {//display working
        if (profileW == "sRGB") { //apply sRGB inverse gamma

            if ( var_R > 0.04045 ) {
                var_R = pow ( ( ( var_R + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_R = var_R / 12.92;
            }

            if ( var_G > 0.04045 ) {
                var_G = pow ( ( ( var_G + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_G = var_G / 12.92;
            }

            if ( var_B > 0.04045 ) {
                var_B = pow ( ( ( var_B + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_B = var_B / 12.92;
            }
        } else if (profileW == "ProPhoto") { // apply inverse gamma 1.8
            var_R = pow ( var_R, 1.8);
            var_G = pow ( var_G, 1.8);
            var_B = pow ( var_B, 1.8);
        } else { // apply inverse gamma 2.2
            var_R = pow ( var_R, 2.2);
            var_G = pow ( var_G, 2.2);
            var_B = pow ( var_B, 2.2);
        }
    } else { //display outout profile

        if (profile == "RT_sRGB" || profile == "RT_Large_gsRGB"  || profile == "RT_Medium_gsRGB") { //apply sRGB inverse gamma
            if ( var_R > 0.04045 ) {
                var_R = pow ( ( ( var_R + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_R = var_R / 12.92;
            }

            if ( var_G > 0.04045 ) {
                var_G = pow ( ( ( var_G + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_G = var_G / 12.92;
            }

            if ( var_B > 0.04045 ) {
                var_B = pow ( ( ( var_B + 0.055 ) / 1.055 ), rtengine::Color::sRGBGammaCurve);
            } else {
                var_B = var_B / 12.92;
            }
        }

        else if (profile == "RT_sRGB_gBT709"  || profile == "RT_Large_gBT709") { //
            if ( var_R > 0.0795 ) {
                var_R = pow ( ( ( var_R + 0.0954 ) / 1.0954 ), 2.2);
            } else {
                var_R = var_R / 4.5;
            }

            if ( var_G > 0.0795 ) {
                var_G = pow ( ( ( var_G + 0.0954 ) / 1.0954 ), 2.2);
            } else {
                var_G = var_G / 4.5;
            }

            if ( var_B > 0.0795 ) {
                var_B = pow ( ( ( var_B + 0.0954 ) / 1.0954 ), 2.2);
            } else {
                var_B = var_B / 4.5;
            }

        } else if (profile == "ProPhoto") { // apply inverse gamma 1.8

            var_R = pow ( var_R, 1.8);
            var_G = pow ( var_G, 1.8);
            var_B = pow ( var_B, 1.8);
        } else if (profile == "RT_sRGB_g10"  || profile == "RT_Large_g10") { // apply inverse gamma 1.8

            var_R = pow ( var_R, 1.);
            var_G = pow ( var_G, 1.);
            var_B = pow ( var_B, 1.);
        }

        else {// apply inverse gamma 2.2
            var_R = pow ( var_R, 2.2);
            var_G = pow ( var_G, 2.2);
            var_B = pow ( var_B, 2.2);
        }
    }

    // TMatrix wprof = rtengine::ICCStore::getInstance()->workingSpaceMatrix (profileW);

    TMatrix wprof = rtengine::ICCStore::getInstance()->workingSpaceMatrix (profileCalc);

    for (int m = 0; m < 3; m++)
        for (int n = 0; n < 3; n++) {
            xyz_rgb[m][n] = wprof[m][n];
        }

    double varxx, varyy, varzz;
    double var_X = ( xyz_rgb[0][0] * var_R + xyz_rgb[0][1] * var_G + xyz_rgb[0][2] * var_B ) / Color::D50x;
    double var_Y = ( xyz_rgb[1][0] * var_R + xyz_rgb[1][1] * var_G + xyz_rgb[1][2] * var_B ) ;
    double var_Z = ( xyz_rgb[2][0] * var_R + xyz_rgb[2][1] * var_G + xyz_rgb[2][2] * var_B ) / Color::D50z;

    varxx = var_X > ep ? cbrt(var_X) : ( ka * var_X  +  16.0) / 116.0 ;
    varyy = var_Y > ep ? cbrt(var_Y) : ( ka * var_Y  +  16.0) / 116.0 ;
    varzz = var_Z > ep ? cbrt(var_Z) : ( ka * var_Z  +  16.0) / 116.0 ;
    LAB_l = ( 116 * varyy ) - 16;
    LAB_a = 500 * ( varxx - varyy );
    LAB_b = 200 * ( varyy - varzz );

}

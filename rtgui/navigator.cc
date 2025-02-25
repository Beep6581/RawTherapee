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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <iomanip>
#include "navigator.h"
#include "previewwindow.h"
#include "toolpanel.h"
#include "rtengine/color.h"
#include "rtengine/improcfun.h"
#include "rtengine/rt_math.h"
#include "options.h"

using namespace rtengine;

namespace
{

const rtengine::procparams::ColorManagementParams DEFAULT_CMP;

}

Navigator::Navigator() :
    pointer_moved_delayed_call(50, 100),
    currentRGBUnit(options.navRGBUnit),
    currentHSVUnit(options.navHSVUnit)
{
    pointer_moved_delayed_call.setFunction(
        [this](bool validPos, const rtengine::procparams::ColorManagementParams *cmp, int x, int y, int r, int g, int b, bool isRaw)
        {
            if (!validPos) {
                setInvalid (x, y);
            } else {
                Glib::ustring s1, s2, s3;

                position->set_text (Glib::ustring::compose ("x: %1, y: %2", x, y));

                getRGBText (r, g, b, s1, s2, s3, isRaw);
                R->set_text (s1);
                G->set_text (s2);
                B->set_text (s3);
                if (isRaw) {
                    H->set_text ("--");
                    S->set_text ("--");
                    V->set_text ("--");
                    LAB_L->set_text ("--");
                    LAB_A->set_text ("--");
                    LAB_B->set_text ("--");
                } else {
                    float h, s, v;
                    float LAB_a, LAB_b, LAB_l;
                    Color::rgb2hsv01(r / 255.f, g / 255.f, b / 255.f, h, s, v);
                    getHSVText (h, s, v, s1, s2, s3);
                    H->set_text (s1);
                    S->set_text (s2);
                    V->set_text (s3);

                    ImProcFunctions::rgb2lab(
                        static_cast<std::uint8_t>(r),
                        static_cast<std::uint8_t>(g),
                        static_cast<std::uint8_t>(b),
                        LAB_l, LAB_a, LAB_b,
                        cmp != nullptr ? *cmp : DEFAULT_CMP,
                        true);
                    LAB_l /= 327.68f;
                    LAB_a /= 327.68f;
                    LAB_b /= 327.68f;
                    getLABText (LAB_l, LAB_a, LAB_b, s1, s2, s3);
                    LAB_L->set_text (s1);
                    LAB_A->set_text (s2);
                    LAB_B->set_text (s3);
                }
            }
        }
    );

    set_label (M("MAIN_MSG_NAVIGATOR"));
    set_name("Navigator");
    Gtk::Box* mbox = Gtk::manage (new Gtk::Box(Gtk::ORIENTATION_VERTICAL));
    previewWindow = Gtk::manage (new PreviewWindow ());
    mbox->pack_start (*previewWindow, Gtk::PACK_EXPAND_WIDGET, 2);
    dimension = Gtk::manage (new Gtk::Label ());
    mbox->pack_start (*dimension, Gtk::PACK_SHRINK, 2);
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
    
    // expand labels
    lR->set_hexpand();
    lG->set_hexpand();
    lB->set_hexpand();
    lH->set_hexpand();
    lS->set_hexpand();
    lV->set_hexpand();
    lLAB_A->set_hexpand();
    lLAB_B->set_hexpand();
    lLAB_L->set_hexpand();

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
    /*
    Glib::ustring fontname;

#ifdef _WIN32
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
    */

    // setup the tables
    Gtk::Grid* table0 = Gtk::manage (new Gtk::Grid()); //rows, cols The main table container
    // let's pack tables1,2-3 into table0


    // RGB
    Gtk::EventBox *evBox1 = Gtk::manage (new Gtk::EventBox());
    Gtk::Box* hbox1 = Gtk::manage (new Gtk::Box ());
    Gtk::Grid* table1 = Gtk::manage (new Gtk::Grid());
    
    table1->attach(*lR, 0, 0, 1, 1);
    table1->attach(*R, 1, 0, 1, 1);
    table1->attach(*lG, 0, 1, 1, 1);
    table1->attach(*G, 1, 1, 1, 1);
    table1->attach(*lB, 0, 2, 1, 1);
    table1->attach(*B, 1, 2, 1, 1);

    evBox1->add (*table1);
    evBox1->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &Navigator::cycleUnitsRGB));

    hbox1->pack_start (*evBox1, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox1->pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_VERTICAL)), Gtk::PACK_SHRINK, 4);
    table0->attach(*hbox1, 0, 0, 1, 1);

    // HSV
    Gtk::EventBox *evBox2 = Gtk::manage (new Gtk::EventBox());
    Gtk::Box* hbox2 = Gtk::manage (new Gtk::Box ());
    Gtk::Grid* table2 = Gtk::manage (new Gtk::Grid());

    table2->attach(*lH, 0, 0, 1, 1);
    table2->attach(*H, 1, 0, 1, 1);
    table2->attach(*lS, 0, 1, 1, 1);
    table2->attach(*S, 1, 1, 1, 1);
    table2->attach(*lV, 0, 2, 1, 1);
    table2->attach(*V, 1, 2, 1, 1);

    evBox2->add (*table2);
    evBox2->signal_button_release_event().connect_notify( sigc::mem_fun(*this, &Navigator::cycleUnitsHSV));

    hbox2->pack_start (*evBox2, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox2->pack_start (*Gtk::manage (new Gtk::Separator(Gtk::ORIENTATION_VERTICAL)), Gtk::PACK_SHRINK, 4);
    table0->attach(*hbox2, 1, 0, 1, 1);

    // LAB
    Gtk::Box* hbox3 = Gtk::manage (new Gtk::Box ());
    Gtk::Grid* table3 = Gtk::manage (new Gtk::Grid());

    table3->attach(*lLAB_L, 0, 0, 1, 1);
    table3->attach(*LAB_L, 1, 0, 1, 1);
    table3->attach(*lLAB_A, 0, 1, 1, 1);
    table3->attach(*LAB_A, 1, 1, 1, 1);
    table3->attach(*lLAB_B, 0, 2, 1, 1);
    table3->attach(*LAB_B, 1, 2, 1, 1);

    hbox3->pack_start (*table3, Gtk::PACK_EXPAND_WIDGET, 4);
    hbox3->pack_start (*Gtk::manage (new  Gtk::Box()), Gtk::PACK_SHRINK, 2);
    table0->attach(*hbox3, 2, 0, 1, 1);

    table0->set_column_homogeneous(true); // all cells will have equal width

    mbox->pack_start (*table0, Gtk::PACK_SHRINK, 2);
    add (*mbox);

    setInvalid ();
    show_all ();
}

Navigator::~Navigator()
{
    pointer_moved_delayed_call.cancel();
}

void Navigator::setInvalid (int fullWidth, int fullHeight)
{
    if (fullWidth > 0 && fullHeight > 0) {
        dimension->set_text (Glib::ustring::compose (M("NAVIGATOR_XY_FULL"), fullWidth, fullHeight));
    }
    position->set_text (M("NAVIGATOR_XY_NA"));

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

void Navigator::getRGBText (int r, int g, int b, Glib::ustring &sR, Glib::ustring &sG, Glib::ustring &sB, bool isRaw)
{
    if (isRaw) {
        sR = Glib::ustring::format(std::fixed, std::setprecision(0), r);
        sG = Glib::ustring::format(std::fixed, std::setprecision(0), g);
        sB = Glib::ustring::format(std::fixed, std::setprecision(0), b);
    } else {
        switch (currentRGBUnit) {
        case (Options::NavigatorUnit::R0_1):
            sR = Glib::ustring::format(std::fixed, std::setprecision(4), r / 255.f);
            sG = Glib::ustring::format(std::fixed, std::setprecision(4), g / 255.f);
            sB = Glib::ustring::format(std::fixed, std::setprecision(4), b / 255.f);
            break;
        case (Options::NavigatorUnit::R0_255):
            sR = Glib::ustring::format(std::fixed, std::setprecision(0), r);
            sG = Glib::ustring::format(std::fixed, std::setprecision(0), g);
            sB = Glib::ustring::format(std::fixed, std::setprecision(0), b);
            break;
        case (Options::NavigatorUnit::PERCENT):
        default:
            sR = Glib::ustring::format(std::fixed, std::setprecision(1), r * 100.f / 255.f) + Glib::ustring("%");
            sG = Glib::ustring::format(std::fixed, std::setprecision(1), g * 100.f / 255.f) + Glib::ustring("%");
            sB = Glib::ustring::format(std::fixed, std::setprecision(1), b * 100.f / 255.f) + Glib::ustring("%");
        }
    }
}

void Navigator::getHSVText (float h, float s, float v, Glib::ustring &sH, Glib::ustring &sS, Glib::ustring &sV)
{
    switch (currentHSVUnit) {
    case (Options::NavigatorUnit::R0_1):
        sH = Glib::ustring::format(std::fixed, std::setprecision(4), h);
        sS = Glib::ustring::format(std::fixed, std::setprecision(4), s);
        sV = Glib::ustring::format(std::fixed, std::setprecision(4), v);
        break;
    case (Options::NavigatorUnit::R0_255):
        sH = Glib::ustring::format(std::fixed, std::setprecision(0), h * 255);
        sS = Glib::ustring::format(std::fixed, std::setprecision(0), s * 255);
        sV = Glib::ustring::format(std::fixed, std::setprecision(0), v * 255);
        break;
    case (Options::NavigatorUnit::PERCENT):
    default:
        sH = Glib::ustring::format(std::fixed, std::setprecision(1), h * 360.f) + Glib::ustring("\xc2\xb0");
        sS = Glib::ustring::format(std::fixed, std::setprecision(1), s * 100.f) + Glib::ustring("%");
        sV = Glib::ustring::format(std::fixed, std::setprecision(1), v * 100.f) + Glib::ustring("%");
    }
}

void Navigator::getLABText (float l, float a, float b, Glib::ustring &sL, Glib::ustring &sA, Glib::ustring &sB)
{
    sL = Glib::ustring::format(std::fixed, std::setprecision(1), l);
    sA = Glib::ustring::format(std::fixed, std::setprecision(1), a);
    sB = Glib::ustring::format(std::fixed, std::setprecision(1), b);
}

// if !validPos then x/y contain the full image size
void Navigator::pointerMoved (bool validPos, const rtengine::procparams::ColorManagementParams &cmp, int x, int y, int r, int g, int b, bool isRaw)
{
    pointer_moved_delayed_call(validPos, &cmp, x, y, r, g, b, isRaw);
}

void Navigator::cycleUnitsRGB (GdkEventButton *event) {
    uint16_t v = (uint16_t)currentRGBUnit;
    ++v;
    if (v == (uint16_t)Options::NavigatorUnit::_COUNT) {
        v = 0;
    }
    options.navRGBUnit = currentRGBUnit = (Options::NavigatorUnit)v;

    switch (currentRGBUnit) {
    case Options::NavigatorUnit::R0_1:
        R->set_text ("[0-1]");
        G->set_text ("[0-1]");
        B->set_text ("[0-1]");
        break;
    case Options::NavigatorUnit::R0_255:
        R->set_text ("[0-255]");
        G->set_text ("[0-255]");
        B->set_text ("[0-255]");
        break;
    case Options::NavigatorUnit::PERCENT:
    default:
        R->set_text ("[%]");
        G->set_text ("[%]");
        B->set_text ("[%]");
        break;
    }
    sig_cycle_rgb.emit();
}

void Navigator::cycleUnitsHSV (GdkEventButton *event) {
    uint16_t v = (uint16_t)currentHSVUnit;
    ++v;
    if (v == (uint16_t)Options::NavigatorUnit::_COUNT) {
        v = 0;
    }
    options.navHSVUnit = currentHSVUnit = (Options::NavigatorUnit)v;

    switch (currentHSVUnit) {
    case Options::NavigatorUnit::R0_1:
        H->set_text ("[0-1]");
        S->set_text ("[0-1]");
        V->set_text ("[0-1]");
        break;
    case Options::NavigatorUnit::R0_255:
        H->set_text ("[0-255]");
        S->set_text ("[0-255]");
        V->set_text ("[0-255]");
        break;
    case Options::NavigatorUnit::PERCENT:
    default:
        H->set_text ("[\xc2\xb0]");
        S->set_text ("[%]");
        V->set_text ("[%]");
        break;
    }
    sig_cycle_hsv.emit();
}

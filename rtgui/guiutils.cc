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
#include "../rtengine/rt_math.h"

#include "guiutils.h"
#include "options.h"
#include "../rtengine/rt_math.h"
#include "../rtengine/utils.h"
#include "../rtengine/safegtk.h"
#include "rtimage.h"
#include "multilangmgr.h"

#include <assert.h>

using namespace std;

#if TRACE_MYRWMUTEX==1 && !defined NDEBUG
unsigned int MyReaderLock::readerLockCounter = 0;
unsigned int MyWriterLock::writerLockCounter = 0;
#endif

Glib::RefPtr<Gdk::Pixbuf> MyExpander::inconsistentPBuf;
Glib::RefPtr<Gdk::Pixbuf> MyExpander::enabledPBuf;
Glib::RefPtr<Gdk::Pixbuf> MyExpander::disabledPBuf;
Glib::RefPtr<Gdk::Pixbuf> MyExpander::openedPBuf;
Glib::RefPtr<Gdk::Pixbuf> MyExpander::closedPBuf;

Glib::ustring escapeHtmlChars(const Glib::ustring &src)
{

    // Sources chars to be escaped
    static const Glib::ustring srcChar("&<>");

    // Destination strings, in the same order than the source
    static std::vector<Glib::ustring> dstChar(3);
    dstChar.at(0) = "&amp;";
    dstChar.at(1) = "&lt;";
    dstChar.at(2) = "&gt;";

    // Copying the original string, that will be modified
    Glib::ustring dst(src);

    // Iterating all chars of the copy of the source string
    for (size_t i = 0; i < dst.length();) {

        // Looking out if it's part of the characters to be escaped
        size_t pos = srcChar.find_first_of(dst.at(i), 0);

        if (pos != Glib::ustring::npos) {
            // If yes, replacing the char in the destination string
            dst.replace(i, 1, dstChar.at(pos));
            // ... and going forward  by the length of the new string
            i += dstChar.at(pos).length();
        } else {
            ++i;
        }
    }

    return dst;
}

bool removeIfThere (Gtk::Container* cont, Gtk::Widget* w, bool increference)
{

    Glib::ListHandle<Gtk::Widget*> list = cont->get_children ();
    Glib::ListHandle<Gtk::Widget*>::iterator i = list.begin ();

    for (; i != list.end() && *i != w; ++i);

    if (i != list.end()) {
        if (increference) {
            w->reference ();
        }

        cont->remove (*w);
        return true;
    } else {
        return false;
    }
}

void thumbInterp (const unsigned char* src, int sw, int sh, unsigned char* dst, int dw, int dh)
{

    if (options.thumbInterp == 0) {
        rtengine::nearestInterp (src, sw, sh, dst, dw, dh);
    } else if (options.thumbInterp == 1) {
        rtengine::bilinearInterp (src, sw, sh, dst, dw, dh);
    }
}

Glib::ustring removeExtension (const Glib::ustring& filename)
{

    Glib::ustring bname = Glib::path_get_basename(filename);
    size_t lastdot = bname.find_last_of ('.');
    size_t lastwhitespace = bname.find_last_of (" \t\f\v\n\r");

    if (lastdot != bname.npos && (lastwhitespace == bname.npos || lastdot > lastwhitespace)) {
        return filename.substr (0, filename.size() - (bname.size() - lastdot));
    } else {
        return filename;
    }
}

Glib::ustring getExtension (const Glib::ustring& filename)
{

    Glib::ustring bname = Glib::path_get_basename(filename);
    size_t lastdot = bname.find_last_of ('.');
    size_t lastwhitespace = bname.find_last_of (" \t\f\v\n\r");

    if (lastdot != bname.npos && (lastwhitespace == bname.npos || lastdot > lastwhitespace)) {
        return filename.substr (filename.size() - (bname.size() - lastdot) + 1, filename.npos);
    } else {
        return "";
    }
}

bool confirmOverwrite (Gtk::Window& parent, const std::string& filename)
{
    bool safe = true;

    if (safe_file_test (filename, Glib::FILE_TEST_EXISTS)) {
        Glib::ustring msg_ = Glib::ustring ("<b>\"") + Glib::path_get_basename (filename) + "\": "
                             + M("MAIN_MSG_ALREADYEXISTS") + "</b>\n" + M("MAIN_MSG_QOVERWRITE");
        Gtk::MessageDialog msgd (parent, msg_, true, Gtk::MESSAGE_WARNING, Gtk::BUTTONS_YES_NO, true);
        safe = (msgd.run () == Gtk::RESPONSE_YES);
    }

    return safe;
}

void writeFailed (Gtk::Window& parent, const std::string& filename)
{
    Glib::ustring msg_ = Glib::ustring::compose(M("MAIN_MSG_WRITEFAILED"), filename);
    Gtk::MessageDialog msgd (parent, msg_, true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
    msgd.run ();
}

void drawCrop (Cairo::RefPtr<Cairo::Context> cr, int imx, int imy, int imw, int imh, int startx, int starty, double scale, const rtengine::procparams::CropParams& cparams, bool drawGuide, bool useBgColor, bool fullImageVisible)
{

    cr->set_line_width (0.);
    cr->rectangle (imx, imy, imw, imh);
    cr->clip ();

    double c1x = (cparams.x - startx) * scale;
    double c1y = (cparams.y - starty) * scale;
    double c2x = (cparams.x + cparams.w - startx) * scale - (fullImageVisible ? 0.0 : 1.0);
    double c2y = (cparams.y + cparams.h - starty) * scale - (fullImageVisible ? 0.0 : 1.0);

    // crop overlay color, linked with crop windows background
    if (options.bgcolor == 0 || !useBgColor) {
        cr->set_source_rgba (options.cutOverlayBrush[0], options.cutOverlayBrush[1], options.cutOverlayBrush[2], options.cutOverlayBrush[3]);
    } else if (options.bgcolor == 1) {
        cr->set_source_rgb (0, 0, 0);
    } else if (options.bgcolor == 2) {
        cr->set_source_rgb (1, 1, 1);
    }


    cr->rectangle (imx, imy, imw + 0.5, round(c1y) + 0.5);
    cr->rectangle (imx, round(imy + c2y) + 0.5, imw + 0.5, round(imh - c2y) + 0.5);
    cr->rectangle (imx, round(imy + c1y) + 0.5, round(c1x) + 0.5, round(c2y - c1y + 1) + 0.5);
    cr->rectangle (round(imx + c2x) + 0.5, round(imy + c1y) + 0.5, round(imw - c2x) + 0.5, round(c2y - c1y + 1) + 0.5);
    cr->fill ();

    // rectangle around the cropped area and guides
    if (cparams.guide != "None" && drawGuide) {
        double rectx1 = round(c1x) + imx + 0.5;
        double recty1 = round(c1y) + imy + 0.5;
        double rectx2 = round(c2x) + imx + 0.5;
        double recty2 = round(c2y) + imy + 0.5;

        if(fullImageVisible) {
            rectx2 = min(rectx2, imx + imw - 0.5);
            recty2 = min(recty2, imy + imh - 0.5);
        }

        cr->set_line_width (1.0);
        cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
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

        if (cparams.guide != "Rule of diagonals" && cparams.guide != "Golden Triangle 1" && cparams.guide != "Golden Triangle 2") {
            // draw guide lines
            std::vector<double> horiz_ratios;
            std::vector<double> vert_ratios;

            if (cparams.guide == "Rule of thirds") {
                horiz_ratios.push_back (1.0 / 3.0);
                horiz_ratios.push_back (2.0 / 3.0);
                vert_ratios.push_back (1.0 / 3.0);
                vert_ratios.push_back (2.0 / 3.0);
            } else if (!strncmp(cparams.guide.data(), "Harmonic means", 14)) {
                horiz_ratios.push_back (1.0 - 0.618);
                horiz_ratios.push_back (0.618);
                vert_ratios.push_back (0.618);
                vert_ratios.push_back (1.0 - 0.618);
            } else if (cparams.guide == "Grid") {
                // To have even distribution, normalize it a bit
                const int longSideNumLines = 10;

                int w = rectx2 - rectx1, h = recty2 - recty1;

                if (w > longSideNumLines && h > longSideNumLines) {
                    if (w > h) {
                        for (int i = 1; i < longSideNumLines; i++) {
                            vert_ratios.push_back ((double)i / longSideNumLines);
                        }

                        int shortSideNumLines = (int)round(h * (double)longSideNumLines / w);

                        for (int i = 1; i < shortSideNumLines; i++) {
                            horiz_ratios.push_back ((double)i / shortSideNumLines);
                        }
                    } else {
                        for (int i = 1; i < longSideNumLines; i++) {
                            horiz_ratios.push_back ((double)i / longSideNumLines);
                        }

                        int shortSideNumLines = (int)round(w * (double)longSideNumLines / h);

                        for (int i = 1; i < shortSideNumLines; i++) {
                            vert_ratios.push_back ((double)i / shortSideNumLines);
                        }
                    }
                }
            } else if (cparams.guide == "ePassport") {
                /* Official measurements do not specify exact ratios, just min/max measurements within which the eyes and chin-crown distance must lie. I averaged those measurements to produce these guides.
                 * The first horizontal guide is for the crown, the second is rougly for the nostrils, the third is for the chin.
                 * http://www.homeoffice.gov.uk/agencies-public-bodies/ips/passports/information-photographers/
                 * "(...) the measurement of the face from the bottom of the chin to the crown (ie the top of the head, not the top of the hair) is between 29mm and 34mm."
                 */
                horiz_ratios.push_back (7.0 / 45.0);
                horiz_ratios.push_back (26.0 / 45.0);
                horiz_ratios.push_back (37.0 / 45.0);
                vert_ratios.push_back (0.5);
            }

            // Horizontals
            for (size_t i = 0; i < horiz_ratios.size(); i++) {
                cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
                cr->move_to (rectx1, recty1 + round((recty2 - recty1) * horiz_ratios[i]));
                cr->line_to (rectx2, recty1 + round((recty2 - recty1) * horiz_ratios[i]));
                cr->stroke ();
                cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (rectx1, recty1 + round((recty2 - recty1) * horiz_ratios[i]));
                cr->line_to (rectx2, recty1 + round((recty2 - recty1) * horiz_ratios[i]));
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }

            // Verticals
            for (size_t i = 0; i < vert_ratios.size(); i++) {
                cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
                cr->move_to (rectx1 + round((rectx2 - rectx1) * vert_ratios[i]), recty1);
                cr->line_to (rectx1 + round((rectx2 - rectx1) * vert_ratios[i]), recty2);
                cr->stroke ();
                cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (rectx1 + round((rectx2 - rectx1) * vert_ratios[i]), recty1);
                cr->line_to (rectx1 + round((rectx2 - rectx1) * vert_ratios[i]), recty2);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }
        } else if (cparams.guide == "Rule of diagonals") {
            double corners_from[4][2];
            double corners_to[4][2];
            int mindim = min(rectx2 - rectx1, recty2 - recty1);
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

            for (int i = 0; i < 4; i++) {
                cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_dash (ds, 0);
                cr->move_to (corners_from[i][0], corners_from[i][1]);
                cr->line_to (corners_to[i][0], corners_to[i][1]);
                cr->stroke ();
                ds.resize (0);
                cr->set_dash (ds, 0);
            }
        } else if (cparams.guide == "Golden Triangle 1" || cparams.guide == "Golden Triangle 2") {
            // main diagonal
            if(cparams.guide == "Golden Triangle 2") {
                std::swap(rectx1, rectx2);
            }

            cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
            cr->move_to (rectx1, recty1);
            cr->line_to (rectx2, recty2);
            cr->stroke ();
            cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
            std::valarray<double> ds (1);
            ds[0] = 4;
            cr->set_dash (ds, 0);
            cr->move_to (rectx1, recty1);
            cr->line_to (rectx2, recty2);
            cr->stroke ();
            ds.resize (0);
            cr->set_dash (ds, 0);

            double height = recty2 - recty1;
            double width = rectx2 - rectx1;
            double d = sqrt(height * height + width * width);
            double alpha = asin(width / d);
            double beta = asin(height / d);
            double a = sin(beta) * height;
            double b = sin(alpha) * height;

            double x = (a * b) / height;
            double y = height - (b * (d - a)) / width;
            cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
            cr->move_to (rectx1, recty2);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
            ds.resize (1);
            ds[0] = 4;
            cr->set_dash (ds, 0);
            cr->move_to (rectx1, recty2);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            ds.resize (0);
            cr->set_dash (ds, 0);

            x = width - (a * b) / height;
            y = (b * (d - a)) / width;
            cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
            cr->move_to (rectx2, recty1);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
            ds.resize (1);
            ds[0] = 4;
            cr->set_dash (ds, 0);
            cr->move_to (rectx2, recty1);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            ds.resize (0);
            cr->set_dash (ds, 0);
        }
    }

    cr->reset_clip ();
}

bool ExpanderBox::on_expose_event(GdkEventExpose* event)
{
    bool retVal = Gtk::EventBox::on_expose_event(event);

    if (!options.useSystemTheme) {
        Glib::RefPtr<Gdk::Window> window = get_window();
        Glib::RefPtr<Gtk::Style> style = get_style ();
        Cairo::RefPtr<Cairo::Context> cr = window->create_cairo_context();

        int x_, y_, w_, h_, foo;
        window->get_geometry(x_, y_, w_, h_, foo);
        double x = 0.;
        double y = 0.;
        double w = double(w_);
        double h = double(h_);

        cr->set_antialias (Cairo::ANTIALIAS_NONE);

        // draw a frame
        cr->set_line_width (1.0);
        Gdk::Color c = style->get_fg (Gtk::STATE_NORMAL);
        cr->set_source_rgb (c.get_red_p(), c.get_green_p(), c.get_blue_p());
        cr->move_to(x + 0.5, y + 0.5);
        cr->line_to(x + w, y + 0.5);
        cr->line_to(x + w, y + h);
        cr->line_to(x + 0.5, y + h);
        cr->line_to(x + 0.5, y + 0.5);
        cr->stroke ();
    }

    return retVal;
}

ExpanderBox::ExpanderBox( Gtk::Container *p): pC(p)
{
    set_name ("ExpanderBox");
    updateStyle();
}

void ExpanderBox::on_style_changed (const Glib::RefPtr<Gtk::Style>& style)
{
    updateStyle();
}

void ExpanderBox::updateStyle()
{
    set_border_width(options.slimUI ? 2 : 8);  // Outer space around the tool's frame 2:7
}

void ExpanderBox::show_all()
{
    // ask childs to show themselves, but not us (remain unchanged)
    Gtk::Container::show_all_children(true);
}

void ExpanderBox::showBox()
{
    Gtk::EventBox::show();
}

void ExpanderBox::hideBox()
{
    Gtk::EventBox::hide();
}

void MyExpander::init()
{
    inconsistentPBuf = Gdk::Pixbuf::create_from_file(RTImage::findIconAbsolutePath("expanderInconsistent.png"));
    enabledPBuf = Gdk::Pixbuf::create_from_file(RTImage::findIconAbsolutePath("expanderEnabled.png"));
    disabledPBuf = Gdk::Pixbuf::create_from_file(RTImage::findIconAbsolutePath("expanderDisabled.png"));
    openedPBuf = Gdk::Pixbuf::create_from_file(RTImage::findIconAbsolutePath("expanderOpened.png"));
    closedPBuf = Gdk::Pixbuf::create_from_file(RTImage::findIconAbsolutePath("expanderClosed.png"));
}

MyExpander::MyExpander(bool useEnabled, Gtk::Widget* titleWidget) :
    enabled(false), inconsistent(false), flushEvent(false), expBox(NULL),
    child(NULL), headerWidget(NULL), statusImage(NULL),
    label(NULL), useEnabled(useEnabled)
{
    set_spacing(options.slimUI ? 0 : 2);
    set_name("MyExpander");
    set_can_focus(false);

    headerHBox = Gtk::manage( new Gtk::HBox());
    headerHBox->set_can_focus(false);

    if (useEnabled) {
        statusImage = Gtk::manage(new Gtk::Image(disabledPBuf));
        imageEvBox = Gtk::manage(new Gtk::EventBox());
        imageEvBox->add(*statusImage);
        imageEvBox->set_above_child(true);
        imageEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_enabled_change) );
        imageEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        imageEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        headerHBox->pack_start(*imageEvBox, Gtk::PACK_SHRINK, 0);
    } else {
        statusImage = Gtk::manage(new Gtk::Image(openedPBuf));
        headerHBox->pack_start(*statusImage, Gtk::PACK_SHRINK, 0);
    }

    statusImage->set_can_focus(false);

    if (titleWidget) {
        headerHBox->pack_start(*titleWidget, Gtk::PACK_EXPAND_WIDGET, 0);
        headerWidget = titleWidget;
    }

    titleEvBox = Gtk::manage(new Gtk::EventBox());
    titleEvBox->set_name("MyExpanderTitle");
    titleEvBox->add(*headerHBox);
    titleEvBox->set_above_child(false);  // this is the key! By making it below the child, they will get the events first.
    titleEvBox->set_can_focus(false);

    pack_start(*titleEvBox, Gtk::PACK_EXPAND_WIDGET, 0);

    updateStyle();
    titleEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_toggle) );
    titleEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_title), false);
    titleEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_title), false);
}

MyExpander::MyExpander(bool useEnabled, Glib::ustring titleLabel) :
    enabled(false), inconsistent(false), flushEvent(false), expBox(NULL),
    child(NULL), headerWidget(NULL), statusImage(NULL),
    label(NULL), useEnabled(useEnabled)
{
    set_spacing(options.slimUI ? 0 : 2);
    set_name("MyExpander");
    set_can_focus(false);

    headerHBox = Gtk::manage( new Gtk::HBox());
    headerHBox->set_can_focus(false);


    if (useEnabled) {
        statusImage = Gtk::manage(new Gtk::Image(disabledPBuf));
        imageEvBox = Gtk::manage(new Gtk::EventBox());
        imageEvBox->add(*statusImage);
        imageEvBox->set_above_child(true);
        imageEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_enabled_change) );
        imageEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        imageEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        headerHBox->pack_start(*imageEvBox, Gtk::PACK_SHRINK, 0);
    } else {
        statusImage = Gtk::manage(new Gtk::Image(openedPBuf));
        headerHBox->pack_start(*statusImage, Gtk::PACK_SHRINK, 0);
    }

    statusImage->set_can_focus(false);

    Glib::ustring str("-");

    if (!titleLabel.empty()) {
        str = titleLabel;
    }

    label = Gtk::manage(new Gtk::Label());
    label->set_alignment(Gtk::ALIGN_LEFT, Gtk::ALIGN_CENTER);
    label->set_markup(Glib::ustring("<b>") + escapeHtmlChars(titleLabel) + Glib::ustring("</b>"));
    headerHBox->pack_start(*label, Gtk::PACK_EXPAND_WIDGET, 0);

    titleEvBox = Gtk::manage(new Gtk::EventBox());
    titleEvBox->set_name("MyExpanderTitle");
    titleEvBox->add(*headerHBox);
    titleEvBox->set_above_child(false);  // this is the key! By make it below the child, they will get the events first.
    titleEvBox->set_can_focus(false);

    pack_start(*titleEvBox, Gtk::PACK_EXPAND_WIDGET, 0);

    updateStyle();
    titleEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_toggle));
    titleEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_title), false);
    titleEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_title), false);
}

bool MyExpander::on_enter_leave_title (GdkEventCrossing* event)
{
    if (is_sensitive()) {
        if (event->type == GDK_ENTER_NOTIFY) {
            titleEvBox->set_state(Gtk::STATE_PRELIGHT);
            queue_draw();
        } else if (event->type == GDK_LEAVE_NOTIFY) {
            titleEvBox->set_state(Gtk::STATE_NORMAL);
            queue_draw();
        }
    }

    return true;
}

bool MyExpander::on_enter_leave_enable (GdkEventCrossing* event)
{
    if (is_sensitive()) {
        if (event->type == GDK_ENTER_NOTIFY) {
            imageEvBox->set_state(Gtk::STATE_PRELIGHT);
            queue_draw();
        } else if (event->type == GDK_LEAVE_NOTIFY) {
            imageEvBox->set_state(Gtk::STATE_NORMAL);
            queue_draw();
        }
    }

    return true;
}

void MyExpander::updateStyle()
{
    headerHBox->set_spacing(options.slimUI ? 2 : 5);
    headerHBox->set_border_width(options.slimUI ? 1 : 2);
    set_spacing(0);
    set_border_width(options.slimUI ? 0 : 1);

    if (expBox) {
        expBox->updateStyle();
    }
}

void MyExpander::setLabel (Glib::ustring newLabel)
{
    if (label) {
        label->set_markup(Glib::ustring("<b>") + escapeHtmlChars(newLabel) + Glib::ustring("</b>"));
    }
}

void MyExpander::setLabel (Gtk::Widget *newWidget)
{
    if (headerWidget) {
        removeIfThere(headerHBox, headerWidget, false);
        headerHBox->pack_start(*newWidget, Gtk::PACK_EXPAND_WIDGET, 0);
    }
}

bool MyExpander::get_inconsistent()
{
    return inconsistent;
}

void MyExpander::set_inconsistent(bool isInconsistent)
{
    if (inconsistent != isInconsistent) {
        inconsistent = isInconsistent;

        if (useEnabled) {
            if (isInconsistent) {
                statusImage->set(inconsistentPBuf);
            } else {
                if (enabled) {
                    statusImage->set(enabledPBuf);
                } else {
                    statusImage->set(disabledPBuf);
                }
            }
        }

    }
}

bool MyExpander::getUseEnabled()
{
    return useEnabled;
}

bool MyExpander::getEnabled()
{
    return enabled;
}

void MyExpander::setEnabled(bool isEnabled)
{
    if (isEnabled != enabled) {
        if (useEnabled) {
            if (enabled) {
                enabled = false;

                if (!inconsistent) {
                    statusImage->set(disabledPBuf);
                    message.emit();
                }
            } else {
                enabled = true;

                if (!inconsistent) {
                    statusImage->set(enabledPBuf);
                    message.emit();
                }
            }
        }
    }
}

void MyExpander::setEnabledTooltipMarkup(Glib::ustring tooltipMarkup)
{
    if (useEnabled) {
        statusImage->set_tooltip_markup(tooltipMarkup);
    }
}

void MyExpander::setEnabledTooltipText(Glib::ustring tooltipText)
{
    if (useEnabled) {
        statusImage->set_tooltip_text(tooltipText);
    }
}

void MyExpander::set_expanded( bool expanded )
{
    if (!expBox) {
        return;
    }

    bool isVisible = expBox->is_visible();

    if (isVisible == expanded) {
        return;
    }

    if (!useEnabled) {
        if (expanded ) {
            statusImage->set(openedPBuf);
        } else {
            statusImage->set(closedPBuf);
        }
    }

    if (expanded) {
        expBox->showBox();
    } else {
        expBox->hideBox();
    }
}

bool MyExpander::get_expanded()
{
    return expBox ? expBox->get_visible() : false;
}

void MyExpander::add  (Gtk::Container& widget)
{
    child = &widget;
    expBox = Gtk::manage (new ExpanderBox (child));
    expBox->add (*child);
    pack_start(*expBox, Gtk::PACK_SHRINK, 0);
    child->show();
    expBox->hideBox();
}

bool MyExpander::on_toggle(GdkEventButton* event)
{
    if (flushEvent) {
        flushEvent = false;
        return false;
    }

    if (!expBox || event->button != 1) {
        return false;
    }

    bool isVisible = expBox->is_visible();

    if (!useEnabled) {
        if (isVisible) {
            statusImage->set(closedPBuf);
        } else {
            statusImage->set(openedPBuf);
        }
    }

    if (isVisible) {
        expBox->hideBox();
    } else {
        expBox->showBox();
    }

    return false;
}

Gtk::Container* MyExpander::getChild()
{
    return child;
}

// used to connect a function to the enabled_toggled signal
MyExpander::type_signal_enabled_toggled MyExpander::signal_enabled_toggled()
{
    return message;
}

// internal use ; when the user clicks on the toggle button, it calls this method that will emit an enabled_change event
bool MyExpander::on_enabled_change(GdkEventButton* event)
{
    if (event->button == 1) {
        if (enabled) {
            enabled = false;
            statusImage->set(disabledPBuf);
        } else {
            enabled = true;
            statusImage->set(enabledPBuf);
        }

        message.emit();
        flushEvent = true;
    }

    return false;
}

/*
 *
 * Derived class of some widgets to properly handle the scroll wheel ;
 * the user has to use the Shift key to be able to change the widget's value,
 * otherwise the mouse wheel will scroll the editor's tabs content.
 *
 */
MyScrolledWindow::MyScrolledWindow ()
{
    set_size_request(-1, 30);
}

bool MyScrolledWindow::on_scroll_event (GdkEventScroll* event)
{
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
            value2 = value + step;

            if (value2 > upper) {
                value2 = upper;
            }

            if (value2 != value) {
                scroll->set_value(value2);
            }
        } else {
            value2 = value - step;

            if (value2 < lower) {
                value2 = lower;
            }

            if (value2 != value) {
                scroll->set_value(value2);
            }
        }
    }

    return true;
}

MyComboBoxText::MyComboBoxText ()
{
    set_size_request(40, -1);
}

bool MyComboBoxText::on_scroll_event (GdkEventScroll* event)
{

    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::ComboBoxText::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

MyComboBox::MyComboBox ()
{
    set_size_request(40, -1);
}

bool MyComboBox::on_scroll_event (GdkEventScroll* event)
{

    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::ComboBox::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

MySpinButton::MySpinButton ()
{
    Gtk::Border border;
    border.bottom = 0;
    border.top = 0;
    border.left = 3;
    border.right = 3;
    set_inner_border(border);
    set_numeric(true);
    set_wrap(false);
    set_alignment(Gtk::ALIGN_RIGHT);
}

void MySpinButton::updateSize()
{
    double vMin, vMax;
    int maxAbs;
    unsigned int digits, digits2;
    unsigned int maxLen;

    get_range(vMin, vMax);

    digits = get_digits();
    maxAbs = (int)(fmax(fabs(vMin), fabs(vMax)) + 0.000001);

    if (maxAbs == 0) {
        digits2 = 1;
    } else {
        digits2 = (int)(log10(double(maxAbs)) + 0.000001);
        digits2++;
    }

    maxLen = digits + digits2 + (vMin < 0 ? 1 : 0) + (digits > 0 ? 1 : 0);
    set_max_length(maxLen);
    set_width_chars(maxLen);
}

bool MySpinButton::on_key_press_event (GdkEventKey* event)
{
    double vMin, vMax;
    get_range(vMin, vMax);

    if ( (event->string[0] >= 'a' && event->string[0] <= 'z')
            || (event->string[0] >= 'A' && event->string[0] <= 'Z')
            || event->string[0] == '+' || (event->string[0] == '-' && vMin >= 0)
            || event->string[0] == '=' || event->string[0] == '_'
       ) {
        return false;
    } else {
        if(event->string[0] == ',') {
            event->keyval = GDK_period;
            event->string[0] = '.';
        }

        return Gtk::Widget::on_key_press_event(event);
    }
}

bool MySpinButton::on_scroll_event (GdkEventScroll* event)
{
    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::SpinButton::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

bool MyHScale::on_scroll_event (GdkEventScroll* event)
{

    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::HScale::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

bool MyHScale::on_key_press_event (GdkEventKey* event)
{

    if ( event->string[0] == '+' || event->string[0] == '-' ) {
        return false;
    } else {
        return Gtk::Widget::on_key_press_event(event);
    }
}

MyFileChooserButton::MyFileChooserButton (const Glib::ustring& title, Gtk::FileChooserAction action) : Gtk::FileChooserButton(title, action)
{
    set_size_request(20, -1);
}

// For an unknown reason (a bug ?), it doesn't work when action = FILE_CHOOSER_ACTION_SELECT_FOLDER !
bool MyFileChooserButton::on_scroll_event (GdkEventScroll* event)
{

    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::FileChooserButton::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

void bindCurrentFolder (Gtk::FileChooser& chooser, Glib::ustring& variable)
{
    chooser.signal_selection_changed ().connect ([&]()
    {
        const auto current_folder = chooser.get_current_folder ();

        if (!current_folder.empty ())
            variable = current_folder;
    });

    if (!variable.empty ())
        chooser.set_current_folder (variable);
}

TextOrIcon::TextOrIcon (Glib::ustring fname, Glib::ustring labelTx, Glib::ustring tooltipTx, TOITypes type)
{

    imgIcon = 0;
    label = 0;
    filename = fname;
    labelText = labelTx;
    tooltipText = tooltipTx;

    switchTo(type);
}

TextOrIcon::~TextOrIcon ()
{
    if (imgIcon) {
        delete imgIcon;
    }

    if (label) {
        delete label;
    }
}

void TextOrIcon::switchTo(TOITypes type)
{
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

BackBuffer::BackBuffer() : x(0), y(0), w(0), h(0), offset(0, 0), dirty(true) {}

void BackBuffer::setSrcOffset(int x, int y)
{
    // values will be clamped when used...
    offset.x = x;
    offset.y = y;
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Glib::RefPtr<Gdk::Window> window, int newX, int newY, int newW, int newH, bool updateBackBufferSize)
{
    assert(newW && newH);

    bool newSize = w != newW || h != newH;

    x = newX;
    y = newY;
    w = newW;
    h = newH;

    // WARNING: we're assuming that the surface type won't change during all the execution time of RT. I guess it may be wrong when the user change the gfx card display settings!?
    if (updateBackBufferSize && newSize && window) {
        // allocate a new Surface
        surface.clear();  // ... don't know if this is necessary?
        surface = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, w, h);
        dirty = true;
    }

    return dirty;
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Cairo::Format format, int newX, int newY, int newW, int newH, bool updateBackBufferSize)
{
    assert(!newW && !newH);

    bool newSize = w != newW || h != newH;

    x = newX;
    y = newY;
    w = newW;
    h = newH;

    // WARNING: we're assuming that the surface type won't change during all the execution time of RT. I guess it may be wrong when the user change the gfx card display settings!?
    if (updateBackBufferSize && newSize) {
        // allocate a new Surface
        surface.clear();  // ... don't know if this is necessary?
        surface = Cairo::ImageSurface::create(format, w, h);
        dirty = true;
    }

    return dirty;
}

/*
 * Copy the backbuffer to a Gdk::Window
 */
void BackBuffer::copySurface(Glib::RefPtr<Gdk::Window> window, GdkRectangle *rectangle)
{
    if (surface && window) {
        // TODO: look out if window can be different on each call, and if not, store a reference to the window
        Cairo::RefPtr<Cairo::Context> crSrc = window->create_cairo_context();
        Cairo::RefPtr<Cairo::Surface> destSurface = crSrc->get_target();

        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        Cairo::RefPtr<Cairo::Context> crDest = Cairo::Context::create(destSurface);
        crDest->set_source(surface, x - offsetX, y - offsetY);
        crDest->set_line_width(0.);

        if (rectangle) {
            crDest->rectangle(rectangle->x, rectangle->y, rectangle->width, rectangle->height);
        } else {
            crDest->rectangle(x, y, w, h);
        }

        crDest->fill();
    }
}

/*
 * Copy the BackBuffer to another BackBuffer
 */
void BackBuffer::copySurface(BackBuffer *destBackBuffer, GdkRectangle *rectangle)
{
    if (surface && destBackBuffer) {
        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        Cairo::RefPtr<Cairo::Context> crDest = Cairo::Context::create(destBackBuffer->getSurface());
        crDest->set_source(surface, x - offsetX, y - offsetY);
        crDest->set_line_width(0.);

        if (rectangle) {
            crDest->rectangle(rectangle->x, rectangle->y, rectangle->width, rectangle->height);
        } else {
            crDest->rectangle(x, y, w, h);
        }

        crDest->fill();
    }
}

/*
 * Copy the BackBuffer to another Cairo::Surface
 */
void BackBuffer::copySurface(Cairo::RefPtr<Cairo::ImageSurface> destSurface, GdkRectangle *rectangle)
{
    if (surface && destSurface) {
        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        Cairo::RefPtr<Cairo::Context> crDest = Cairo::Context::create(destSurface);
        crDest->set_source(surface, x - offsetX, y - offsetY);
        crDest->set_line_width(0.);

        if (rectangle) {
            crDest->rectangle(rectangle->x, rectangle->y, rectangle->width, rectangle->height);
        } else {
            crDest->rectangle(x, y, w, h);
        }

        crDest->fill();
    }
}


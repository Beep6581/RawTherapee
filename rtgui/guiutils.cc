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
#include <cairomm/cairomm.h>
#include "../rtengine/rt_math.h"

#include "guiutils.h"
#include "options.h"
#include "../rtengine/rt_math.h"
#include "../rtengine/utils.h"
#include "../rtengine/procparams.h"
#include "rtimage.h"
#include "multilangmgr.h"

#include <assert.h>

//extern Glib::Threads::Thread* mainThread;

using namespace std;

Glib::RefPtr<RTImage> MyExpander::inconsistentImage;
Glib::RefPtr<RTImage> MyExpander::enabledImage;
Glib::RefPtr<RTImage> MyExpander::disabledImage;
Glib::RefPtr<RTImage> MyExpander::openedImage;
Glib::RefPtr<RTImage> MyExpander::closedImage;

IdleRegister::~IdleRegister()
{
    destroy();
}

void IdleRegister::add(std::function<bool ()> function, gint priority)
{
    const auto dispatch =
        [](gpointer data) -> gboolean
        {
            DataWrapper* const data_wrapper = static_cast<DataWrapper*>(data);

            if (!data_wrapper->function()) {
                data_wrapper->self->mutex.lock();
                data_wrapper->self->ids.erase(data_wrapper);
                data_wrapper->self->mutex.unlock();

                delete data_wrapper;
                return FALSE;
            }

            return TRUE;
        };

    DataWrapper* const data_wrapper = new DataWrapper{
        this,
        std::move(function)
    };

    mutex.lock();
    ids[data_wrapper] = gdk_threads_add_idle_full(priority, dispatch, data_wrapper, nullptr);
    mutex.unlock();
}

void IdleRegister::destroy()
{
    mutex.lock();
    for (const auto& id : ids) {
        g_source_remove(id.second);
        delete id.first;
    }
    ids.clear();
    mutex.unlock();
}

/*
gboolean giveMeAGo(void* data) {
    GThreadLock *threadMutex = static_cast<GThreadLock*>(data);
    printf("A\n");
    Glib::Threads::Mutex::Lock GUILock(threadMutex->GUI);
    printf("B\n");
    {
    Glib::Threads::Mutex::Lock operationLock(threadMutex->operation);
    printf("C\n");

    threadMutex->operationCond.signal();
    printf("D\n");
    operationLock.release();  // because we're not sure that "lock" destructor happens here...
    }
    threadMutex->GUICond.wait(threadMutex->GUI);
    printf("E\n");

    GUILock.release();

    return false;
}

GThreadLock::GThreadLock() : sameThread(false) {
    if (Glib::Threads::Thread::self() == mainThread) {
        sameThread = true;
        return;
    }

    printf("10\n");
    {
    Glib::Threads::Mutex::Lock operationLock(operation);

    printf("20\n");
    gdk_threads_add_idle(giveMeAGo, this);

    printf("30\n");
    operationCond.wait(operation);
    printf("40\n");
    operationLock.release();
    }
}

GThreadLock::~GThreadLock() {
    if (!sameThread) {
        printf("50\n");
        Glib::Threads::Mutex::Lock lock(GUI);
        printf("60\n");
        GUICond.signal();
        printf("Fin\n");
    }
}
*/

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

void setExpandAlignProperties(Gtk::Widget *widget, bool hExpand, bool vExpand, enum Gtk::Align hAlign, enum Gtk::Align vAlign)
{
    widget->set_hexpand(hExpand);
    widget->set_vexpand(vExpand);
    widget->set_halign(hAlign);
    widget->set_valign(vAlign);
}

Gtk::Border getPadding(const Glib::RefPtr<Gtk::StyleContext> style)
{
    Gtk::Border padding;
    if (!style) {
        return padding;
    }

    int s = (double)RTScalable::getScale();
    padding = style->get_padding();
    if (s > 1) {
        padding.set_left(padding.get_left() * s);
        padding.set_right(padding.get_right() * s);
        padding.set_top(padding.get_top() * s);
        padding.set_bottom(padding.get_bottom() * s);
    }

    return padding;
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

bool confirmOverwrite (Gtk::Window& parent, const std::string& filename)
{
    bool safe = true;

    if (Glib::file_test (filename, Glib::FILE_TEST_EXISTS)) {
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
    } else if (options.bgcolor == 3) {
        cr->set_source_rgb (0.467, 0.467, 0.467);
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
        cr->set_dash (std::valarray<double>({4}), 0);
        cr->move_to (rectx1, recty1);
        cr->line_to (rectx2, recty1);
        cr->line_to (rectx2, recty2);
        cr->line_to (rectx1, recty2);
        cr->line_to (rectx1, recty1);
        cr->stroke ();
        cr->set_dash (std::valarray<double>(), 0);

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
                 * The first horizontal guide is for the crown, the second is roughly for the nostrils, the third is for the chin.
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
            cr->set_dash (std::valarray<double>({4}), 0);
            cr->move_to (rectx1, recty1);
            cr->line_to (rectx2, recty2);
            cr->stroke ();
            cr->set_dash (std::valarray<double>(), 0);

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
            cr->set_dash (std::valarray<double>({4}), 0);
            cr->move_to (rectx1, recty2);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            cr->set_dash (std::valarray<double>(), 0);

            x = width - (a * b) / height;
            y = (b * (d - a)) / width;
            cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
            cr->move_to (rectx2, recty1);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
            cr->set_dash (std::valarray<double>({4}), 0);
            cr->move_to (rectx2, recty1);
            cr->line_to (rectx1 + x, recty1 + y);
            cr->stroke ();
            cr->set_dash (std::valarray<double>(), 0);
        }
    }

    cr->reset_clip ();
}

/*
bool ExpanderBox::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) {

    if (!options.useSystemTheme) {
        Glib::RefPtr<Gdk::Window> window = get_window();
        Glib::RefPtr<Gtk::StyleContext> style = get_style_context ();

        int x_, y_, w_, h_;
        window->get_geometry(x_, y_, w_, h_);
        double x = 0.;
        double y = 0.;
        double w = double(w_);
        double h = double(h_);

        cr->set_antialias (Cairo::ANTIALIAS_NONE);

        // draw a frame
        style->render_background(cr, x, y, w, h);
        / *
        cr->set_line_width (1.0);
        Gdk::RGBA c = style->get_color (Gtk::STATE_FLAG_NORMAL);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->move_to(x+0.5, y+0.5);
        cr->line_to(x+w, y+0.5);
        cr->line_to(x+w, y+h);
        cr->line_to(x+0.5, y+h);
        cr->line_to(x+0.5, y+0.5);
        cr->stroke ();
        * /
    }
    return Gtk::EventBox::on_draw(cr);
}
*/

ExpanderBox::ExpanderBox( Gtk::Container *p): pC(p)
{
    set_name ("ExpanderBox");
//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    set_border_width(2);
#endif
//GTK318
}

void ExpanderBox::setLevel(int level)
{
    if (level <= 1) {
        set_name("ExpanderBox");
    } else if (level == 2) {
        set_name("ExpanderBox2");
    } else if (level >= 3) {
        set_name("ExpanderBox3");
    }
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
    if (!inconsistentImage) {  // if one is null, all are null
        inconsistentImage = Glib::RefPtr<RTImage>(new RTImage("power-inconsistent-small.png"));
        enabledImage = Glib::RefPtr<RTImage>(new RTImage("power-on-small.png"));
        disabledImage = Glib::RefPtr<RTImage>(new RTImage("power-off-small.png"));
        openedImage = Glib::RefPtr<RTImage>(new RTImage("expander-open-small.png"));
        closedImage = Glib::RefPtr<RTImage>(new RTImage("expander-closed-small.png"));
    }
}

void MyExpander::cleanup()
{
    inconsistentImage.reset();
    enabledImage.reset();
    disabledImage.reset();
    openedImage.reset();
    closedImage.reset();
}

MyExpander::MyExpander(bool useEnabled, Gtk::Widget* titleWidget) :
    enabled(false), inconsistent(false), flushEvent(false), expBox(nullptr),
    child(nullptr), headerWidget(nullptr), statusImage(nullptr),
    label(nullptr), useEnabled(useEnabled)
{
    set_spacing(0);
    set_name("MyExpander");
    set_can_focus(false);
    setExpandAlignProperties(this, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    headerHBox = Gtk::manage( new Gtk::HBox());
    headerHBox->set_can_focus(false);
    setExpandAlignProperties(headerHBox, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    if (useEnabled) {
        statusImage = Gtk::manage(new RTImage(disabledImage));
        imageEvBox = Gtk::manage(new Gtk::EventBox());
        imageEvBox->add(*statusImage);
        imageEvBox->set_above_child(true);
        imageEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_enabled_change) );
        imageEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        imageEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        headerHBox->pack_start(*imageEvBox, Gtk::PACK_SHRINK, 0);
    } else {
        statusImage = Gtk::manage(new RTImage(openedImage));
        headerHBox->pack_start(*statusImage, Gtk::PACK_SHRINK, 0);
    }

    statusImage->set_can_focus(false);

    if (titleWidget) {
        setExpandAlignProperties(titleWidget, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
        headerHBox->pack_start(*titleWidget, Gtk::PACK_EXPAND_WIDGET, 0);
        headerWidget = titleWidget;
    }

    titleEvBox = Gtk::manage(new Gtk::EventBox());
    titleEvBox->set_name("MyExpanderTitle");
    titleEvBox->set_border_width(2);
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
    enabled(false), inconsistent(false), flushEvent(false), expBox(nullptr),
    child(nullptr), headerWidget(nullptr),
    label(nullptr), useEnabled(useEnabled)
{
    set_spacing(0);
    set_name("MyExpander");
    set_can_focus(false);
    setExpandAlignProperties(this, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    headerHBox = Gtk::manage( new Gtk::HBox());
    headerHBox->set_can_focus(false);
    setExpandAlignProperties(headerHBox, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);


    if (useEnabled) {
        statusImage = Gtk::manage(new RTImage(disabledImage));
        imageEvBox = Gtk::manage(new Gtk::EventBox());
        imageEvBox->set_name("MyExpanderStatus");
        imageEvBox->add(*statusImage);
        imageEvBox->set_above_child(true);
        imageEvBox->signal_button_release_event().connect( sigc::mem_fun(this, & MyExpander::on_enabled_change) );
        imageEvBox->signal_enter_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        imageEvBox->signal_leave_notify_event().connect( sigc::mem_fun(this, & MyExpander::on_enter_leave_enable), false );
        headerHBox->pack_start(*imageEvBox, Gtk::PACK_SHRINK, 0);
    } else {
        statusImage = Gtk::manage(new RTImage(openedImage));
        headerHBox->pack_start(*statusImage, Gtk::PACK_SHRINK, 0);
    }

    statusImage->set_can_focus(false);

    label = Gtk::manage(new Gtk::Label());
    setExpandAlignProperties(label, true, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    label->set_markup(Glib::ustring("<b>") + escapeHtmlChars(titleLabel) + Glib::ustring("</b>"));
    headerHBox->pack_start(*label, Gtk::PACK_EXPAND_WIDGET, 0);

    titleEvBox = Gtk::manage(new Gtk::EventBox());
    titleEvBox->set_name("MyExpanderTitle");
    titleEvBox->set_border_width(2);
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
    updateVScrollbars(options.hideTPVScrollbar);

//GTK318
#if GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION < 20
    headerHBox->set_spacing(2);
    headerHBox->set_border_width(1);
    set_spacing(0);
    set_border_width(0);
#endif
//GTK318
}

void MyExpander::updateVScrollbars(bool hide)
{
    if (hide) {
        get_style_context()->remove_class("withScrollbar");
    } else {
        get_style_context()->add_class("withScrollbar");
    }
}

void MyExpander::setLevel (int level)
{
    if (expBox) {
        expBox->setLevel(level);
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
                statusImage->set(inconsistentImage->get_surface());
            } else {
                if (enabled) {
                    statusImage->set(enabledImage->get_surface());
                } else {
                    statusImage->set(disabledImage->get_surface());
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
                    statusImage->set(disabledImage->get_surface());
                    message.emit();
                }
            } else {
                enabled = true;

                if (!inconsistent) {
                    statusImage->set(enabledImage->get_surface());
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
            statusImage->set(openedImage->get_surface());
        } else {
            statusImage->set(closedImage->get_surface());
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

void MyExpander::add  (Gtk::Container& widget, bool setChild)
{
    if(setChild) {
        child = &widget;
    }
    expBox = Gtk::manage (new ExpanderBox (child));
    expBox->add (widget);
    pack_start(*expBox, Gtk::PACK_SHRINK, 0);
    widget.show();
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
            statusImage->set(closedImage->get_surface());
        } else {
            statusImage->set(openedImage->get_surface());
        }
    }

    if (isVisible) {
        expBox->hideBox();
    } else {
        expBox->showBox();
    }

    return false;
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
            statusImage->set(disabledImage->get_surface());
        } else {
            enabled = true;
            statusImage->set(enabledImage->get_surface());
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
}

bool MyScrolledWindow::on_scroll_event (GdkEventScroll* event)
{
    if (!options.hideTPVScrollbar) {
        Gtk::ScrolledWindow::on_scroll_event (event);
        return true;
    }

    Glib::RefPtr<Gtk::Adjustment> adjust = get_vadjustment();
    Gtk::Scrollbar *scroll = get_vscrollbar();

    if (adjust && scroll) {
        const double upperBound = adjust->get_upper();
        const double lowerBound = adjust->get_lower();
        double value = adjust->get_value();
        double step  = adjust->get_step_increment();
        double value2 = 0.;

        if (event->direction == GDK_SCROLL_DOWN) {
            value2 = rtengine::min<double>(value + step, upperBound);

            if (value2 != value) {
                scroll->set_value(value2);
            }
        } else if (event->direction == GDK_SCROLL_UP) {
            value2 = rtengine::max<double>(value - step, lowerBound);

            if (value2 != value) {
                scroll->set_value(value2);
            }
        } else if (event->direction == GDK_SCROLL_SMOOTH) {
            value2 = rtengine::LIM<double>(value + event->delta_y * step, lowerBound, upperBound);

            if (value2 != value) {
                scroll->set_value(value2);
            }
        }
    }

    return true;
}

void MyScrolledWindow::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    natural_width = minimum_width = 100 * RTScalable::getScale();
}

void MyScrolledWindow::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    natural_height = minimum_height = 50 * RTScalable::getScale();
}

void MyScrolledWindow::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    natural_height = minimum_height = 50 * RTScalable::getScale();
}

/*
 *
 * Derived class of some widgets to properly handle the scroll wheel ;
 * the user has to use the Shift key to be able to change the widget's value,
 * otherwise the mouse wheel will scroll the toolbar.
 *
 */
MyScrolledToolbar::MyScrolledToolbar ()
{
    set_policy (Gtk::POLICY_EXTERNAL, Gtk::POLICY_NEVER);
    get_style_context()->add_class("scrollableToolbar");

    // Works fine with Gtk 3.22, but a a custom made get_preferred_height had to be created as a workaround
    // taken from the official Gtk3.22 source code
    //set_propagate_natural_height(true);
}

bool MyScrolledToolbar::on_scroll_event (GdkEventScroll* event)
{
    Glib::RefPtr<Gtk::Adjustment> adjust = get_hadjustment();
    Gtk::Scrollbar *scroll = get_hscrollbar();

    if (adjust && scroll) {
        const double upperBound = adjust->get_upper();
        const double lowerBound = adjust->get_lower();
        double value = adjust->get_value();
        double step  = adjust->get_step_increment() * 2;
        double value2 = 0.;

//        printf("MyScrolledToolbar::on_scroll_event / delta_x=%.5f, delta_y=%.5f, direction=%d, type=%d, send_event=%d\n",
//                event->delta_x, event->delta_y, (int)event->direction, (int)event->type, event->send_event);

        if (event->direction == GDK_SCROLL_DOWN) {
            value2 = rtengine::min<double>(value + step, upperBound);
            if (value2 != value) {
                scroll->set_value(value2);
            }
        } else if (event->direction == GDK_SCROLL_UP) {
            value2 = rtengine::max<double>(value - step, lowerBound);
            if (value2 != value) {
                scroll->set_value(value2);
            }
        } else if (event->direction == GDK_SCROLL_SMOOTH) {
            if (event->delta_x) {  // if the user use a pad, it can scroll horizontally
                value2 = rtengine::LIM<double>(value + (event->delta_x > 0 ? 30 : -30), lowerBound, upperBound);
            } else if (event->delta_y) {
                value2 = rtengine::LIM<double>(value + (event->delta_y > 0 ? 30 : -30), lowerBound, upperBound);
            }
            if (value2 != value) {
                scroll->set_value(value2);
            }
        }
    }

    return true;
}

void MyScrolledToolbar::get_preferred_height_vfunc (int &minimumHeight, int &naturalHeight) const
{
    int currMinHeight = 0;
    int currNatHeight = 0;
    std::vector<const Widget*> childs = get_children();
    minimumHeight = naturalHeight = 0;

    for (auto child : childs)
    {
        if(child->is_visible()) {
            child->get_preferred_height(currMinHeight, currNatHeight);
            minimumHeight = rtengine::max(currMinHeight, minimumHeight);
            naturalHeight = rtengine::max(currNatHeight, naturalHeight);
        }
    }
}

MyComboBoxText::MyComboBoxText (bool has_entry) : Gtk::ComboBoxText(has_entry)
{
    minimumWidth = naturalWidth = 70;
    Gtk::CellRendererText* cellRenderer = dynamic_cast<Gtk::CellRendererText*>(get_first_cell());
    cellRenderer->property_ellipsize() = Pango::ELLIPSIZE_MIDDLE;
    add_events(Gdk::SCROLL_MASK|Gdk::SMOOTH_SCROLL_MASK);
}

bool MyComboBoxText::on_scroll_event (GdkEventScroll* event)
{

//    printf("MyComboboxText::on_scroll_event / delta_x=%.5f, delta_y=%.5f, direction=%d, type=%d, send_event=%d\n",
//            event->delta_x, event->delta_y, (int)event->direction, (int)event->type, event->send_event);
    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::ComboBoxText::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

void MyComboBoxText::setPreferredWidth (int minimum_width, int natural_width)
{
    if (natural_width == -1 && minimum_width == -1) {
        naturalWidth = minimumWidth = 70 * RTScalable::getScale();
    } else if (natural_width == -1) {
        naturalWidth =  minimumWidth = minimum_width;
    } else if (minimum_width == -1) {
        naturalWidth = natural_width;
        minimumWidth = rtengine::max(naturalWidth / 2, 20);
        minimumWidth = rtengine::min(naturalWidth, minimumWidth);
    } else {
        naturalWidth = natural_width;
        minimumWidth = minimum_width;
    }
}

void MyComboBoxText::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    natural_width = rtengine::max(naturalWidth, 10 * RTScalable::getScale());
    minimum_width = rtengine::max(minimumWidth, 10 * RTScalable::getScale());
}
void MyComboBoxText::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    natural_width = rtengine::max(naturalWidth, 10 * RTScalable::getScale());
    minimum_width = rtengine::max(minimumWidth, 10 * RTScalable::getScale());
}


MyComboBox::MyComboBox ()
{
    minimumWidth = naturalWidth = 70 * RTScalable::getScale();
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

void MyComboBox::setPreferredWidth (int minimum_width, int natural_width)
{
    if (natural_width == -1 && minimum_width == -1) {
        naturalWidth = minimumWidth = 70 * RTScalable::getScale();
    } else if (natural_width == -1) {
        naturalWidth =  minimumWidth = minimum_width;
    } else if (minimum_width == -1) {
        naturalWidth = natural_width;
        minimumWidth = rtengine::max(naturalWidth / 2, 20);
        minimumWidth = rtengine::min(naturalWidth, minimumWidth);
    } else {
        naturalWidth = natural_width;
        minimumWidth = minimum_width;
    }
}

void MyComboBox::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    natural_width = rtengine::max(naturalWidth, 10 * RTScalable::getScale());
    minimum_width = rtengine::max(minimumWidth, 10 * RTScalable::getScale());
}
void MyComboBox::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    natural_width = rtengine::max(naturalWidth, 10 * RTScalable::getScale());
    minimum_width = rtengine::max(minimumWidth, 10 * RTScalable::getScale());
}

MySpinButton::MySpinButton ()
{
    Gtk::Border border;
    border.set_bottom(0);
    border.set_top(0);
    border.set_left(3);
    border.set_right(3);
    set_inner_border(border);
    set_numeric(true);
    set_wrap(false);
    set_alignment(Gtk::ALIGN_END);
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
    set_max_width_chars(maxLen);
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
            event->keyval = GDK_KEY_period;
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

//    printf("MyHScale::on_scroll_event / delta_x=%.5f, delta_y=%.5f, direction=%d, type=%d, send_event=%d\n",
//            event->delta_x, event->delta_y, (int)event->direction, (int)event->type, event->send_event);
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

MyFileChooserButton::MyFileChooserButton(const Glib::ustring &title, Gtk::FileChooserAction action):
    title_(title),
    action_(action),
    lbl_("", Gtk::ALIGN_START),
    show_hidden_(false)
{
    lbl_.set_ellipsize(Pango::ELLIPSIZE_MIDDLE);
    lbl_.set_justify(Gtk::JUSTIFY_LEFT);
    set_none();
    box_.pack_start(lbl_, true, true);
    Gtk::Image *img = Gtk::manage(new Gtk::Image());
    img->set_from_icon_name("folder-open", Gtk::ICON_SIZE_BUTTON);
    box_.pack_start(*Gtk::manage(new Gtk::VSeparator()), false, false, 5);
    box_.pack_start(*img, false, false);
    box_.show_all_children();
    add(box_);
    signal_clicked().connect(sigc::mem_fun(*this, &MyFileChooserButton::show_chooser));

    if (GTK_MINOR_VERSION < 20) {
        set_border_width(2); // margin doesn't work on GTK < 3.20
    }

    set_name("MyFileChooserButton");
}


void MyFileChooserButton::show_chooser()
{
    Gtk::FileChooserDialog dlg(getToplevelWindow(this), title_, action_);
    dlg.add_button(M("GENERAL_CANCEL"), Gtk::RESPONSE_CANCEL);
    dlg.add_button(M(action_ == Gtk::FILE_CHOOSER_ACTION_SAVE ? "GENERAL_SAVE" : "GENERAL_OPEN"), Gtk::RESPONSE_OK);
    dlg.set_filename(filename_);
    for (auto &f : file_filters_) {
        dlg.add_filter(f);
    }
    if (cur_filter_) {
        dlg.set_filter(cur_filter_);
    }
    for (auto &f : shortcut_folders_) {
        dlg.add_shortcut_folder(f);
    }
    if (!current_folder_.empty()) {
        dlg.set_current_folder(current_folder_);
    }
    dlg.set_show_hidden(show_hidden_);
    int res = dlg.run();
    if (res == Gtk::RESPONSE_OK) {
        filename_ = dlg.get_filename();
        current_folder_ = dlg.get_current_folder();
        lbl_.set_label(Glib::path_get_basename(filename_));
        selection_changed_.emit();
    }
}


sigc::signal<void> &MyFileChooserButton::signal_selection_changed()
{
    return selection_changed_;
}


sigc::signal<void> &MyFileChooserButton::signal_file_set()
{
    return selection_changed_;
}


std::string MyFileChooserButton::get_filename() const
{
    return filename_;
}


bool MyFileChooserButton::set_filename(const std::string &filename)
{
    filename_ = filename;
    if (Glib::file_test(filename_, Glib::FILE_TEST_EXISTS)) {
        lbl_.set_label(Glib::path_get_basename(filename_));
    } else {
        set_none();
    }
    return true;
}


void MyFileChooserButton::add_filter(const Glib::RefPtr<Gtk::FileFilter> &filter)
{
    file_filters_.push_back(filter);
}


void MyFileChooserButton::remove_filter(const Glib::RefPtr<Gtk::FileFilter> &filter)
{
    auto it = std::find(file_filters_.begin(), file_filters_.end(), filter);
    if (it != file_filters_.end()) {
        file_filters_.erase(it);
    }
}


void MyFileChooserButton::set_filter(const Glib::RefPtr<Gtk::FileFilter> &filter)
{
    cur_filter_ = filter;
}


std::vector<Glib::RefPtr<Gtk::FileFilter>> MyFileChooserButton::list_filters()
{
    return file_filters_;
}


bool MyFileChooserButton::set_current_folder(const std::string &filename)
{
    current_folder_ = filename;
    if (action_ == Gtk::FILE_CHOOSER_ACTION_SELECT_FOLDER) {
        set_filename(filename);
    }
    return true;
}

std::string MyFileChooserButton::get_current_folder() const
{
    return current_folder_;
}


bool MyFileChooserButton::add_shortcut_folder(const std::string &folder)
{
    shortcut_folders_.push_back(folder);
    return true;
}


bool MyFileChooserButton::remove_shortcut_folder(const std::string &folder)
{
    auto it = std::find(shortcut_folders_.begin(), shortcut_folders_.end(), folder);
    if (it != shortcut_folders_.end()) {
        shortcut_folders_.erase(it);
    }
    return true;
}


void MyFileChooserButton::unselect_all()
{
    filename_ = "";
    set_none();
}


void MyFileChooserButton::unselect_filename(const std::string &filename)
{
    if (filename_ == filename) {
        unselect_all();
    }
}


void MyFileChooserButton::set_show_hidden(bool yes)
{
    show_hidden_ = yes;
}


void MyFileChooserButton::set_none()
{
    lbl_.set_label(Glib::ustring("(") + M("GENERAL_NONE") + ")");
}

// For an unknown reason (a bug ?), it doesn't work when action = FILE_CHOOSER_ACTION_SELECT_FOLDER !
bool MyFileChooserButton::on_scroll_event (GdkEventScroll* event)
{

    // If Shift is pressed, the widget is modified
    if (event->state & GDK_SHIFT_MASK) {
        Gtk::Button::on_scroll_event(event);
        return true;
    }

    // ... otherwise the scroll event is sent back to an upper level
    return false;
}

void MyFileChooserButton::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = natural_width = 35 * RTScalable::getScale();
}
void MyFileChooserButton::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    minimum_width = natural_width = 35 * RTScalable::getScale();
}



TextOrIcon::TextOrIcon (const Glib::ustring &fname, const Glib::ustring &labelTx, const Glib::ustring &tooltipTx)
{

    RTImage *img = Gtk::manage(new RTImage(fname));
    pack_start(*img, Gtk::PACK_SHRINK, 0);
    set_tooltip_markup("<span font_size=\"large\" font_weight=\"bold\">" + labelTx  + "</span>\n" + tooltipTx);

    set_name("TextOrIcon");
    show_all();

}

MyImageMenuItem::MyImageMenuItem(Glib::ustring label, Glib::ustring imageFileName)
{
    box = Gtk::manage (new Gtk::Grid());
    this->label = Gtk::manage( new Gtk::Label(label));
    box->set_orientation(Gtk::ORIENTATION_HORIZONTAL);

    if (!imageFileName.empty()) {
        image = Gtk::manage( new RTImage(imageFileName) );
        box->attach_next_to(*image, Gtk::POS_LEFT, 1, 1);
    } else {
        image = nullptr;
    }

    box->attach_next_to(*this->label, Gtk::POS_RIGHT, 1, 1);
    box->set_column_spacing(4);
    box->set_row_spacing(0);
    add(*box);
}

const RTImage *MyImageMenuItem::getImage () const
{
    return image;
}

const Gtk::Label* MyImageMenuItem::getLabel () const
{
    return label;
}

MyProgressBar::MyProgressBar(int width) : w(rtengine::max(width, 10 * RTScalable::getScale())) {}
MyProgressBar::MyProgressBar() : w(200 * RTScalable::getScale()) {}

void MyProgressBar::setPreferredWidth(int width)
{
    w = rtengine::max(width, 10 * RTScalable::getScale());
}

void MyProgressBar::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = rtengine::max(w / 2, 50 * RTScalable::getScale());
    natural_width = rtengine::max(w, 50 * RTScalable::getScale());
}

void MyProgressBar::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

BackBuffer::BackBuffer() : x(0), y(0), w(0), h(0), offset(0, 0), dirty(true) {}
BackBuffer::BackBuffer(int width, int height, Cairo::Format format) : x(0), y(0), w(width), h(height), offset(0, 0), dirty(true)
{
    if (w > 0 && h > 0) {
        surface = Cairo::ImageSurface::create(format, w, h);
    } else {
        w = h = 0;
    }
}

void BackBuffer::setDestPosition(int x, int y)
{
    // values will be clamped when used...
    this->x = x;
    this->y = y;
}

void BackBuffer::setSrcOffset(int x, int y)
{
    // values will be clamped when used...
    offset.set(x, y);
}

void BackBuffer::setSrcOffset(const rtengine::Coord &newOffset)
{
    // values will be clamped when used...
    offset = newOffset;
}

void BackBuffer::getSrcOffset(int &x, int &y)
{
    // values will be clamped when used...
    offset.get(x, y);
}

void BackBuffer::getSrcOffset(rtengine::Coord &offset)
{
    // values will be clamped when used...
    offset = this->offset;
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Glib::RefPtr<Gdk::Window> window, Gdk::Rectangle &rectangle, bool updateBackBufferSize)
{
    return setDrawRectangle(window, rectangle.get_x(), rectangle.get_y(), rectangle.get_width(), rectangle.get_height(), updateBackBufferSize);
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Glib::RefPtr<Gdk::Window> window, int newX, int newY, int newW, int newH, bool updateBackBufferSize)
{
    assert(newW && newH);

    bool newSize = (newW > 0 && w != newW) || (newH > 0 && h != newH);

    x = newX;
    y = newY;
    if (newW > 0) {
        w = newW;
    }
    if (newH > 0) {
        h = newH;
    }

    // WARNING: we're assuming that the surface type won't change during all the execution time of RT. I guess it may be wrong when the user change the gfx card display settings!?
    if (((updateBackBufferSize && newSize) || !surface) && window) {
        // allocate a new Surface
        surface.clear();  // ... don't know if this is necessary?
        surface = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, w, h);
        dirty = true;
    }

    return dirty;
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Cairo::Format format, Gdk::Rectangle &rectangle, bool updateBackBufferSize)
{
    return setDrawRectangle(format, rectangle.get_x(), rectangle.get_y(), rectangle.get_width(), rectangle.get_height(), updateBackBufferSize);
}

// Note: newW & newH must be > 0
bool BackBuffer::setDrawRectangle(Cairo::Format format, int newX, int newY, int newW, int newH, bool updateBackBufferSize)
{
    assert(newW && newH);

    bool newSize = (newW > 0 && w != newW) || (newH > 0 && h != newH);

    x = newX;
    y = newY;
    if (newW > 0) {
        w = newW;
    }
    if (newH > 0) {
        h = newH;
    }

    // WARNING: we're assuming that the surface type won't change during all the execution time of RT. I guess it may be wrong when the user change the gfx card display settings!?
    if ((updateBackBufferSize && newSize) || !surface) {
        // allocate a new Surface
        surface.clear();  // ... don't know if this is necessary?
        surface = Cairo::ImageSurface::create(format, w, h);
        dirty = true;
    }

    return dirty;
}

/*
 * Copy uint8 RGB raw data to an ImageSurface. We're assuming that the source contains enough data for the given srcX, srcY, srcW, srcH -> no error checking!
 */
void BackBuffer::copyRGBCharData(const unsigned char *srcData, int srcX, int srcY, int srcW, int srcH, int srcRowStride, int dstX, int dstY)
{
    unsigned char r, g, b;

    if (!surface) {
        return;
    }

    //printf("copyRGBCharData:    src: (X:%d Y:%d, W:%d H:%d)  /  dst: (X: %d Y:%d)\n", srcX, srcY, srcW, srcH, dstX, dstY);

    unsigned char *dstData = surface->get_data();
    int surfW = surface->get_width();
    int surfH = surface->get_height();

    if (!srcData || dstX >= surfW || dstY >= surfH || srcW <= 0 || srcH <= 0 || srcX < 0 || srcY < 0) {
        return;
    }

    for (int i = 0; i < srcH; ++i) {
        if (dstY + i >= surfH) {
            break;
        }

        const unsigned char *src = srcData + i * srcRowStride;
        unsigned char *dst = dstData + ((dstY + i) * surfW + dstX) * 4;

        for (int j = 0; j < srcW; ++j) {
            if (dstX + j >= surfW) {
                break;
            }

            r = *(src++);
            g = *(src++);
            b = *(src++);

            rtengine::poke255_uc(dst, r, g, b);
        }
    }

    surface->mark_dirty();

}

/*
 * Copy the backbuffer to a Gdk::Window
 */
void BackBuffer::copySurface(Glib::RefPtr<Gdk::Window> window, Gdk::Rectangle *destRectangle)
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
        crDest->set_line_width(0.);

        if (destRectangle) {
            crDest->set_source(surface, -offsetX + destRectangle->get_x(), -offsetY + destRectangle->get_y());
            int w_ = destRectangle->get_width() > 0 ? destRectangle->get_width() : w;
            int h_ = destRectangle->get_height() > 0 ? destRectangle->get_height() : h;
            //printf("BackBuffer::copySurface / rectangle1(%d, %d, %d, %d)\n", destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            crDest->rectangle(destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            //printf("BackBuffer::copySurface / rectangle1\n");
        } else {
            crDest->set_source(surface, -offsetX + x, -offsetY + y);
            //printf("BackBuffer::copySurface / rectangle2(%d, %d, %d, %d)\n", x, y, w, h);
            crDest->rectangle(x, y, w, h);
            //printf("BackBuffer::copySurface / rectangle2\n");
        }

        crDest->fill();
    }
}

/*
 * Copy the BackBuffer to another BackBuffer
 */
void BackBuffer::copySurface(BackBuffer *destBackBuffer, Gdk::Rectangle *destRectangle)
{
    if (surface && destBackBuffer) {
        Cairo::RefPtr<Cairo::ImageSurface> destSurface = destBackBuffer->getSurface();

        if (!destSurface) {
            return;
        }

        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        Cairo::RefPtr<Cairo::Context> crDest = Cairo::Context::create(destSurface);
        crDest->set_line_width(0.);

        if (destRectangle) {
            crDest->set_source(surface, -offsetX + destRectangle->get_x(), -offsetY + destRectangle->get_y());
            int w_ = destRectangle->get_width() > 0 ? destRectangle->get_width() : w;
            int h_ = destRectangle->get_height() > 0 ? destRectangle->get_height() : h;
            //printf("BackBuffer::copySurface / rectangle3(%d, %d, %d, %d)\n", destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            crDest->rectangle(destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            //printf("BackBuffer::copySurface / rectangle3\n");
        } else {
            crDest->set_source(surface, -offsetX + x, -offsetY + y);
            //printf("BackBuffer::copySurface / rectangle4(%d, %d, %d, %d)\n", x, y, w, h);
            crDest->rectangle(x, y, w, h);
            //printf("BackBuffer::copySurface / rectangle4\n");
        }

        crDest->fill();
    }
}

/*
 * Copy the BackBuffer to another Cairo::Surface
 */
void BackBuffer::copySurface(Cairo::RefPtr<Cairo::ImageSurface> destSurface, Gdk::Rectangle *destRectangle)
{
    if (surface && destSurface) {
        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        Cairo::RefPtr<Cairo::Context> crDest = Cairo::Context::create(destSurface);
        crDest->set_line_width(0.);

        if (destRectangle) {
            crDest->set_source(surface, -offsetX + destRectangle->get_x(), -offsetY + destRectangle->get_y());
            int w_ = destRectangle->get_width() > 0 ? destRectangle->get_width() : w;
            int h_ = destRectangle->get_height() > 0 ? destRectangle->get_height() : h;
            //printf("BackBuffer::copySurface / rectangle5(%d, %d, %d, %d)\n", destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            crDest->rectangle(destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            //printf("BackBuffer::copySurface / rectangle5\n");
        } else {
            crDest->set_source(surface, -offsetX + x, -offsetY + y);
            //printf("BackBuffer::copySurface / rectangle6(%d, %d, %d, %d)\n", x, y, w, h);
            crDest->rectangle(x, y, w, h);
            //printf("BackBuffer::copySurface / rectangle6\n");
        }

        crDest->fill();
    }
}

/*
 * Copy the BackBuffer to another Cairo::Surface
 */
void BackBuffer::copySurface(Cairo::RefPtr<Cairo::Context> crDest, Gdk::Rectangle *destRectangle)
{
    if (surface && crDest) {
        // compute the source offset
        int offsetX = rtengine::LIM<int>(offset.x, 0, surface->get_width());
        int offsetY = rtengine::LIM<int>(offset.y, 0, surface->get_height());

        // now copy the off-screen Surface to the destination Surface
        // int srcSurfW = surface->get_width();
        // int srcSurfH = surface->get_height();
        //printf("srcSurf:  w: %d, h: %d\n", srcSurfW, srcSurfH);
        crDest->set_line_width(0.);

        if (destRectangle) {
            crDest->set_source(surface, -offsetX + destRectangle->get_x(), -offsetY + destRectangle->get_y());
            int w_ = destRectangle->get_width() > 0 ? destRectangle->get_width() : w;
            int h_ = destRectangle->get_height() > 0 ? destRectangle->get_height() : h;
            //printf("BackBuffer::copySurface / rectangle7(%d, %d, %d, %d)\n", destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            crDest->rectangle(destRectangle->get_x(), destRectangle->get_y(), w_, h_);
            //printf("BackBuffer::copySurface / rectangle7\n");
        } else {
            crDest->set_source(surface, -offsetX + x, -offsetY + y);
            //printf("BackBuffer::copySurface / rectangle8(%d, %d, %d, %d)\n", x, y, w, h);
            crDest->rectangle(x, y, w, h);
            //printf("BackBuffer::copySurface / rectangle8\n");
        }

        crDest->fill();
    }
}

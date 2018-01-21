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
#include "rawcrop.h"
#include "options.h"
#include "rtimage.h"
#include "eventmapper.h"
#include "../rtengine/refreshmap.h"

using namespace rtengine;
using namespace rtengine::procparams;
using namespace rtedit;

RawCrop::RawCrop():
    FoldableToolPanel(this, "rawcrop", M("TP_RAW_CROP_LABEL"), false, true),
    EditSubscriber(ET_OBJECTS),
    litArea(LitArea::none),
    xSet(false),
    ySet(false),
    wSet(false),
    hSet(false)
{

    maxw = 50000;
    maxh = 50000;

    Gtk::Grid *grid = Gtk::manage (new Gtk::Grid());
    setExpandAlignProperties(grid, false, false, Gtk::ALIGN_FILL, Gtk::ALIGN_START);
    grid->set_column_spacing(4);
    grid->set_row_spacing(2);

    // ---------------- Edit button -------------

    edit = Gtk::manage (new Gtk::ToggleButton());
    setExpandAlignProperties(edit, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    edit->add (*Gtk::manage (new RTImage ("editmodehand.png")));
    edit->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
    grid->attach_next_to(*edit, Gtk::POS_LEFT, 4, 1);

    // ---------------- X -------------

    Gtk::Label *xLabel = Gtk::manage (new Gtk::Label (M("TP_CROP_X") + ":"));
    setExpandAlignProperties(xLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    grid->attach_next_to(*xLabel, *edit, Gtk::POS_BOTTOM, 1, 1);

    x = Gtk::manage (new MySpinButton ());
    x->set_size_request (60, -1);
    setExpandAlignProperties(x, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    grid->attach_next_to(*x, *xLabel, Gtk::POS_RIGHT, 1, 1);

    // ---------------- Y -------------

    Gtk::Label *yLabel = Gtk::manage (new Gtk::Label (M("TP_CROP_Y") + ":"));
    setExpandAlignProperties(yLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    grid->attach_next_to(*yLabel, *x, Gtk::POS_RIGHT, 1, 1);

    y = Gtk::manage (new MySpinButton ());
    y->set_size_request (60, -1);
    setExpandAlignProperties(y, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    grid->attach_next_to(*y, *yLabel, Gtk::POS_RIGHT, 1, 1);

    // ---------------- W -------------

    Gtk::Label *wLabel = Gtk::manage (new Gtk::Label (M("TP_CROP_W") + ":"));
    setExpandAlignProperties(wLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    grid->attach_next_to(*wLabel, *xLabel, Gtk::POS_BOTTOM, 1, 1);

    w = Gtk::manage (new MySpinButton ());
    w->set_size_request (60, -1);
    setExpandAlignProperties(w, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    grid->attach_next_to(*w, *wLabel, Gtk::POS_RIGHT, 1, 1);

    // ---------------- H -------------

    Gtk::Label *hLabel = Gtk::manage (new Gtk::Label (M("TP_CROP_H") + ":"));
    setExpandAlignProperties(hLabel, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    grid->attach_next_to(*hLabel, *w, Gtk::POS_RIGHT, 1, 1);

    h = Gtk::manage (new MySpinButton ());
    h->set_size_request (60, -1);
    setExpandAlignProperties(h, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    grid->attach_next_to(*h, *hLabel, Gtk::POS_RIGHT, 1, 1);

    // --------------------------------

    pack_start (*grid, Gtk::PACK_EXPAND_WIDGET, 2);

    w->set_range (1, maxw);
    h->set_range (1, maxh);
    x->set_range (0, maxw - 1);
    y->set_range (0, maxh - 1);

    x->set_digits (0);
    x->set_increments (1, 100);
    x->set_value (0);

    y->set_digits (0);
    y->set_increments (1, 100);
    y->set_value (0);

    w->set_digits (0);
    w->set_increments (1, 100);
    w->set_value (200);

    h->set_digits (0);
    h->set_increments (1, 100);
    h->set_value (200);

    editConn = edit->signal_toggled().connect( sigc::mem_fun(*this, &RawCrop::editToggled) );
    xconn = x->signal_value_changed().connect ( std::bind(sigc::mem_fun(*this, &RawCrop::spinChanged), x), true);
    yconn = y->signal_value_changed().connect ( std::bind(sigc::mem_fun(*this, &RawCrop::spinChanged), y), true);
    wconn = w->signal_value_changed().connect ( std::bind(sigc::mem_fun(*this, &RawCrop::spinChanged), w), true);
    hconn = h->signal_value_changed().connect ( std::bind(sigc::mem_fun(*this, &RawCrop::spinChanged), h), true);

    // ---------------- ProcEvents -------------

    auto m = ProcEventMapper::getInstance();
    EvRawCropEnabled = m->newEvent(RAWCROP, "TP_RAW_CROP_LABEL");
    EvRawCropDims = m->newEvent(RAWCROP, "HISTORY_MSG_RAW_CROP_DIM");
    EvRawCropDimsOPA = m->newEvent(MINUPDATE, "");

    // ---------------- On preview geometry -------------

    // Instantiating the Editing geometry; positions will be initialized later
    rtengine::Coord coord(0,0);

    // Visible geometry
    centerCircle.center = coord;
    centerCircle.radius = 6;
    centerCircle.datum = Geometry::IMAGE;
    centerCircle.radiusInImageSpace = false;
    centerCircle.filled = true;

    topLeft.points.resize(3);
    topLeft.points[0] = coord;
    topLeft.points[1] = coord;
    topLeft.points[2] = coord;
    topLeft.datum = Geometry::IMAGE;
    topLeft.setHoverable(false);
    topLeft.setAutoColor(true);
    topLeft.filled = false;

    topRight.points.resize(3);
    topRight.points[0] = coord;
    topRight.points[1] = coord;
    topRight.points[2] = coord;
    topRight.datum = Geometry::IMAGE;
    topRight.setHoverable(false);
    topRight.setAutoColor(true);
    topRight.filled = false;

    bottomRight.points.resize(3);
    bottomRight.points[0] = coord;
    bottomRight.points[1] = coord;
    bottomRight.points[2] = coord;
    bottomRight.datum = Geometry::IMAGE;
    bottomRight.setHoverable(false);
    bottomRight.setAutoColor(true);
    bottomRight.filled = false;

    bottomLeft.points.resize(3);
    bottomLeft.points[0] = coord;
    bottomLeft.points[1] = coord;
    bottomLeft.points[2] = coord;
    bottomLeft.datum = Geometry::IMAGE;
    bottomLeft.setHoverable(false);
    bottomLeft.setAutoColor(true);
    bottomLeft.filled = false;

    side.points.resize(4);
    side.points[0] = coord;
    side.points[1] = coord;
    side.points[2] = coord;
    side.points[3] = coord;
    side.datum = Geometry::IMAGE;
    side.setHoverable(false);
    side.setAutoColor(true);
    side.filled = false;
    side.setVisible(false);

    EditSubscriber::visibleGeometry.push_back(&centerCircle);
    EditSubscriber::visibleGeometry.push_back(&topLeft);
    EditSubscriber::visibleGeometry.push_back(&topRight);
    EditSubscriber::visibleGeometry.push_back(&bottomRight);
    EditSubscriber::visibleGeometry.push_back(&bottomLeft);
    EditSubscriber::visibleGeometry.push_back(&side);

    // --------------------------------

    show_all ();
}

RawCrop::~RawCrop()
{
    idle_register.destroy();
}

void RawCrop::read (const ProcParams* pp, const ParamsEdited* pedited)
{

    ConnectionBlocker xBlocker(xconn);
    ConnectionBlocker yBlocker(yconn);
    ConnectionBlocker wBlocker(wconn);
    ConnectionBlocker hBlocker(hconn);

    disableListener ();

    setEnabled(pp->raw.rawCrop);

    // check if the new values are larger than the maximum
    double tmp, maxw2, maxh2;
    w->get_range (tmp, maxw2);
    h->get_range (tmp, maxh2);

    if (pp->raw.rcX + pp->raw.rcWidth > (int)maxw2 || pp->raw.rcY + pp->raw.rcHeight > (int)maxh2) {
        maxw = pp->raw.rcX + pp->raw.rcWidth;
        maxh = pp->raw.rcY + pp->raw.rcHeight;
        w->set_range (1, maxw);
        h->set_range (1, maxh);
        x->set_range (0, maxw - 1);
        y->set_range (0, maxh - 1);
    }

    x->set_value(pp->raw.rcX);
    y->set_value(pp->raw.rcY);
    w->set_value(std::max(pp->raw.rcWidth, 1));
    h->set_value(std::max(pp->raw.rcHeight, 1));

    if (pedited) {
        wSet = pedited->raw.rcWidth;
        hSet = pedited->raw.rcHeight;
        xSet = pedited->raw.rcX;
        ySet = pedited->raw.rcY;

        set_inconsistent (multiImage && !pedited->raw.rawCrop);
    }

    updateGeometry(pp->raw.rcX, pp->raw.rcY, pp->raw.rcWidth, pp->raw.rcHeight);

    enableListener ();
}

void RawCrop::write (ProcParams* pp, ParamsEdited* pedited)
{

    pp->raw.rawCrop = getEnabled ();
    pp->raw.rcX = x->get_value_as_int();
    pp->raw.rcY = y->get_value_as_int();
    pp->raw.rcWidth = w->get_value_as_int();
    pp->raw.rcHeight = h->get_value_as_int();

    if (pedited) {
        pedited->raw.rawCrop = !get_inconsistent();
        pedited->raw.rcWidth = wSet;
        pedited->raw.rcHeight = hSet;
        pedited->raw.rcX = xSet;
        pedited->raw.rcY = ySet;
    }

}

void RawCrop::trim (ProcParams* pp, int ow, int oh)
{

    int xmin = pp->raw.rcX;
    int ymin = pp->raw.rcY;

    if (xmin > ow || ymin > oh) {
        // the crop is completely out of the image, so we disable the crop
        pp->raw.rawCrop = false;
        // and we set the values to the defaults
        pp->raw.rcX = 0;
        pp->raw.rcY = 0;
        pp->raw.rcWidth = ow;
        pp->raw.rcHeight = oh;
    } else {
        if ((xmin + pp->raw.rcWidth) > ow) {
            // crop overflow in the width dimension ; we trim it
            pp->raw.rcWidth = ow - xmin;
        }

        if ((ymin + pp->raw.rcHeight) > oh) {
            // crop overflow in the height dimension ; we trim it
            pp->raw.rcHeight = oh - ymin;
        }
    }
}

void RawCrop::enabledChanged ()
{

    auto m = ProcEventMapper::getInstance();
    m->remapEvent(EvRawCropEnabled, edit->get_active() ? RAWCROP|M_MODE_RAWCROPADJUST : RAWCROP, "TP_RAW_CROP_LABEL");

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvRawCropEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvRawCropEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvRawCropEnabled, M("GENERAL_DISABLED"));
        }
    }
}

/*
 * Set the maximum dimensions of the image. This method can be called with wrong values, then
 * called with the good ones !?
 */
void RawCrop::setDimensions (int mw, int mh)
{

    ConnectionBlocker xBlocker(xconn);
    ConnectionBlocker yBlocker(yconn);
    ConnectionBlocker wBlocker(wconn);
    ConnectionBlocker hBlocker(hconn);

    disableListener ();

    maxw = mw;
    maxh = mh;

    w->set_range (1, maxw);
    h->set_range (1, maxh);
    x->set_range (0, maxw - 1);
    y->set_range (0, maxh - 1);

    enableListener ();
}

void RawCrop::sizeChanged (int x, int y, int ow, int oh)
{
    struct Params {
        RawCrop* rawCrop;
        int x;
        int y;
    };

    Params* const params = new Params{
        this,
        x,
        y
    };

    const auto func = [](gpointer data) -> gboolean {
        Params* const params = static_cast<Params*>(data);
        params->rawCrop->setDimensions(params->x, params->y);
        delete params;
        return FALSE;
    };

    idle_register.add(func, params);
}


void RawCrop::spinChanged (MySpinButton *spinButton)
{

    if (spinButton == x) {
        xSet = true;
    } else if (spinButton == y) {
        ySet = true;
    } else if (spinButton == w) {
        wSet = true;
    } else if (spinButton == h) {
        hSet = true;
    }

    notifyListener (x->get_value_as_int(), y->get_value_as_int(), w->get_value_as_int(), h->get_value_as_int());
}

void RawCrop::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
}

void RawCrop::setEditProvider (EditDataProvider* provider)
{
    EditSubscriber::setEditProvider(provider);
}

void RawCrop::editToggled ()
{
    if (listener) {
        auto m = ProcEventMapper::getInstance();
        m->remapEvent(EvRawCropDims, edit->get_active() ? RAWCROP|M_MODE_RAWCROPADJUST : RAWCROP, "HISTORY_MSG_RAW_CROP_DIM");
        m->remapEvent(EvRawCropDimsOPA, edit->get_active() ? MINUPDATE|M_MODE_RAWCROPADJUST : MINUPDATE, "");

        if (edit->get_active()) {
            // this will refresh the preview with RawCrop disabled
            listener->refreshPreview(EvRawCropDimsOPA);
            subscribe();
        } else {
            // this will refresh the preview with RawCrop enabled, if active
            unsubscribe();
            listener->refreshPreview(EvRawCropDims);
        }
    }
}

CursorShape RawCrop::getCursor(const int objectID)
{
    // We don't care about objectID, its value is -2, i.e. no Mouse Over objects mode
    if (litArea == LitArea::topLeft) {
        return CSResizeTopLeft;
    } else if (litArea == LitArea::topRight) {
        return CSResizeTopRight;
    } else if (litArea == LitArea::bottomLeft) {
        return CSResizeBottomLeft;
    } else if (litArea == LitArea::bottomRight) {
        return CSResizeBottomRight;
    } else if (litArea == LitArea::left || litArea == LitArea::right) {
        return CSResizeWidth;
    } else if (litArea == LitArea::top || litArea == LitArea::bottom) {
        return CSResizeHeight;
    } else if (litArea == LitArea::center) {
        return CSMove2D;
    } else {
        return CSArrow;
    }
}

void RawCrop::updateGeometry(const int x, const int y, const int width, const int height, const int fullWidth, const int fullHeight)
{
    EditDataProvider* dataProvider = getEditProvider();

    if (!dataProvider) {
        return;
    }

    int imW=0;
    int imH=0;
    if (fullWidth != -1 && fullHeight != -1) {
        imW = fullWidth;
        imH = fullHeight;
    } else {
        dataProvider->getImageSize(imW, imH);
        if (!imW || !imH) {
            return;
        }
    }

    int x2 = x;
    int y2 = y;
    int smallCornerW = width / 8;
    int smallCornerH = height / 8;
    int bigCornerW = width / 6;
    int bigCornerH = height / 6;
    int cornerW;
    int cornerH;

    cornerW = litArea == LitArea::topLeft ? bigCornerW : smallCornerW;
    cornerH = litArea == LitArea::topLeft ? bigCornerH : smallCornerH;
    topLeft.points.at(0) = rtengine::Coord(x2, y2 + cornerH);
    topLeft.points.at(1) = rtengine::Coord(x2, y2);
    topLeft.points.at(2) = rtengine::Coord(x2 + cornerW, y2);

    x2 = x + width;
    cornerW = litArea == LitArea::topRight ? bigCornerW : smallCornerW;
    cornerH = litArea == LitArea::topRight ? bigCornerH : smallCornerH;
    topRight.points.at(0) = rtengine::Coord(x2, y2 + cornerH);
    topRight.points.at(1) = rtengine::Coord(x2, y2);
    topRight.points.at(2) = rtengine::Coord(x2 - cornerW, y2);


    y2 = y + height;
    cornerW = litArea == LitArea::bottomRight ? bigCornerW : smallCornerW;
    cornerH = litArea == LitArea::bottomRight ? bigCornerH : smallCornerH;
    bottomRight.points.at(0) = rtengine::Coord(x2, y2 - cornerH);
    bottomRight.points.at(1) = rtengine::Coord(x2, y2);
    bottomRight.points.at(2) = rtengine::Coord(x2 - cornerW, y2);

    x2 = x;
    cornerW = litArea == LitArea::bottomLeft ? bigCornerW : smallCornerW;
    cornerH = litArea == LitArea::bottomLeft ? bigCornerH : smallCornerH;
    bottomLeft.points.at(0) = rtengine::Coord(x2, y2 - cornerH);
    bottomLeft.points.at(1) = rtengine::Coord(x2, y2);
    bottomLeft.points.at(2) = rtengine::Coord(x2 + cornerW, y2);

    cornerW = bigCornerW;
    cornerH = bigCornerH;

    switch (litArea) {
    case LitArea::top:
        side.points.at(0) = rtengine::Coord(x, y + cornerH);
        side.points.at(1) = rtengine::Coord(x, y);
        side.points.at(2) = rtengine::Coord(x + width, y);
        side.points.at(3) = rtengine::Coord(x + width, y + cornerH);
        break;
    case LitArea::left:
        side.points.at(0) = rtengine::Coord(x + cornerW, y);
        side.points.at(1) = rtengine::Coord(x, y);
        side.points.at(2) = rtengine::Coord(x, y + height);
        side.points.at(3) = rtengine::Coord(x + cornerW, y + height);
        break;
    case LitArea::right:
        side.points.at(0) = rtengine::Coord(x + width - cornerW, y);
        side.points.at(1) = rtengine::Coord(x + width, y);
        side.points.at(2) = rtengine::Coord(x + width, y + height);
        side.points.at(3) = rtengine::Coord(x + width - cornerW, y + height);
        break;
    case LitArea::bottom:
        side.points.at(0) = rtengine::Coord(x, y + height - cornerH);
        side.points.at(1) = rtengine::Coord(x, y + height);
        side.points.at(2) = rtengine::Coord(x + width, y + height);
        side.points.at(3) = rtengine::Coord(x + width, y + height -cornerH);
        break;
    default:
        break;
    }

    centerCircle.center.set(x + width / 2, y + height / 2);
    centerCircle.radius = litArea == LitArea::center ? 10 : 6;  // same value as in ctor
}

void RawCrop::updateState()
{
    Geometry::State newState = EditSubscriber::action == ES_ACTION_DRAGGING ? Geometry::DRAGGED : Geometry::PRELIGHT;
    for (auto currGeometry : EditSubscriber::visibleGeometry) {
        currGeometry->setVisible(true);
        currGeometry->state = Geometry::State::NORMAL;
    }

    switch (litArea) {
    case LitArea::none:
        side.setVisible(false);
        break;
    case LitArea::topLeft:
        side.setVisible(false);
        topLeft.state = newState;
        break;
    case LitArea::top:
        topLeft.setVisible(false);
        topRight.setVisible(false);
        side.state = newState;
        break;
    case LitArea::topRight:
        side.setVisible(false);
        topRight.state = newState;
        break;
    case LitArea::left:
        topLeft.setVisible(false);
        bottomLeft.setVisible(false);
        side.state = newState;
        break;
    case LitArea::center:
        side.setVisible(false);
        centerCircle.state = newState;
        break;
    case LitArea::right:
        topRight.setVisible(false);
        bottomRight.setVisible(false);
        side.state = newState;
        break;
    case LitArea::bottomLeft:
        side.setVisible(false);
        bottomLeft.state = newState;
        break;
    case LitArea::bottom:
        bottomLeft.setVisible(false);
        bottomRight.setVisible(false);
        side.state = newState;
        break;
    case LitArea::bottomRight:
        side.setVisible(false);
        bottomRight.state = newState;
        break;
    }
}

bool RawCrop::inArea(const int x1, const int y1, const int x2, const int y2, const rtengine::Coord &mousePos)
{
    return mousePos.x >= x1 && mousePos.x <= x2 && mousePos.y >= y1 && mousePos.y <= y2;
}

bool RawCrop::mouseOver(const int modifierKey)
{
    EditDataProvider* editProvider = getEditProvider();
    LitArea lastLitArea = litArea;


    int x_ = x->get_value_as_int();
    int y_ = y->get_value_as_int();
    int w_ = w->get_value_as_int();
    int h_ = h->get_value_as_int();
    int imgWidth = -1, imgHeight = -1;

    if (!(modifierKey & GDK_CONTROL_MASK)) {
        litArea = LitArea::none;
        if (editProvider) {
            editProvider->getImageSize(imgWidth, imgHeight);
        }
    } else if (editProvider) {
        editProvider->getImageSize(imgWidth, imgHeight);
        int cornerW = w_ / 5;
        int cornerH = h_ / 5;

        if (inArea(0, 0, x_ + cornerW, y_ + cornerH, editProvider->posImage)) {
            litArea = LitArea::topLeft;
        } else if (inArea(x_ + cornerW, 0, x_ + w_ - cornerW, y_ + cornerH, editProvider->posImage)) {
            litArea = LitArea::top;
        } else if (inArea(x_ + w_ - cornerW, 0, imgWidth, y_ + cornerH, editProvider->posImage)) {
            litArea = LitArea::topRight;
        } else if (inArea(0, y_ + cornerH, x_ + cornerW, y_ + h_ - cornerH, editProvider->posImage)) {
            litArea = LitArea::left;
        } else if (inArea(x_ + cornerW, y_ + cornerH, x_ + w_ - cornerW, y_ + h_ - cornerH, editProvider->posImage)) {
            litArea = LitArea::center;
        } else if (inArea(x_ + w_ - cornerW, y_ + cornerH, imgWidth, y_ + h_ - cornerH, editProvider->posImage)) {
            litArea = LitArea::right;
        } else if (inArea(0, y_ - cornerH, x_ + cornerW, imgHeight, editProvider->posImage)) {
            litArea = LitArea::bottomLeft;
        } else if (inArea(x_ + cornerW, y_ - cornerH, x_ + w_ - cornerW, imgHeight, editProvider->posImage)) {
            litArea = LitArea::bottom;
        } else if (inArea(x_ + w_ - cornerW, y_ - cornerH, imgWidth, imgHeight, editProvider->posImage)) {
            litArea = LitArea::bottomRight;
        } else {
            litArea = LitArea::none;
        }
    }

    if (editProvider && lastLitArea != litArea) {
        switch (litArea) {
        case LitArea::topLeft:
            topLeft.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::top:
            side.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::topRight:
            topRight.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::left:
            side.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::center:
            centerCircle.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::right:
            side.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::bottomLeft:
            bottomLeft.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::bottom:
            side.state = Geometry::State::PRELIGHT;
            break;
        case LitArea::bottomRight:
            bottomRight.state = Geometry::State::PRELIGHT;
            break;
        default:
            break;
        }
        updateState();
        updateGeometry(x_, y_, w_, h_, imgWidth, imgHeight);
        return true;
    }

    return false;
}

bool RawCrop::button1Pressed(const int modifierKey)
{
    if (litArea == LitArea::none || modifierKey != GDK_CONTROL_MASK) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    // button press is valid (no modifier key)
    int imW, imH;
    provider->getImageSize(imW, imH);

    oldX[0] = x->get_value_as_int();
    oldY[0] = y->get_value_as_int();
    oldX[1] = rtengine::max<int>(w->get_value_as_int() - 1, 0) + oldX[0];
    oldY[1] = rtengine::max<int>(h->get_value_as_int() - 1, 0) + oldY[0];

    EditSubscriber::action = ES_ACTION_DRAGGING;
    updateState();

    return true;
}

bool RawCrop::button1Released()
{
    EditSubscriber::action = ES_ACTION_NONE;
    centerCircle.radius = 6;  // same value as in ctor
    litArea = LitArea::none;
    updateState();
    return true;
}

bool RawCrop::drag1(const int modifierKey)
{
    ConnectionBlocker xBlocker(xconn);
    ConnectionBlocker yBlocker(yconn);
    ConnectionBlocker wBlocker(wconn);
    ConnectionBlocker hBlocker(hconn);

    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize(imW, imH);
    int X[2];
    int Y[2];

    X[0] = oldX[0];
    Y[0] = oldY[0];
    X[1] = oldX[1];
    Y[1] = oldY[1];

    switch (litArea) {
    case LitArea::topLeft:
        X[0] = rtengine::LIM<int>(oldX[0] + provider->deltaImage.x, 0, X[1]);
        Y[0] = rtengine::LIM<int>(oldY[0] + provider->deltaImage.y, 0, Y[1]);
        break;
    case LitArea::top:
        Y[0] = rtengine::LIM<int>(oldY[0] + provider->deltaImage.y, 0, Y[1]);
        break;
    case LitArea::topRight:
        X[1] = rtengine::LIM<int>(oldX[1] + provider->deltaImage.x, X[0], imW - 1);
        Y[0] = rtengine::LIM<int>(oldY[0] + provider->deltaImage.y, 0, Y[1]);
        break;
    case LitArea::left:
        X[0] = rtengine::LIM<int>(oldX[0] + provider->deltaImage.x, 0, X[1]);
        break;
    case LitArea::center:
        X[0] = oldX[0] + provider->deltaImage.x;
        Y[0] = oldY[0] + provider->deltaImage.y;
        X[1] = oldX[1] + provider->deltaImage.x;
        Y[1] = oldY[1] + provider->deltaImage.y;
        if (X[0] < 0) {
            X[1] -= X[0];
            X[0] = 0;
        } else if (X[1] > (imW - 1)) {
            X[0] -= X[1] - imW - 1;
            X[1] = imW - 1;
        }
        if (Y[0] < 0) {
            Y[1] -= Y[0];
            Y[0] = 0;
        } else if (Y[1] > (imH - 1)) {
            Y[0] -= Y[1] - imH - 1;
            Y[1] = imH - 1;
        }
        break;
    case LitArea::right:
        X[1] = rtengine::LIM<int>(oldX[1] + provider->deltaImage.x, X[0], imW - 1);
        break;
    case LitArea::bottomLeft:
        X[0] = rtengine::LIM<int>(oldX[0] + provider->deltaImage.x, 0, X[1]);
        Y[1] = rtengine::LIM<int>(oldY[1] + provider->deltaImage.y, Y[0], imH - 1);
        break;
    case LitArea::bottom:
        Y[1] = rtengine::LIM<int>(oldY[1] + provider->deltaImage.y, Y[0], imH - 1);
        break;
    case LitArea::bottomRight:
        X[1] = rtengine::LIM<int>(oldX[1] + provider->deltaImage.x, X[0], imW - 1);
        Y[1] = rtengine::LIM<int>(oldY[1] + provider->deltaImage.y, Y[0], imH - 1);
        break;
    default:
        break;
    }

    if (X[0] != oldX[0] || X[1] != oldX[1] || Y[0] != oldY[0] || Y[1] != oldY[1]) {
        int nx = X[0];
        int ny = Y[0];
        int nw = X[1] - X[0] + 1;
        int nh = Y[1] - Y[0] + 1;
        x->set_value(nx);
        y->set_value(ny);
        w->set_value(nw);
        h->set_value(nh);
        updateGeometry (nx, ny, nw, nh, imW, imH);
        notifyListener (nx, ny, nw, nh);
        return true;
    }

    return false;
}

void RawCrop::notifyListener(int nx, int ny, int nw, int nh)
{
    if (listener) {
        auto m = ProcEventMapper::getInstance();
        m->remapEvent(EvRawCropDims, edit->get_active() ? RAWCROP|M_MODE_RAWCROPADJUST : RAWCROP, "HISTORY_MSG_RAW_CROP_DIM");
        m->remapEvent(EvRawCropDimsOPA, edit->get_active() ? MINUPDATE|M_MODE_RAWCROPADJUST : MINUPDATE, "");

        listener->panelChanged (edit->get_active() ? EvRawCropDimsOPA : EvRawCropDims,
                                Glib::ustring::compose ("%1=%2, %3=%4\n%5=%6, %7=%8",
                                        M("TP_CROP_X"),
                                        x->get_value_as_int(),
                                        M("TP_CROP_Y"),
                                        y->get_value_as_int(),
                                        M("TP_CROP_W"),
                                        w->get_value_as_int(),
                                        M("TP_CROP_H"),
                                        h->get_value_as_int())
                                );
    }
}
void RawCrop::switchOffEditMode ()
{
    if (edit->get_active()) {
        // switching off the toggle button
        edit->set_active(false);
    }
}

rtengine::ProcEvent RawCrop::getCropEnabledEvent() const
{
    return EvRawCropEnabled;
}

rtengine::ProcEvent RawCrop::getCropDimsEvent() const
{
    return EvRawCropDims;
}
rtengine::ProcEvent RawCrop::getCropDimsOPAEvent() const
{
    return EvRawCropDimsOPA;
}

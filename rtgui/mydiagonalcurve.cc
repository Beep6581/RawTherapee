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
#include "mydiagonalcurve.h"
#include "../rtengine/curves.h"
#include <cstring>
#include <gdkmm/types.h>

MyDiagonalCurve::MyDiagonalCurve () : activeParam(-1), bghistvalid(false)
{

    graphW = get_allocation().get_width() - RADIUS * 2;
    graphH = get_allocation().get_height() - RADIUS * 2;
    prevGraphW = graphW;
    prevGraphH = graphH;
    grab_point = -1;
    lit_point = -1;
    buttonPressed = false;

    bghist = new unsigned int[256];

    editedPos.resize(2);
    editedPos.at(0) = editedPos.at(1) = 0.0;

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->add_class(GTK_STYLE_CLASS_FRAME);
    //style->add_class(GTK_STYLE_CLASS_NOTEBOOK);
    //style->add_class(GTK_STYLE_CLASS_FLAT);

    signal_event().connect( sigc::mem_fun(*this, &MyDiagonalCurve::handleEvents) );

    curve.type = DCT_Spline;

    curve.x.push_back(0.);
    curve.y.push_back(0.);
    curve.x.push_back(1.);
    curve.y.push_back(1.);
}

MyDiagonalCurve::~MyDiagonalCurve ()
{

    delete [] bghist;
}

std::vector<double> MyDiagonalCurve::get_vector (int veclen)
{

    std::vector<double> vector;
    vector.resize (veclen);

    if (curve.type != DCT_Parametric) {
        // count active points:
        double prev = - 1.0;
        int active = 0;
        int firstact = -1;

        for (int i = 0; i < (int)curve.x.size(); ++i)
            if (curve.x.at(i) > prev) {
                if (firstact < 0) {
                    firstact = i;
                }

                prev = curve.x.at(i);
                ++active;
            }

        // handle degenerate case:
        if (active < 2) {
            double ry;

            if (active > 0) {
                ry = curve.y.at(firstact);
            } else {
                ry = 0.0;
            }

            if (ry < 0.0) {
                ry = 0.0;
            }

            if (ry > 1.0) {
                ry = 1.0;
            }

            for (int x = 0; x < veclen; ++x) {
                vector.at(x) = ry;
            }

            return vector;
        }
    }

    // calculate remaining points
    std::vector<double> curveDescr = getPoints ();
    rtengine::DiagonalCurve rtcurve(curveDescr, veclen * 1.2);
    std::vector<double> t;
    t.resize (veclen);

    for (int i = 0; i < veclen; i++) {
        t[i] = (double) i / (veclen - 1.0);
    }

    rtcurve.getVal (t, vector);
    return vector;
}

void MyDiagonalCurve::get_LUT (LUTf &lut)
{

    int size = lut.getSize();

    if (curve.type != DCT_Parametric) {
        // count active points:
        double prev = - 1.0;
        int active = 0;
        int firstact = -1;

        for (int i = 0; i < (int)curve.x.size(); ++i)
            if (curve.x.at(i) > prev) {
                if (firstact < 0) {
                    firstact = i;
                }

                prev = curve.x.at(i);
                ++active;
            }

        // handle degenerate case:
        if (active < 2) {
            double ry;

            if (active > 0) {
                ry = curve.y.at(firstact);
            } else {
                ry = 0.0;
            }

            if (ry < 0.0) {
                ry = 0.0;
            }

            if (ry > 1.0) {
                ry = 1.0;
            }

            for (int x = 0; x < size; ++x) {
                lut[x] = ry;
            }

            return;
        }
    }

    // calculate remaining points
    std::vector<double> curveDescr = getPoints ();
    rtengine::DiagonalCurve rtcurve(curveDescr, lut.getUpperBound() * 1.2);
    double t;
    double maxVal = double(lut.getUpperBound());

    for (int i = 0; i < size; i++) {
        t = double(i) / maxVal;
        lut[i] = rtcurve.getVal (t);
    }

    return;
}

void MyDiagonalCurve::interpolate ()
{

    prevGraphW = graphW;
    prevGraphH = graphH;
    int nbPoints = rtengine::max(graphW - 2, 201);
    point(nbPoints);
    get_LUT (point);
    upoint.reset();
    lpoint.reset ();

    if (curve.type == DCT_Parametric && activeParam > 0) {
        double tmp = curve.x.at(activeParam - 1);

        if (activeParam >= 4) {
            upoint(nbPoints);
            lpoint(nbPoints);
            curve.x.at(activeParam - 1) = 100;
            get_LUT(upoint);
            curve.x.at(activeParam - 1) = -100;
            get_LUT (lpoint);
            curve.x.at(activeParam - 1) = tmp;
        }
    }

    curveIsDirty = false;
}

void MyDiagonalCurve::draw (int handle)
{
    if (!isDirty()) {
        return;
    }

    Glib::RefPtr<Gdk::Window > win = get_window();

    if (!surfaceCreated() || !win) {
        return;
    }

    // re-calculate curve if dimensions changed
    int currPointSize = point.getUpperBound();

    if (curveIsDirty || /*prevGraphW != graphW || prevGraphH != graphH ||*/ (currPointSize == GRAPH_SIZE && (graphW - 3 > GRAPH_SIZE)) || (currPointSize > GRAPH_SIZE && (graphW - 2 <= GRAPH_SIZE || graphW - 3 != currPointSize))) {
        interpolate ();
    }

    currPointSize = point.getUpperBound();

    Gtk::StateFlags state = !is_sensitive() ? Gtk::STATE_FLAG_INSENSITIVE : Gtk::STATE_FLAG_NORMAL;

    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    Cairo::RefPtr<Cairo::Context> cr = getContext();
    cr->set_line_cap(Cairo::LINE_CAP_SQUARE);

    // clear background
    style->render_background(cr, 0., 0., double(getWidth()), double(getHeight()));

    Gdk::RGBA c;

    // histogram in the background
    if (bghistvalid) {
        // find highest bin
        unsigned int valMax = 0;

        for (int i = 0; i < 256; i++)
            if (bghist[i] > valMax) {
                valMax = bghist[i];
            }

        // draw histogram
        cr->set_line_width (1.0);
        double stepSize = (graphW - 3) / 255.0;
        cr->move_to ( double(graphX + 1), double(graphY - 1) );
        c = style->get_background_color(Gtk::STATE_FLAG_INSENSITIVE);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());

        for (int i = 0; i < 256; i++) {
            double val = double(bghist[i]) * double(graphH - 2) / double(valMax);
            /*
            if (val>graphH-2)
                val = graphH-2;
            */
            //if (i>0)
            cr->line_to (double(graphX) + 1.5 + double(i)*stepSize, double(graphY - 1) - val);
        }

        cr->line_to (double(graphX) + 1.5 + 255.*stepSize, double(graphY - 1));
        cr->close_path();
        cr->fill ();
    }

    // draw the grid lines:
    cr->set_line_width (1.0);
    c = style->get_border_color(state);
    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    cr->set_antialias (Cairo::ANTIALIAS_NONE);

    for (int i = 0; i <= 10; i++) {
        // horizontal lines
        cr->move_to     (double(graphX) + 0.5                                 , double(graphY) - max(0.5, double(graphH * i / 10) - 0.5));
        cr->rel_line_to (double(graphW - 1)                                   , 0.);
        // vertical lines
        cr->move_to     (double(graphX) + max(0.5, double(graphW * i / 10) - 0.5), double(graphY));
        cr->rel_line_to (0.                                                 , double(-graphH + 1));
    }

    cr->stroke ();

    // draw f(x)=x line
    if (snapToElmt == -2) {
        cr->set_source_rgb (1.0, 0.0, 0.0);
    } else {
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    }

    std::valarray<double> ds (1);
    ds[0] = 4;
    cr->set_dash (ds, 0);
    cr->move_to (double(graphX) + 1.5, double(graphY) - 1.5);
    cr->rel_line_to (double(graphW - 3), double(-graphH + 3));
    cr->stroke ();
    cr->unset_dash ();

    cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
    cr->set_line_width (1.0);

    // draw upper and lower bounds
    float graphH_ = float(graphH - 3);
    float graphX_ = float(graphX) + 1.5;
    float graphY_ = float(graphY) - 1.5;

    if (curve.type == DCT_Parametric && activeParam > 0 && lpoint.getUpperBound() > 1 && upoint.getUpperBound() > 1) {
        cr->set_source_rgba (1.0, 1.0, 1.0, 0.1);
        cr->move_to (graphX_, getVal(upoint, 0) * -graphH_ + graphY_);

        for (int i = 1; i < graphW - 2; ++i) {
            cr->line_to (float(i) + graphX_, getVal(upoint, i) * -graphH_ + graphY_);
        }

        for (int i = graphW - 3; i >= 0; --i) {
            cr->line_to (float(i) + graphX_, getVal(lpoint, i) * -graphH_ + graphY_);
        }

        cr->fill ();
    }

    // draw the pipette values
    if (pipetteR > -1.f || pipetteG > -1.f || pipetteB > -1.f) {
        int n = 0;

        if (pipetteR > -1.f) {
            ++n;
        }

        if (pipetteG > -1.f) {
            ++n;
        }

        if (pipetteB > -1.f) {
            ++n;
        }

        if (n > 1) {
            if (pipetteR > -1.f) {
                cr->set_source_rgba (1., 0., 0., 0.5); // WARNING: assuming that red values are stored in pipetteR, which might not be the case!
                cr->move_to (double(graphX) + 1.5 + double(graphW - 3)*pipetteR, double(graphY) - 1.5);
                cr->rel_line_to (0, double(-graphH + 3));
                cr->stroke ();
            }

            if (pipetteG > -1.f) {
                cr->set_source_rgba (0., 1., 0., 0.5); // WARNING: assuming that green values are stored in pipetteG, which might not be the case!
                cr->move_to (double(graphX) + 1.5 + double(graphW - 3)*pipetteG, double(graphY) - 1.5);
                cr->rel_line_to (0, double(-graphH + 3));
                cr->stroke ();
            }

            if (pipetteB > -1.f) {
                cr->set_source_rgba (0., 0., 1., 0.5); // WARNING: assuming that blue values are stored in pipetteB, which might not be the case!
                cr->move_to (double(graphX) + 1.5 + double(graphW - 3)*pipetteB, double(graphY) - 1.5);
                cr->rel_line_to (0, double(-graphH + 3));
                cr->stroke ();
            }
        }

        if (pipetteVal > -1.f) {
            cr->set_line_width (2.);
            c = style->get_color (state);
            cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            cr->move_to (double(graphX) + 1.5 + double(graphW - 3)*pipetteVal, double(graphY) - 1.5);
            cr->rel_line_to (0, double(-graphH + 3));
            cr->stroke ();
            cr->set_line_width (1.);
        }
    }

    c = style->get_color (state);

    // draw the cage of the NURBS curve
    if (curve.type == DCT_NURBS) {
        unsigned int nbPoints;
        std::valarray<double> ch_ds (1);
        ch_ds[0] = 2;
        cr->set_dash (ch_ds, 0);
        cr->set_line_width (0.75);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        std::vector<double> points = getPoints();
        nbPoints = ((int)points.size() - 1) / 2;

        for (unsigned int i = 1; i < nbPoints; i++) {
            int pos = i * 2 + 1;

            double x1 = double(graphX) + 1.5 + double(graphW - 3) * points[pos - 2]; // project (curve.at(i), 0, 1, graphW);
            double y1 = double(graphY) - 1.5 - double(graphH - 3) * points[pos - 1]; // project (curve.y.at(i)i], 0, 1, graphH);
            double x2 = double(graphX) + 0.5 + double(graphW - 3) * points[pos]; // project (curve.at(i), 0, 1, graphW);
            double y2 = double(graphY) - 1.5 - double(graphH - 3) * points[pos + 1]; // project (curve.y.at(i), 0, 1, graphH);

            // set the color of the line when the point is snapped to the cage
            if (curve.x.size() == nbPoints && snapToElmt >= 1000 && ((i == (snapToElmt - 1000)) || (i == (snapToElmt - 999)))) {
                cr->set_source_rgb (1.0, 0.0, 0.0);
            } else {
                cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
            }

            cr->move_to (x1, y1);
            cr->line_to (x2, y2);
            cr->stroke ();
        }

        cr->unset_dash ();
        cr->set_line_width (1.0);
    }

    // draw curve
    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
    cr->move_to (graphX_, getVal(point, 0) * -graphH_ + graphY_);

    for (int i = 1; i < graphW - 2; ++i) {
        cr->line_to (float(i) + graphX_, getVal(point, i) * -graphH_ + graphY_);
    }

    cr->stroke ();

    // draw the left colored bar
    if (leftBar) {
        // first the background
        int bWidth = getBarWidth();
        BackBuffer *bb = this;
        leftBar->setDrawRectangle(win, 1, graphY - graphH + 1, bWidth - 2, graphH - 2);
        leftBar->expose(bb);

        // now the border
        c = style->get_border_color(state);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->rectangle(0.5, graphY - graphH + 0.5, bWidth - 1, graphH - 1);
        cr->stroke();
    }

    // draw the bottom colored bar
    if (bottomBar) {
        // first the background
        int bWidth = getBarWidth();
        BackBuffer *bb = this;
        bottomBar->setDrawRectangle(win, graphX + 1, graphY + CBAR_MARGIN + 1, graphW - 2, bWidth - 2);
        bottomBar->expose(bb);

        // now the border
        c = style->get_border_color (state);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->rectangle(graphX + 0.5, graphY + CBAR_MARGIN + 0.5, graphW - 1, bWidth - 1 );
        cr->stroke();
    }

    // draw bullets
    if (curve.type != DCT_Parametric) {
        c = style->get_color (state);

        for (int i = 0; i < (int)curve.x.size(); ++i) {
            if (curve.x.at(i) == -1) {
                continue;
            }

            if (snapToElmt >= 1000) {
                int pt = snapToElmt - 1000;

                if (i >= (pt - 1) && i <= (pt + 1)) {
                    cr->set_source_rgb(1.0, 0.0, 0.0);
                } else {
                    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
                }
            } else {
                if (i == handle || i == snapToElmt || i == edited_point) {
                    cr->set_source_rgb (1.0, 0.0, 0.0);
                } else {
                    cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
                }
            }

            double x = double(graphX + 1) + double((graphW - 2) * curve.x.at(i)); // project (curve.x.at(i), 0, 1, graphW);
            double y = double(graphY - 1) - double((graphH - 2) * curve.y.at(i)); // project (curve.y.at(i), 0, 1, graphH);

            cr->arc (x, y, RADIUS + 0.5, 0, 2 * M_PI);
            cr->fill ();

            if (i == edited_point) {
                cr->set_line_width(2.);
                cr->arc (x, y, RADIUS + 3.5, 0, 2 * M_PI);
                cr->stroke();
                cr->set_line_width(1.);
            }

        }
    }

    setDirty(false);
    queue_draw();
}

bool MyDiagonalCurve::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{
    Glib::RefPtr<Gdk::Window> win = get_window();
    Gtk::Allocation allocation = get_allocation();
    allocation.set_x(0);
    allocation.set_y(0);

    // setDrawRectangle will allocate the backbuffer Surface
    if (setDrawRectangle(win, allocation)) {
        setDirty(true);

        if (prevGraphW > GRAPH_SIZE || graphW > GRAPH_SIZE) {
            curveIsDirty = true;
        }
    }

    draw (lit_point);
    copySurface(cr);
    return false;
}

/*void MyDiagonalCurve::graphSizeRequest (Gtk::Requisition* req) {
    req->width = getGraphMinSize();
    // The real height request should take care of the presence of the vertical
    // scroll bar and its width
    req->height = sized ? getGraphMinSize() : get_allocation().get_width();
}*/

bool MyDiagonalCurve::handleEvents (GdkEvent* event)
{

    CursorShape new_type = cursor_type;

    bool retval = false;
    int num = (int)curve.x.size();

    /* graphW and graphH are the size of the graph */
    calcDimensions();

    if ((graphW < 0) || (graphH < 0)) {
        return false;
    }

    double minDistanceX = double(MIN_DISTANCE) / double(graphW - 1);
    double minDistanceY = double(MIN_DISTANCE) / double(graphH - 1);

    switch (event->type) {
    case Gdk::BUTTON_PRESS:
        snapToElmt = -100;

        if (curve.type != DCT_Parametric) {
            if (edited_point == -1) {
                if (event->button.button == 1) {
                    std::vector<double>::iterator itx, ity;
                    buttonPressed = true;
                    add_modal_grab ();

                    // get the pointer position
                    getCursorPosition(Gdk::EventType(event->type), event->motion.is_hint != 0, int(event->button.x), int(event->button.y), Gdk::ModifierType(event->button.state));
                    findClosestPoint();

                    new_type = CSMove;

                    if (distanceX > minDistanceX) {
                        if (mod_type & GDK_CONTROL_MASK) {
                            clampedY = point.getVal01(clampedX);
                        }

                        /* insert a new control point */
                        if (num > 0) {
                            if (clampedX > curve.x.at(closest_point)) {
                                ++closest_point;
                            }
                        }

                        itx = curve.x.begin();
                        ity = curve.y.begin();

                        for (int i = 0; i < closest_point; i++) {
                            itx++;
                            ity++;
                        }

                        curve.x.insert (itx, 0);
                        curve.y.insert (ity, 0);
                        num++;

                        // the graph is refreshed only if a new point is created
                        curve.x.at(closest_point) = clampedX;
                        curve.y.at(closest_point) = clampedY;

                        curveIsDirty = true;
                        setDirty(true);
                        draw (closest_point);
                        notifyListener ();
                    }

                    grab_point = closest_point;
                    lit_point = closest_point;
                    ugpX = curve.x.at(closest_point);
                    ugpY = curve.y.at(closest_point);
                } else if (event->button.button == 3) {
                    if (lit_point > -1 && grab_point == -1) {
                        if (!coordinateAdjuster->is_visible()) {
                            coordinateAdjuster->showMe(this);
                        }

                        edited_point = lit_point;
                        std::vector<CoordinateAdjuster::Boundaries> newBoundaries(2);
                        unsigned int size = curve.x.size();

                        if      (edited_point == 0)      {
                            newBoundaries.at(0).minVal = 0.;
                            newBoundaries.at(0).maxVal = curve.x.at(1);
                        } else if (edited_point == size - 1) {
                            newBoundaries.at(0).minVal = curve.x.at(edited_point - 1);
                            newBoundaries.at(0).maxVal = 1.;
                        } else if (curve.x.size() > 2)     {
                            newBoundaries.at(0).minVal = curve.x.at(edited_point - 1);
                            newBoundaries.at(0).maxVal = curve.x.at(edited_point + 1);
                        }

                        newBoundaries.at(1).minVal = 0.;
                        newBoundaries.at(1).maxVal = 1.;
                        editedPos.at(0) = curve.x.at(edited_point);
                        editedPos.at(1) = curve.y.at(edited_point);
                        coordinateAdjuster->setPos(editedPos);
                        coordinateAdjuster->startNumericalAdjustment(newBoundaries);
                        setDirty(true);
                        draw (lit_point);
                        new_type = CSArrow;
                        retval = true;
                    }
                }

                if (buttonPressed) {
                    retval = true;
                }
            } else { // if (edited_point > -1)
                if (event->button.button == 3) {
                    // do we edit another point?
                    if (edited_point > -1 && grab_point == -1) {
                        /*  get the pointer position  */
                        getCursorPosition(Gdk::EventType(event->type), event->motion.is_hint != 0, int(event->button.x), int(event->button.y), Gdk::ModifierType(event->button.state));
                        findClosestPoint();

                        if (cursorX >= 0 && cursorX <= graphW && cursorY >= 0 && cursorY <= graphH) {
                            if (distanceX <= minDistanceX) {
                                // the cursor is close to an existing point
                                lit_point = closest_point;

                                if (lit_point != edited_point) {
                                    edited_point = lit_point;
                                    curveIsDirty = true;
                                    setDirty(true);
                                    draw (lit_point);
                                    std::vector<CoordinateAdjuster::Boundaries> newBoundaries;
                                    newBoundaries.resize(2);
                                    unsigned int size = curve.x.size();

                                    if      (edited_point == 0)      {
                                        newBoundaries.at(0).minVal = 0.;
                                        newBoundaries.at(0).maxVal = curve.x.at(1);
                                    } else if (edited_point == size - 1) {
                                        newBoundaries.at(0).minVal = curve.x.at(edited_point - 1);
                                        newBoundaries.at(0).maxVal = 1.;
                                    } else if (curve.x.size() > 2)     {
                                        newBoundaries.at(0).minVal = curve.x.at(edited_point - 1);
                                        newBoundaries.at(0).maxVal = curve.x.at(edited_point + 1);
                                    }

                                    newBoundaries.at(1).minVal = 0.;
                                    newBoundaries.at(1).maxVal = 1.;
                                    retval = true;
                                    editedPos.at(0) = curve.x.at(edited_point);
                                    editedPos.at(1) = curve.y.at(edited_point);
                                    coordinateAdjuster->switchAdjustedPoint(editedPos, newBoundaries);
                                }
                            } else {
                                // the cursor is inside the graph but away from existing points
                                new_type = CSPlus;
                                curveIsDirty = true;
                                stopNumericalAdjustment();
                            }
                        }

                        retval = true;
                    }
                }
            }

            retval = true;
        }

        break;

    case Gdk::BUTTON_RELEASE:
        snapToElmt = -100;

        if (curve.type != DCT_Parametric && edited_point == -1) {
            if (buttonPressed && event->button.button == 1) {
                std::vector<double>::iterator itx, ity;
                int src, dst;
                buttonPressed = false;
                /*  get the pointer position  */
                getCursorPosition(Gdk::EventType(event->type), event->motion.is_hint != 0, int(event->button.x), int(event->button.y), Gdk::ModifierType(event->button.state));
                findClosestPoint();

                remove_modal_grab ();
                int previous_lit_point = lit_point;
                /* delete inactive points: */
                itx = curve.x.begin();
                ity = curve.y.begin();

                for (src = dst = 0; src < num; ++src)
                    if (curve.x.at(src) >= 0.0) {
                        curve.x.at(dst) = curve.x.at(src);
                        curve.y.at(dst) = curve.y.at(src);
                        ++dst;
                        ++itx;
                        ++ity;
                    }

                if (dst < src) {
                    curve.x.erase (itx, curve.x.end());
                    curve.y.erase (ity, curve.y.end());

                    if (curve.x.empty()) {
                        curve.x.push_back (0);
                        curve.y.push_back (0);
                        curveIsDirty = true;
                        setDirty(true);
                        draw (lit_point);
                    }
                }

                if (distanceX <= minDistanceX) {
                    new_type = CSMove;
                    lit_point = closest_point;
                } else {
                    new_type = CSPlus;
                    lit_point = -1;
                }

                if (lit_point != previous_lit_point) {
                    setDirty(true);
                    draw (lit_point);
                }

                grab_point = -1;
                retval = true;
                notifyListener ();
            }
        }

        break;

    case Gdk::LEAVE_NOTIFY:

        // Pointer can LEAVE even when dragging the point, so we don't modify the cursor in this case
        // The cursor will have to LEAVE another time after the drag...
        if (!buttonPressed)
            if (grab_point == -1) {
                new_type = CSArrow;
                lit_point = -1;
                pipetteR = pipetteG = pipetteB = -1.f;
                setDirty(true);
                draw (lit_point);
            }

        retval = true;
        break;

    case Gdk::MOTION_NOTIFY:
        snapToElmt = -100;

        if (curve.type == DCT_Linear || curve.type == DCT_Spline || curve.type == DCT_NURBS) {

            snapToMinDistY = snapToMinDistX = 10.;
            snapToValY = snapToValX = 0.;
            snapToElmt = -100;

            // get the pointer position
            getCursorPosition(Gdk::EventType(event->type), event->motion.is_hint != 0, int(event->button.x), int(event->button.y), Gdk::ModifierType(event->button.state));

            if (grab_point == -1) {
                if (edited_point == -1) {
                    // there's no point currently being moved
                    int previous_lit_point = lit_point;
                    findClosestPoint();

                    if (cursorX < 0 || cursorX > graphW || cursorY < 0 || cursorY > graphH) {
                        // the cursor has left the graph area
                        new_type = CSArrow;
                        lit_point = -1;
                    } else if (distanceX <= minDistanceX) {
                        // the cursor is close to an existing point
                        new_type = CSMove;
                        lit_point = closest_point;
                    } else {
                        // the cursor is inside the graph but away from existing points
                        new_type = CSPlus;
                        lit_point = -1;
                    }

                    if (lit_point != previous_lit_point) {
                        setDirty(true);
                        draw (lit_point);

                        if (lit_point > -1) {
                            editedPos.at(0) = curve.x.at(lit_point);
                            editedPos.at(1) = curve.y.at(lit_point);
                        }

                        coordinateAdjuster->setPos(editedPos);
                    }

                    if (lit_point == -1 && new_type == CSPlus) {
                        editedPos.at(0) = clampedX;
                        editedPos.at(1) = clampedY;
                        coordinateAdjuster->setPos(editedPos);
                    }
                } else { // if (edited_point > -1)
                    // there's no point currently being moved
                    int previous_lit_point = lit_point;
                    findClosestPoint();

                    if (distanceX <= minDistanceX) {
                        // the cursor is close to an existing point
                        lit_point = closest_point;
                    } else {
                        // the cursor is outside the graph or inside the graph but away from existing points
                        lit_point = -1;
                    }

                    if (lit_point != previous_lit_point) {
                        setDirty(true);
                        draw (lit_point);
                    }
                }
            } else {
                // a point is being moved

                // bounds of the grabbed point
                double leftBound         = (grab_point == 0    ) ? 0. : curve.x.at(grab_point - 1);
                double rightBound        = (grab_point == num - 1) ? 1. : curve.x.at(grab_point + 1);
                double const bottomBound = 0.;
                double const topBound    = 1.;

                double leftDeletionBound   = leftBound   - minDistanceX;
                double rightDeletionBound  = rightBound  + minDistanceX;
                double bottomDeletionBound = bottomBound - minDistanceY;
                double topDeletionBound    = topBound    + minDistanceY;

                // we memorize the previous position of the point, for optimization purpose
                double prevPosX = curve.x.at(grab_point);
                double prevPosY = curve.y.at(grab_point);

                // we memorize the previous position of the point, for optimization purpose
                ugpX += deltaX;
                ugpY += deltaY;

                // the unclamped grabbed point is brought back in the range when snapTo is active
                if (snapTo) {
                    ugpY = CLAMP(ugpY, 0.0, 1.0);
                }

                // handling limitations along X axis
                if (ugpX >= rightDeletionBound && (grab_point > 0 && grab_point < (num - 1))) {
                    curve.x.at(grab_point) = -1.;
                } else if (ugpX <= leftDeletionBound && (grab_point > 0 && grab_point < (num - 1))) {
                    curve.x.at(grab_point) = -1.;
                } else
                    // nextPosX is in bounds
                {
                    curve.x.at(grab_point) = CLAMP(ugpX, leftBound, rightBound);
                }

                // Handling limitations along Y axis
                if (ugpY >= topDeletionBound && grab_point != 0 && grab_point != num - 1) {
                    curve.x.at(grab_point) = -1.;
                } else if (ugpY <= bottomDeletionBound && grab_point != 0 && grab_point != num - 1) {
                    curve.x.at(grab_point) = -1.;
                } else {
                    // snapping point to specific values
                    if (snapTo && curve.x.at(grab_point) != -1.) {
                        if (grab_point > 0 && grab_point < (curve.y.size() - 1)) {
                            double prevX = curve.x.at(grab_point - 1);
                            double prevY = curve.y.at(grab_point - 1);
                            double nextX = curve.x.at(grab_point + 1);
                            double nextY = curve.y.at(grab_point + 1);

                            double ratio = (curve.x.at(grab_point) - prevX) / (nextX - prevX);
                            double y = (nextY - prevY) * ratio + prevY;

                            if (snapCoordinateY(y, ugpY)) {
                                snapToElmt = 1000 + grab_point;
                            }
                        }

                        if (grab_point > 0) {
                            int prevP = grab_point - 1;

                            if (snapCoordinateY(curve.y.at(prevP), ugpY)) {
                                snapToElmt = prevP;
                            }
                        }

                        if (grab_point < (curve.y.size() - 1)) {
                            int nextP = grab_point + 1;

                            if (snapCoordinateY(curve.y.at(nextP), ugpY)) {
                                snapToElmt = nextP;
                            }
                        }

                        if (snapCoordinateY(1.0,                    ugpY)) {
                            snapToElmt = -3;
                        }

                        if (snapCoordinateY(curve.x.at(grab_point), ugpY)) {
                            snapToElmt = -2;
                        }

                        if (snapCoordinateY(0.0,                    ugpY)) {
                            snapToElmt = -1;
                        }

                        curve.y.at(grab_point) = snapToValY;
                    } else {
                        // nextPosY is in the bounds
                        curve.y.at(grab_point) = CLAMP(ugpY, 0.0, 1.0);
                    }
                }

                if (curve.x.at(grab_point) != prevPosX || curve.y.at(grab_point) != prevPosY) {
                    // we recalculate the curve only if we have to
                    curveIsDirty = true;
                    setDirty(true);
                    draw (lit_point);
                    notifyListener ();

                    if (coordinateAdjuster->is_visible()) {
                        editedPos.at(0) = curve.x.at(grab_point);
                        editedPos.at(1) = curve.y.at(grab_point);
                        coordinateAdjuster->setPos(editedPos);
                    }
                }
            }
        }

        retval = true;
        break;

    default:
        break;
    }

    if (new_type != cursor_type) {
        cursor_type = new_type;
        CursorManager::setCursorOfMainWindow(get_window(), cursor_type);
    }

    return retval;
}

CursorShape MyDiagonalCurve::motionNotify(CursorShape type, double minDistanceX, double minDistanceY, int num)
{
    CursorShape new_type = type;

    return new_type;
}

void MyDiagonalCurve::pipetteMouseOver (CurveEditor *ce, EditDataProvider *provider, int modifierKey)
{
    if (!provider) {
        // occurs when leaving the preview area -> cleanup the curve editor
        pipetteR = pipetteG = pipetteB = -1.f;
        lit_point = -1;
        return;
    }

    pipetteR = provider->pipetteVal[0];
    pipetteG = provider->pipetteVal[1];
    pipetteB = provider->pipetteVal[2];
    pipetteVal = 0.f;

    if (listener) {
        pipetteVal = listener->blendPipetteValues(ce, pipetteR, pipetteG, pipetteB);
    } else {
        int n = 0;

        if (pipetteR != -1.f) {
            pipetteVal += pipetteR;
            ++n;
        }

        if (pipetteG != -1.f) {
            pipetteVal += pipetteG;
            ++n;
        }

        if (pipetteB != -1.f) {
            pipetteVal += pipetteB;
            ++n;
        }

        if (n > 1) {
            pipetteVal /= n;
        } else if (!n) {
            pipetteVal = -1.f;
        }
    }

    /* graphW and graphH are the size of the graph */
    calcDimensions();

    if ((graphW < 0) || (graphH < 0)) {
        return;
    }

    double minDistanceX = double(MIN_DISTANCE) / double(graphW - 1);

    if (curve.type == DCT_Linear || curve.type == DCT_Spline || curve.type == DCT_NURBS) {
        // get the pointer position
        getCursorPositionFromCurve(pipetteVal);

        if (edited_point == -1) {
            if (grab_point == -1) {
                // there's no point currently being moved
                int previous_lit_point = lit_point;
                findClosestPoint();

                if (cursorX < 0 || cursorX > graphW || cursorY < 0 || cursorY > graphH) {
                    // the cursor has left the graph area
                    lit_point = -1;
                } else if (distanceX <= minDistanceX) {
                    lit_point = closest_point;
                } else {
                    lit_point = -1;
                }

                if (lit_point != previous_lit_point) {
                    setDirty(true);
                    draw (lit_point);
                }
            }
        } else {
            draw(lit_point);
        }

        if (edited_point == -1) {
            editedPos.at(0) = pipetteVal;
            editedPos.at(1) = point.getVal01(pipetteVal);
            coordinateAdjuster->setPos(editedPos);
        }
    }
}


void MyDiagonalCurve::pipetteButton1Pressed(EditDataProvider *provider, int modifierKey)
{
    if (edited_point > 1) {
        return;
    }

    int num = (int)curve.x.size();

    /* graphW and graphH are the size of the graph */
    calcDimensions();

    double minDistanceX = double(MIN_DISTANCE) / double(graphW - 1);

    if ((graphW < 0) || (graphH < 0)) {
        return;
    }

    snapToElmt = -100;

    if (curve.type != DCT_Parametric) {
        std::vector<double>::iterator itx, ity;
        buttonPressed = true;

        // get the pointer position
        getCursorPositionFromCurve(pipetteVal);
        findClosestPoint();

        if (distanceX > minDistanceX) {
            /* insert a new control point */
            if (num > 0) {
                if (clampedX > curve.x.at(closest_point)) {
                    ++closest_point;
                }
            }

            itx = curve.x.begin();
            ity = curve.y.begin();

            for (int i = 0; i < closest_point; i++) {
                itx++;
                ity++;
            }

            lit_point = closest_point;
            curve.x.insert (itx, 0);
            curve.y.insert (ity, 0);

            // the graph is refreshed only if a new point is created (snapped to a pixel)
            curve.x.at(lit_point) = clampedX;
            curve.y.at(lit_point) = clampedY;

            if (lit_point > -1 && grab_point == -1 && coordinateAdjuster->is_visible()) {
                std::vector<double> position;
                position.resize(2);
                position.at(0) = clampedX;
                position.at(1) = clampedY;
                coordinateAdjuster->setPos(position);
            }

            curveIsDirty = true;
            setDirty(true);
            draw (lit_point);
            notifyListener ();
        }

        grab_point = closest_point;
        lit_point = closest_point;
        ugpX = curve.x.at(closest_point);
        ugpY = curve.y.at(closest_point);
    }
}

void MyDiagonalCurve::pipetteButton1Released(EditDataProvider *provider)
{
    if (edited_point > 1) {
        return;
    }

    /* graphW and graphH are the size of the graph */
    calcDimensions();

    double minDistanceX = double(MIN_DISTANCE) / double(graphW - 1);
    //double minDistanceY = double(MIN_DISTANCE) / double(graphH-1);

    if ((graphW < 0) || (graphH < 0)) {
        return;
    }

    snapToElmt = -100;

    if (curve.type != DCT_Parametric) {
        std::vector<double>::iterator itx, ity;
        buttonPressed = false;
        /*  get the pointer position  */
        getCursorPosition(Gdk::EventType(Gdk::BUTTON_RELEASE), false, graphY, 0, Gdk::ModifierType(0));
        findClosestPoint();

        int previous_lit_point = lit_point;

        if (distanceX <= minDistanceX) {
            lit_point = closest_point;
        } else {
            lit_point = -1;
        }

        if (lit_point != previous_lit_point) {
            setDirty(true);
            draw (lit_point);
        }

        grab_point = -1;
        //notifyListener ();
    }
}

void MyDiagonalCurve::pipetteDrag(EditDataProvider *provider, int modifierKey)
{
    if (edited_point > -1 || curve.type == DCT_Parametric || graphW < 0 || graphH < 0) {
        return;
    }

    snapToMinDistY = snapToMinDistX = 10.;
    snapToValY = snapToValX = 0.;
    snapToElmt = -100;

    /* graphW and graphH are the size of the graph */
    calcDimensions();

    getCursorPosition(Gdk::MOTION_NOTIFY, false, cursorX + graphX, graphY - cursorY + provider->deltaPrevScreen.y, Gdk::ModifierType(modifierKey));

    // we memorize the previous position of the point, for optimization purpose
    double prevPosX = curve.x.at(grab_point);
    double prevPosY = curve.y.at(grab_point);

    // we memorize the previous position of the point, for optimization purpose
    ugpX += deltaX;
    ugpY += deltaY;

    // the unclamped grabbed point is brought back in the range
    ugpY = CLAMP(ugpY, 0.0, 1.0);

    // snapping point to specific values
    if (snapTo && curve.x.at(grab_point) != -1.) {
        if (grab_point > 0 && grab_point < (curve.y.size() - 1)) {
            double prevX = curve.x.at(grab_point - 1);
            double prevY = curve.y.at(grab_point - 1);
            double nextX = curve.x.at(grab_point + 1);
            double nextY = curve.y.at(grab_point + 1);

            double ratio = (curve.x.at(grab_point) - prevX) / (nextX - prevX);
            double y = (nextY - prevY) * ratio + prevY;

            if (snapCoordinateY(y, ugpY)) {
                snapToElmt = 1000 + grab_point;
            }
        }

        if (grab_point > 0) {
            int prevP = grab_point - 1;

            if (snapCoordinateY(curve.y.at(prevP), ugpY)) {
                snapToElmt = prevP;
            }
        }

        if (grab_point < (curve.y.size() - 1)) {
            int nextP = grab_point + 1;

            if (snapCoordinateY(curve.y.at(nextP), ugpY)) {
                snapToElmt = nextP;
            }
        }

        if (snapCoordinateY(1.0,                    ugpY)) {
            snapToElmt = -3;
        }

        if (snapCoordinateY(curve.x.at(grab_point), ugpY)) {
            snapToElmt = -2;
        }

        if (snapCoordinateY(0.0,                    ugpY)) {
            snapToElmt = -1;
        }

        curve.y.at(grab_point) = snapToValY;
    } else {
        // nextPosY is in the bounds
        curve.y.at(grab_point) = ugpY;
    }

    if (curve.x.at(grab_point) != prevPosX || curve.y.at(grab_point) != prevPosY) {
        // we recalculate the curve only if we have to
        curveIsDirty = true;
        setDirty(true);
        draw (lit_point);
        notifyListener ();

        if (lit_point > -1 && coordinateAdjuster->is_visible()) {
            std::vector<double> position;
            position.resize(2);
            position.at(0) = curve.x.at(grab_point);
            position.at(1) = curve.y.at(grab_point);
            coordinateAdjuster->setPos(position);
        }

    }
}

void MyDiagonalCurve::getCursorPositionFromCurve(float x)
{

    // the graph is refreshed only if a new point is created (snaped to a pixel)
    clampedX = x;
    clampedY = point.getVal01(x);

    cursorX = int(clampedX * float(graphW - 3)) + graphX + 1.5;
    cursorY = graphY - int(clampedY * float(graphH - 3));
}

// x = cursor position found in the event
void MyDiagonalCurve::getCursorPositionFromCurve(int x)
{

    // the graph is refreshed only if a new point is created (snaped to a pixel)
    cursorX = x - graphX;
    clampedX = (float(cursorX) - 1.5) / float(graphW - 3);
    clampedY = point.getVal01(clampedX);
    cursorY = graphY - int(float(1. - clampedY) * float(graphH - 3));
}

void MyDiagonalCurve::getCursorPosition(Gdk::EventType evType, bool isHint, int evX, int evY, Gdk::ModifierType modifierKey)
{
    int tx, ty;
    int prevCursorX, prevCursorY;
    double incrementX = 1. / double(graphW);
    double incrementY = 1. / double(graphH);

    // getting the cursor position
    switch (evType) {
    case (Gdk::MOTION_NOTIFY) :
        if (isHint) {
            get_window()->get_pointer (tx, ty, mod_type);
        } else {
            tx = evX;
            ty = evY;
            mod_type = modifierKey;
        }

        break;

    case (Gdk::BUTTON_PRESS) :
    case (Gdk::BUTTON_RELEASE) :
        tx = evX;
        ty = evY;
        mod_type = modifierKey;
        break;

    default :
        // The cursor position is not available
        return;
        break;
    }

    if (grab_point != -1) {
        prevCursorX = cursorX;
        prevCursorY = cursorY;
    }

    cursorX = tx - graphX;
    cursorY = graphY - ty;

    snapTo = ST_None;

    // update deltaX/Y if the user drags a point
    if (grab_point != -1) {
        // set the dragging factor
        int control_key = mod_type & GDK_CONTROL_MASK;
        int shift_key = mod_type & GDK_SHIFT_MASK;

        // the increment get smaller if modifier key are used, and "snap to" may be enabled
        if (control_key) {
            incrementX *= 0.05;
            incrementY *= 0.05;
        }

        if (shift_key)   {
            snapTo = true;
        }

        deltaX = double(cursorX - prevCursorX) * incrementX;
        deltaY = double(cursorY - prevCursorY) * incrementY;
    }
    // otherwise set the position of the new point (modifier keys has no effect here)
    else {
        double tempCursorX = cursorX * incrementX;
        double tempCursorY = cursorY * incrementY;
        clampedX = CLAMP (tempCursorX, 0., 1.);  // X position of the pointer from the origin of the graph
        clampedY = CLAMP (tempCursorY, 0., 1.); // Y position of the pointer from the origin of the graph
    }

}

void MyDiagonalCurve::findClosestPoint()
{
    distanceX = 10.0;
    distanceY = 10.0;
    closest_point = -1;

    if (curve.type != DCT_Parametric) {
        for (int i = 0; i < (int)curve.x.size(); i++) {
            double dX = curve.x.at(i) - clampedX;
            double dY = curve.y.at(i) - clampedY;
            double currDistX = dX < 0. ? -dX : dX; //abs (dX);
            double currDistY = dY < 0. ? -dY : dY; //abs (dY);

            if (currDistX < distanceX) {
                distanceX = currDistX;
                distanceY = currDistY;
                closest_point = i;
            } else if (currDistX == distanceX && currDistY < distanceY) {
                // there is more than 1 point for that X coordinate, we select the closest point to the cursor
                distanceY = currDistY;
                closest_point = i;
            }
        }
    }
}

std::vector<double> MyDiagonalCurve::getPoints ()
{
    std::vector<double> result;

    if (curve.type == DCT_Parametric) {
        result.push_back ((double)(DCT_Parametric));

        for (int i = 0; i < (int)curve.x.size(); i++) {
            result.push_back (curve.x.at(i));
        }
    } else {
        // the first value gives the type of the curve
        if (curve.type == DCT_Linear) {
            result.push_back (double(DCT_Linear));
        } else if (curve.type == DCT_Spline) {
            result.push_back (double(DCT_Spline));
        } else if (curve.type == DCT_NURBS) {
            result.push_back (double(DCT_NURBS));
        }

        // then we push all the points coordinate
        for (int i = 0; i < (int)curve.x.size(); i++) {
            if (curve.x.at(i) >= 0) {
                result.push_back (curve.x.at(i));
                result.push_back (curve.y.at(i));
            }
        }
    }

    return result;
}

void MyDiagonalCurve::setPoints (const std::vector<double>& p)
{
    int ix = 0;
    stopNumericalAdjustment();
    DiagonalCurveType t = (DiagonalCurveType)p[ix++];
    curve.type = t;

    if (t == DCT_Parametric) {
        curve.x.clear ();
        curve.y.clear ();

        for (size_t i = 1; i < p.size(); i++) {
            curve.x.push_back (p[ix++]);
        }
    } else {
        curve.x.clear ();
        curve.y.clear ();

        for (size_t i = 0; i < p.size() / 2; i++) {
            curve.x.push_back (p[ix++]);
            curve.y.push_back (p[ix++]);
        }

        activeParam = -1;
    }

    curveIsDirty = true;
    setDirty(true);
    queue_draw ();
}

void MyDiagonalCurve::setPos(double pos, int chanIdx)
{
    assert (edited_point > -1);

    if (chanIdx == 0) {
        curve.x.at(edited_point) = pos;
    } else if (chanIdx == 1) {
        curve.y.at(edited_point) = pos;
    }

    curveIsDirty = true;
    setDirty(true);
    draw(lit_point);
    notifyListener ();
}

void MyDiagonalCurve::stopNumericalAdjustment()
{
    if (edited_point > -1) {
        edited_point = grab_point = lit_point = -1;
        coordinateAdjuster->stopNumericalAdjustment();
        setDirty(true);
        draw(lit_point);
    }
}

void MyDiagonalCurve::setType (DiagonalCurveType t)
{

    curve.type = t;
    setDirty(true);
}

void MyDiagonalCurve::setActiveParam (int ac)
{

    activeParam = ac;
    setDirty(true);
    queue_draw ();
}

int diagonalmchistupdateUI (void* data)
{

    MyCurveIdleHelper* mcih = static_cast<MyCurveIdleHelper*>(data);

    if (mcih->destroyed) {
        if (mcih->pending == 1) {
            delete mcih;
        } else {
            mcih->pending--;
        }

        return 0;
    }

    mcih->clearPixmap ();
    mcih->myCurve->queue_draw ();

    mcih->pending--;

    return 0;
}

void MyDiagonalCurve::updateBackgroundHistogram (LUTu & hist)
{

    if (hist) {
        //memcpy (bghist, hist, 256*sizeof(unsigned int));
        for (int i = 0; i < 256; i++) {
            bghist[i] = hist[i];
        }

        //hist = bghist;
        bghistvalid = true;
    } else {
        bghistvalid = false;
    }

    mcih->pending++;
    // Can be done outside of the GUI thread, so we're using g_idle_add instead of add_idle
    g_idle_add (diagonalmchistupdateUI, mcih);

}

void MyDiagonalCurve::reset(const std::vector<double> &resetCurve, double identityValue)
{

    stopNumericalAdjustment();

    if (!resetCurve.empty()) {
        setPoints(resetCurve);
        return;
    }

    switch (curve.type) {
    case DCT_Spline :
    case DCT_NURBS :
        curve.x.resize(2);
        curve.y.resize(2);
        curve.x.at(0) = 0.;
        curve.y.at(0) = 0.;
        curve.x.at(1) = 1.;
        curve.y.at(1) = 1.;
        grab_point = -1;
        lit_point = -1;
        curveIsDirty = true;
        break;

    case DCT_Parametric :
        curve.x.resize(7);
        curve.y.clear();
        // the SHCSelector values doesn't really matter for the identity curve display
        curve.x.at(0) = 0.25;
        curve.x.at(1) = 0.50;
        curve.x.at(2) = 0.75;
        curve.x.at(3) = 0.00;
        curve.x.at(4) = 0.00;
        curve.x.at(5) = 0.00;
        curve.x.at(6) = 0.00;
        grab_point = -1;  // not sure that it's necessary
        lit_point = -1;   // not sure that it's necessary
        curveIsDirty = true;
        break;

    default:
        break;
    }

    setDirty(true);
    draw(-1);
}

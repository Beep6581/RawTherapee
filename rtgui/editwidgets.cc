/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Jean-Christophe FRISCH <natureh.510@gmail.com>
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

#include "editwidgets.h"

#include "editbuffer.h"
#include "editcallbacks.h"
#include "rtsurface.h"

const std::vector<double> Geometry::dash = {3., 1.5};

#define INNERGEOM_OPACITY 1.
#define OUTERGEOM_OPACITY 0.7

RGBColor Geometry::getInnerLineColor ()
{
    RGBColor color;

    if (flags & F_AUTO_COLOR) {
        if      (state == NORMAL)   {
            color.setColor (1., 1., 1.);           // White
        } else if (state == ACTIVE) {
            color.setColor (1., 1., 0.);           // Yellow
        } else if (state == PRELIGHT) {
            color.setColor (1., 100. / 255., 0.);  // Orange
        } else if (state == DRAGGED) {
            color.setColor (1., 0., 0.);           // Red
        }
    } else {
        color = innerLineColor;
    }

    return color;
}

RGBColor Geometry::getOuterLineColor ()
{
    RGBColor color;

    if (flags & F_AUTO_COLOR) {
        /*
        if      (state == NORMAL)   { color.setColor (0., 0., 0.); }  // Black
        else if (state == PRELIGHT) { color.setColor (0., 0., 0.); }  // Black
        else if (state == DRAGGED)  { color.setColor (1., 0., 0.); }  // Black
        */
        color.setColor (0., 0., 0.);  // Black
    } else {
        color = outerLineColor;
    }

    return color;
}

#ifdef GUIVERSION

void Circle::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    double lineWidth = getOuterLineWidth();
    if ((flags & F_VISIBLE) && state != INSENSITIVE && lineWidth > 0. && innerLineWidth > 0.) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgba (color.getR(), color.getG(), color.getB(), OUTERGEOM_OPACITY * rtengine::min(innerLineWidth / 2.f, 1.f));
        cr->set_line_width (lineWidth);
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);
        cr->stroke();
    }
}

void Circle::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgba (color.getR(), color.getG(), color.getB(), INNERGEOM_OPACITY);
        }

        cr->set_line_width(innerLineWidth);
        cr->set_line_cap(flags & F_DASHED ? Cairo::LINE_CAP_BUTT : Cairo::LINE_CAP_ROUND);

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (flags & F_DASHED) {
            cr->set_dash(dash, 0.);
        }

        if (filled) {
            cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);
            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.f) {
            cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);
            cr->stroke();
        }

        if (flags & F_DASHED) {
            cr->unset_dash();
        }
}
}

void Circle::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);
        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        // setting the color to the objet's ID
        if (objectBuffer->getObjectMode() == OM_255) {
            cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
        } else {
            cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
        }
        cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0, 2.*rtengine::RT_PI);

        if (filled) {
            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else {
            cr->stroke();
        }
    }
}

void Line::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    double lineWidth = getOuterLineWidth();
    if ((flags & F_VISIBLE) && state != INSENSITIVE && lineWidth > 0. && innerLineWidth > 0.) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgba (color.getR(), color.getG(), color.getB(), OUTERGEOM_OPACITY * rtengine::min(innerLineWidth / 2.f, 1.f));
        cr->set_line_width (lineWidth);
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);

        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);
        cr->stroke();
    }
}

void Line::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && innerLineWidth > 0.f) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgba (color.getR(), color.getG(), color.getB(), INNERGEOM_OPACITY);
        }

        cr->set_line_width(innerLineWidth);
        cr->set_line_cap(flags & F_DASHED ? Cairo::LINE_CAP_BUTT : Cairo::LINE_CAP_ROUND);

        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (flags & F_DASHED) {
            cr->set_dash(dash, 0.);
        }

        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);
        cr->stroke();

        if (flags & F_DASHED) {
            cr->unset_dash();
        }
    }
}

void Line::drawToMOChannel(Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);
        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToCropCanvas (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        // setting the color to the objet's ID
        if (objectBuffer->getObjectMode() == OM_255) {
            cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
        } else {
            cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
        }
        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);
        cr->stroke();
    }
}

void Polyline::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    double lineWidth = getOuterLineWidth();
    if ((flags & F_VISIBLE) && state != INSENSITIVE && points.size() > 1 && lineWidth > 0.) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgba (color.getR(), color.getG(), color.getB(), OUTERGEOM_OPACITY * rtengine::min(innerLineWidth / 2.f, 1.f));
        cr->set_line_width (lineWidth);
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);

        rtengine::Coord currPos;

        for (unsigned int i = 0; i < points.size(); ++i) {
            currPos  = points.at(i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToScreen (points.at(i).x, points.at(i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += objectBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            }

            if (!i) {
                cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
            } else {
                cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
            }
        }

        if (filled) {
            cr->fill_preserve();
            cr->stroke();
        } else {
            cr->stroke();
        }
    }
}

void Polyline::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && points.size() > 1) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgba (color.getR(), color.getG(), color.getB(), INNERGEOM_OPACITY);
        }

        cr->set_line_width(innerLineWidth);
        cr->set_line_cap(flags & F_DASHED ? Cairo::LINE_CAP_BUTT : Cairo::LINE_CAP_ROUND);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);

        if (flags & F_DASHED) {
            cr->set_dash(dash, 0.);
        }

        if (filled && state != INSENSITIVE) {
            rtengine::Coord currPos;

            for (unsigned int i = 0; i < points.size(); ++i) {
                currPos  = points.at(i);

                if      (datum == IMAGE) {
                    coordSystem.imageCoordToScreen (points.at(i).x, points.at(i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += objectBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
                }
            }

            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.f) {
            rtengine::Coord currPos;

            for (unsigned int i = 0; i < points.size(); ++i) {
                currPos  = points.at(i);

                if (datum == IMAGE) {
                    coordSystem.imageCoordToScreen (points.at(i).x, points.at(i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += objectBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
                }
            }
            cr->stroke();
        }

        if (flags & F_DASHED) {
            cr->unset_dash();
        }
    }
}

void Polyline::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_HOVERABLE) && points.size() > 1) {
        rtengine::Coord currPos;

        // setting the color to the objet's ID
        if (objectBuffer->getObjectMode() == OM_255) {
            cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
        } else {
            cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
        }

        cr->set_line_width( getMouseOverLineWidth() );
        cr->set_line_cap(Cairo::LINE_CAP_ROUND);
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);

        for (unsigned int i = 0; i < points.size(); ++i) {
            currPos  = points.at(i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToCropCanvas (points.at(i).x, points.at(i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += objectBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            }

            if (!i) {
                cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
            } else {
                cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
            }
        }

        if (filled) {
            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else {
            cr->stroke();
        }
    }
}

void EditRectangle::setXYWH(int left, int top, int width, int height)
{
    topLeft.set(left, top);
    bottomRight.set(left + width, top + height);
}

void EditRectangle::setXYXY(int left, int top, int right, int bottom)
{
    topLeft.set(left, top);
    bottomRight.set(right, bottom);
}

void EditRectangle::setXYWH(rtengine::Coord topLeft, rtengine::Coord widthHeight)
{
    this->topLeft = topLeft;
    this->bottomRight = topLeft + widthHeight;
}

void EditRectangle::setXYXY(rtengine::Coord topLeft, rtengine::Coord bottomRight)
{
    this->topLeft = topLeft;
    this->bottomRight = bottomRight;
}

void EditRectangle::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    double lineWidth = getOuterLineWidth();
    if ((flags & F_VISIBLE) && state != INSENSITIVE && lineWidth > 0. && innerLineWidth > 0.) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgba (color.getR(), color.getG(), color.getB(), OUTERGEOM_OPACITY * rtengine::min(innerLineWidth / 2.f, 1.f));
        cr->set_line_width (lineWidth);
        cr->set_line_join(Cairo::LINE_JOIN_BEVEL);

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen (topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen (bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

        if (filled) {
            cr->fill_preserve();
            cr->stroke();
        } else {
            cr->stroke();
        }
    }
}

void EditRectangle::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgba (color.getR(), color.getG(), color.getB(), INNERGEOM_OPACITY);
        }

        cr->set_line_width(innerLineWidth);
        cr->set_line_join(Cairo::LINE_JOIN_BEVEL);

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen (topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen (bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (flags & F_DASHED) {
            cr->set_dash(dash, 0.);
        }

        if (filled) {
            cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.f) {
            cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);
            cr->stroke();
        }

        if (flags & F_DASHED) {
            cr->unset_dash();
        }
    }
}

void EditRectangle::drawToMOChannel(Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );
        cr->set_line_join(Cairo::LINE_JOIN_ROUND);

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        // setting the color to the objet's ID
        if (objectBuffer->getObjectMode() == OM_255) {
            cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
        } else {
            cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
        }
        cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

        if (filled) {
            if (innerLineWidth > 0.f) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else {
            cr->stroke();
        }
    }
}

void Ellipse::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgba (color.getR(), color.getG(), color.getB(), opacity / 100.);
        cr->set_line_width ( getOuterLineWidth() );

        rtengine::Coord center_ = center;
        double radYT_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radYT)) : double (radYT);
        double radY_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radY)) : double (radY);
        double radXL_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radXL)) : double (radXL);
        double radX_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radX)) : double (radX);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (radYT_ > 0 && radY_ > 0 && radXL_ > 0 && radX_ > 0) {
            // To have an ellipse with radius of (radX, radX), a circle of radius 1. shall be twisted with a scale
            // of radX for x-axis, radY for y-axis
            // Center of coordinates (x, y) in previous coordinates system becomes (X, Y) = (radX * x, radY * y) in new one
            // To go back to previous location, center shall be translated to tx = -X * (1 - 1 / radX) in x-axis (x = tx + X)
            // and ty = -Y * (1 - 1 / radY) in y-axis (y = ty + Y)
            cr->save();

            // Drawing bottom-right part
            cr->scale (radX_, radY_);
            cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radY_));
            cr->arc (center_.x, center_.y, 1.0, 0.0, rtengine::RT_PI_2);

            cr->restore ();
            cr->save();

            // Drawing bottom-left part
            cr->scale (radXL_, radY_);
            cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radY_));
            cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI_2, rtengine::RT_PI);
            cr->scale (radXL_, radY_);

            cr->restore ();
            cr->save();

            // Drawing top-left part
            cr->scale (radXL_, radYT_);
            cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radYT_));
            cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI, 3. * rtengine::RT_PI_2);

            cr->restore ();
            cr->save();

            // Drawing top-right part
            cr->scale (radX_, radYT_);
            cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radYT_));
            cr->arc (center_.x, center_.y, 1.0, 3. * rtengine::RT_PI_2, 2. * rtengine::RT_PI);

            cr->restore ();
            cr->stroke ();
        }
    }
}

void Ellipse::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgba (color.getR(), color.getG(), color.getB(), opacity / 100.);
        }

        cr->set_line_width ( innerLineWidth );

        rtengine::Coord center_ = center;
        double radYT_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radYT)) : double (radYT);
        double radY_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radY)) : double (radY);
        double radXL_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radXL)) : double (radXL);
        double radX_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radX)) : double (radX);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (filled && state != INSENSITIVE) {
            if (radYT_ > 0 && radY_ > 0 && radXL_ > 0 && radX_ > 0) {
                // To have an ellipse with radius of (radX, radX), a circle of radius 1. shall be twisted with a scale
                // of radX for x-axis, radY for y-axis
                // Center of coordinates (x, y) in previous coordinates system becomes (X, Y) = (radX * x, radY * y) in new one
                // To go back to previous location, center shall be translated to tx = -X * (1 - 1 / radX) in x-axis (x = tx + X)
                // and ty = -Y * (1 - 1 / radY) in y-axis (y = ty + Y)
                cr->save();

                // Drawing bottom-right part
                cr->scale (radX_, radY_);
                cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radY_));
                cr->arc (center_.x, center_.y, 1.0, 0.0, rtengine::RT_PI_2);

                cr->restore ();
                cr->save();

                // Drawing bottom-left part
                cr->scale (radXL_, radY_);
                cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radY_));
                cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI_2, rtengine::RT_PI);
                cr->scale (radXL_, radY_);

                cr->restore ();
                cr->save();

                // Drawing top-left part
                cr->scale (radXL_, radYT_);
                cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radYT_));
                cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI, 3. * rtengine::RT_PI_2);

                cr->restore ();
                cr->save();

                // Drawing top-right part
                cr->scale (radX_, radYT_);
                cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radYT_));
                cr->arc (center_.x, center_.y, 1.0, 3. * rtengine::RT_PI_2, 2. * rtengine::RT_PI);

                cr->restore ();
                cr->stroke ();
            }

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            if (radYT_ > 0 && radY_ > 0 && radXL_ > 0 && radX_ > 0) {
                // To have an ellipse with radius of (radX, radX), a circle of radius 1. shall be twisted with a scale
                // of radX for x-axis, radY for y-axis
                // Center of coordinates (x, y) in previous coordinates system becomes (X, Y) = (radX * x, radY * y) in new one
                // To go back to previous location, center shall be translated to tx = -X * (1 - 1 / radX) in x-axis (x = tx + X)
                // and ty = -Y * (1 - 1 / radY) in y-axis (y = ty + Y)
                cr->save();

                // Drawing bottom-right part
                cr->scale (radX_, radY_);
                cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radY_));
                cr->arc (center_.x, center_.y, 1.0, 0.0, rtengine::RT_PI_2);

                cr->restore ();
                cr->save();

                // Drawing bottom-left part
                cr->scale (radXL_, radY_);
                cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radY_));
                cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI_2, rtengine::RT_PI);
                cr->scale (radXL_, radY_);

                cr->restore ();
                cr->save();

                // Drawing top-left part
                cr->scale (radXL_, radYT_);
                cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radYT_));
                cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI, 3. * rtengine::RT_PI_2);

                cr->restore ();
                cr->save();

                // Drawing top-right part
                cr->scale (radX_, radYT_);
                cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radYT_));
                cr->arc (center_.x, center_.y, 1.0, 3. * rtengine::RT_PI_2, 2. * rtengine::RT_PI);

                cr->restore ();
                cr->stroke ();
            }

            if (state == INSENSITIVE) {
                std::valarray<double> ds (1);
                ds[0] = 4;
                cr->set_source_rgba (1.0, 1.0, 1.0, 0.618);
                cr->stroke_preserve();
                cr->set_source_rgba (0.0, 0.0, 0.0, 0.618);
                cr->set_dash (ds, 0);
                cr->stroke();
                ds.resize (0);
                cr->set_dash (ds, 0);
            } else {
                cr->stroke();
            }
        }
    }
}

void Ellipse::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );

        rtengine::Coord center_ = center;
        double radYT_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radYT)) : double (radYT);
        double radY_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radY)) : double (radY);
        double radXL_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radXL)) : double (radXL);
        double radX_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radX)) : double (radX);

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (radYT_ > 0 && radY_ > 0 && radXL_ > 0 && radX_ > 0) {
            // To have an ellipse with radius of (radX, radX), a circle of radius 1. shall be twisted with a scale
            // of radX for x-axis, radY for y-axis
            // Center of coordinates (x, y) in previous coordinates system becomes (X, Y) = (radX * x, radY * y) in new one
            // To go back to previous location, center shall be translated to tx = -X * (1 - 1 / radX) in x-axis (x = tx + X)
            // and ty = -Y * (1 - 1 / radY) in y-axis (y = ty + Y)
            cr->save();

            // Drawing bottom-right part
            cr->scale (radX_, radY_);
            cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radY_));
            cr->arc (center_.x, center_.y, 1.0, 0.0, rtengine::RT_PI_2);

            cr->restore ();
            cr->save();

            // Drawing bottom-left part
            cr->scale (radXL_, radY_);
            cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radY_));
            cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI_2, rtengine::RT_PI);
            cr->scale (radXL_, radY_);

            cr->restore ();
            cr->save();

            // Drawing top-left part
            cr->scale (radXL_, radYT_);
            cr->translate(- center_.x * (1 - 1 / radXL_), - center_.y * (1 - 1 / radYT_));
            cr->arc (center_.x, center_.y, 1.0, rtengine::RT_PI, 3. * rtengine::RT_PI_2);

            cr->restore ();
            cr->save();

            // Drawing top-right part
            cr->scale (radX_, radYT_);
            cr->translate(- center_.x * (1 - 1 / radX_), - center_.y * (1 - 1 / radYT_));
            cr->arc (center_.x, center_.y, 1.0, 3. * rtengine::RT_PI_2, 2. * rtengine::RT_PI);

            cr->restore ();
            cr->stroke ();
        }

        if (filled) {
            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else {
            cr->stroke();
        }
    }
}

void OPIcon::drivenPointToRectangle(const rtengine::Coord &pos,
                                    rtengine::Coord &topLeft, rtengine::Coord &bottomRight, int W, int H)
{
    switch (drivenPoint) {
    case (DP_CENTERCENTER):
        topLeft.x = pos.x - W / 2;
        topLeft.y = pos.y - H / 2;
        break;

    case (DP_TOPLEFT):
        topLeft.x = pos.x;
        topLeft.y = pos.y;
        break;

    case (DP_TOPCENTER):
        topLeft.x = pos.x - W / 2;
        topLeft.y = pos.y;
        break;

    case (DP_TOPRIGHT):
        topLeft.x = pos.x - W;
        topLeft.y = pos.y;
        break;

    case (DP_CENTERRIGHT):
        topLeft.x = pos.x - W;
        topLeft.y = pos.y - H / 2;
        break;

    case (DP_BOTTOMRIGHT):
        topLeft.x = pos.x - W;
        topLeft.y = pos.y - H;
        break;

    case (DP_BOTTOMCENTER):
        topLeft.x = pos.x - W / 2;
        topLeft.y = pos.y - H;
        break;

    case (DP_BOTTOMLEFT):
        topLeft.x = pos.x;
        topLeft.y = pos.y - H;
        break;

    case (DP_CENTERLEFT):
        topLeft.x = pos.x;
        topLeft.y = pos.y - H / 2;
        break;
    }

    bottomRight.x = topLeft.x + W - 1;
    bottomRight.y = topLeft.y + H - 1;
}

OPIcon::OPIcon(const Cairo::RefPtr<RTSurface> &normal,
               const Cairo::RefPtr<RTSurface> &active,
               const Cairo::RefPtr<RTSurface> &prelight,
               const Cairo::RefPtr<RTSurface> &dragged,
               const Cairo::RefPtr<RTSurface> &insensitive,
               DrivenPoint drivenPoint) :
    drivenPoint(drivenPoint)
{
    if (normal) {
        normalImg = normal;
    }

    if (prelight) {
        prelightImg = prelight;
    }

    if (active) {
        activeImg = active;
    }

    if (dragged) {
        draggedImg = dragged;
    }

    if (insensitive) {
        insensitiveImg = insensitive;
    }
}

OPIcon::OPIcon(Glib::ustring normalImage, Glib::ustring activeImage, Glib::ustring prelightImage,
               Glib::ustring  draggedImage, Glib::ustring insensitiveImage, DrivenPoint drivenPoint) : drivenPoint(drivenPoint)
{
    if (!normalImage.empty()) {
        normalImg = Cairo::RefPtr<RTSurface>(new RTSurface(normalImage));
    }

    if (!prelightImage.empty()) {
        prelightImg = Cairo::RefPtr<RTSurface>(new RTSurface(prelightImage));
    }

    if (!activeImage.empty()) {
        activeImg = Cairo::RefPtr<RTSurface>(new RTSurface(activeImage));
    }

    if (!draggedImage.empty()) {
        draggedImg = Cairo::RefPtr<RTSurface>(new RTSurface(draggedImage));
    }

    if (!insensitiveImage.empty()) {
        insensitiveImg = Cairo::RefPtr<RTSurface>(new RTSurface(insensitiveImage));
    }
}

const Cairo::RefPtr<RTSurface> OPIcon::getNormalImg()
{
    return normalImg;
}
const Cairo::RefPtr<RTSurface> OPIcon::getPrelightImg()
{
    return prelightImg;
}
const Cairo::RefPtr<RTSurface> OPIcon::getActiveImg()
{
    return activeImg;
}
const Cairo::RefPtr<RTSurface> OPIcon::getDraggedImg()
{
    return draggedImg;
}
const Cairo::RefPtr<RTSurface> OPIcon::getInsensitiveImg()
{
    return insensitiveImg;
}

void OPIcon::drawImage(Cairo::RefPtr<RTSurface> &img,
                       Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer,
                       EditCoordSystem &coordSystem)
{
    int imgW = img->getWidth();
    int imgH = img->getHeight();

    rtengine::Coord pos;

    if (datum == IMAGE) {
        coordSystem.imageCoordToScreen(position.x, position.y, pos.x, pos.y);
    } else if (datum == CLICKED_POINT) {
        pos = position + objectBuffer->getDataProvider()->posScreen;
    } else if (datum == CURSOR)
        pos = position + objectBuffer->getDataProvider()->posScreen
              + objectBuffer->getDataProvider()->deltaScreen;

    rtengine::Coord tl, br; // Coordinate of the rectangle in the CropBuffer coordinate system
    drivenPointToRectangle(pos, tl, br, imgW, imgH);

    cr->set_source(img->get(), tl.x, tl.y);
    cr->set_line_width(0.);
    cr->rectangle(tl.x, tl.y, imgW, imgH);
    cr->fill();
}

void OPIcon::drawMOImage(Cairo::RefPtr<RTSurface> &img, Cairo::RefPtr<Cairo::Context> &cr,
                         unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    // test of F_HOVERABLE has already been done

    int imgW = img->getWidth();
    int imgH = img->getHeight();

    rtengine::Coord pos;

    if (datum == IMAGE)
        coordSystem.imageCoordToCropCanvas (position.x, position.y, pos.x, pos.y);
    else if (datum == CLICKED_POINT) {
        pos = position + objectBuffer->getDataProvider()->posScreen;
    } else if (datum == CURSOR)
        pos = position + objectBuffer->getDataProvider()->posScreen
              + objectBuffer->getDataProvider()->deltaScreen;

    rtengine::Coord tl, br; // Coordinate of the rectangle in the CropBuffer coordinate system
    drivenPointToRectangle(pos, tl, br, imgW, imgH);

    // drawing the lower byte's value
    if (objectBuffer->getObjectMode() == OM_255) {
        cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
    } else {
        cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
    }
    cr->set_line_width(0.);
    cr->rectangle(tl.x, tl.y, imgW, imgH);
    cr->fill();
}

void OPIcon::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr,
                               ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem) {}

void OPIcon::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr,
                               ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        // Here we will handle fall-back solutions

        State tmpState = state;  // can be updated through the successive test

        if (tmpState == INSENSITIVE) {
            if (!insensitiveImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawImage(insensitiveImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == DRAGGED) {
            if (!draggedImg) {
                tmpState = ACTIVE;
            } else {
                OPIcon::drawImage(draggedImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == ACTIVE) {
            if (!activeImg) {
                tmpState = PRELIGHT;
            } else {
                OPIcon::drawImage(activeImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == PRELIGHT) {
            if (!prelightImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawImage(prelightImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == NORMAL && normalImg) {
            OPIcon::drawImage(normalImg, cr, objectBuffer, coordSystem);
        }
    }
}

void OPIcon::drawToMOChannel(Cairo::RefPtr<Cairo::Context> &cr, unsigned short id,
                             ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        // Here we will handle fallback solutions
        State tmpState = state;

        if (tmpState == INSENSITIVE) {
            if (!insensitiveImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawMOImage(insensitiveImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == DRAGGED) {
            if (!draggedImg) {
                tmpState = ACTIVE;
            } else {
                OPIcon::drawMOImage(draggedImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == ACTIVE) {
            if (!activeImg) {
                tmpState = PRELIGHT;
            } else {
                OPIcon::drawMOImage(activeImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == PRELIGHT) {
            if (!prelightImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawMOImage(prelightImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == NORMAL && normalImg) {
            OPIcon::drawMOImage(normalImg, cr, id, objectBuffer, coordSystem);
        }
    }
}

#endif

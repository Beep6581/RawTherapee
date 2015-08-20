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

#include "edit.h"
#include "../rtengine/editbuffer.h"
#include "rtimage.h"

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

void Circle::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width( getOuterLineWidth() );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToScreen(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen(center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*M_PI);
        cr->stroke();
    }
}

void Circle::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgb(color.getR(), color.getG(), color.getB());
        }

        cr->set_line_width( innerLineWidth );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToScreen(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen(center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        if (filled && state != INSENSITIVE) {
            cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*M_PI);

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*M_PI);

            if (state == INSENSITIVE) {
                std::valarray<double> ds(1);
                ds[0] = 4;
                cr->set_source_rgba(1.0, 1.0, 1.0, 0.618);
                cr->stroke_preserve();
                cr->set_source_rgba(0.0, 0.0, 0.0, 0.618);
                cr->set_dash(ds, 0);
                cr->stroke();
                ds.resize(0);
                cr->set_dash(ds, 0);
            } else {
                cr->stroke();
            }
        }
    }
}

void Circle::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );
        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToScreen(double(radius)) : double(radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropBuffer(center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        // drawing the lower byte's value
        unsigned short a = (id + 1) & 0xFF;
        cr->set_source_rgba (0., 0., 0., double(a) / 255.);
        cr->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0, 2.*M_PI);

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

        // drawing the higher byte's value
        if (editBuffer->getObjectMode() == OM_65535) {
            a = (id + 1) >> 8;
            cr2->set_source_rgba (0., 0., 0., double(a) / 255.);
            cr2->arc(center_.x + 0.5, center_.y + 0.5, radius_, 0, 2.*M_PI);

            if (filled) {
                if (innerLineWidth > 0.) {
                    cr2->fill_preserve();
                    cr2->stroke();
                } else {
                    cr2->fill();
                }
            } else {
                cr2->stroke();
            }
        }
    }
}

void Line::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width( getOuterLineWidth() );

        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen(begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen(end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += editBuffer->getDataProvider()->posScreen;
            end_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
            end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);
        cr->stroke();
    }
}

void Line::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && innerLineWidth > 0.) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        }

        cr->set_line_width(innerLineWidth);

        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen(begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen(end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += editBuffer->getDataProvider()->posScreen;
            end_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
            end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);

        if (state == INSENSITIVE) {
            std::valarray<double> ds(1);
            ds[0] = 4;
            cr->set_source_rgba(1.0, 1.0, 1.0, 0.618);
            cr->stroke_preserve();
            cr->set_source_rgba(0.0, 0.0, 0.0, 0.618);
            cr->set_dash(ds, 0);
            cr->stroke();
            ds.resize(0);
            cr->set_dash(ds, 0);
        } else {
            cr->stroke();
        }
    }
}

void Line::drawToMOChannel(Cairo::RefPtr<Cairo::Context> &cr,
                           Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id,
                           rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );
        rtengine::Coord begin_ = begin;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropBuffer(begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToCropBuffer(end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += editBuffer->getDataProvider()->posScreen;
            end_ += editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
            end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        // drawing the lower byte's value
        unsigned short a = (id + 1) & 0xFF;
        cr->set_source_rgba (0., 0., 0., double(a) / 255.);
        cr->move_to(begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to(end_.x + 0.5, end_.y + 0.5);
        cr->stroke();

        // drawing the higher byte's value
        if (editBuffer->getObjectMode() == OM_65535) {
            a = (id + 1) >> 8;
            cr2->set_source_rgba (0., 0., 0., double(a) / 255.);
            cr2->move_to(begin_.x + 0.5, begin_.y + 0.5);
            cr2->line_to(end_.x + 0.5, end_.y + 0.5);
            cr2->stroke();
        }
    }
}

void Polyline::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE && points.size() > 1) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width( getOuterLineWidth() );

        rtengine::Coord currPos;

        for (unsigned int i = 0; i < points.size(); ++i) {
            currPos  = points.at(i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += editBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
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

void Polyline::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && points.size() > 1) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        }

        cr->set_line_width( innerLineWidth );

        if (filled && state != INSENSITIVE) {
            rtengine::Coord currPos;

            for (unsigned int i = 0; i < points.size(); ++i) {
                currPos  = points.at(i);

                if      (datum == IMAGE) {
                    coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += editBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
                }
            }

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            rtengine::Coord currPos;

            for (unsigned int i = 0; i < points.size(); ++i) {
                currPos  = points.at(i);

                if (datum == IMAGE) {
                    coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += editBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
                }
            }

            if (state == INSENSITIVE) {
                std::valarray<double> ds(1);
                ds[0] = 4;
                cr->set_source_rgba(1.0, 1.0, 1.0, 0.618);
                cr->stroke_preserve();
                cr->set_source_rgba(0.0, 0.0, 0.0, 0.618);
                cr->set_dash(ds, 0);
                cr->stroke();
                ds.resize(0);
                cr->set_dash(ds, 0);
            } else {
                cr->stroke();
            }
        }
    }
}

void Polyline::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_HOVERABLE) && points.size() > 1) {
        rtengine::Coord currPos;

        // drawing the lower byte's value
        unsigned short a = (id + 1) & 0xFF;
        cr->set_source_rgba (0., 0., 0., double(a) / 255.);

        for (unsigned int i = 0; i < points.size(); ++i) {
            cr->set_line_width( getMouseOverLineWidth() );
            currPos  = points.at(i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToCropBuffer(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += editBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
            }

            if (!i) {
                cr->move_to(currPos.x + 0.5, currPos.y + 0.5);
            } else {
                cr->line_to(currPos.x + 0.5, currPos.y + 0.5);
            }
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

        // drawing the higher byte's value
        if (editBuffer->getObjectMode() == OM_65535) {
            a = (id + 1) >> 8;
            cr2->set_source_rgba (0., 0., 0., double(a) / 255.);

            for (unsigned int i = 0; i < points.size(); ++i) {
                cr2->set_line_width( getMouseOverLineWidth() );
                currPos  = points.at(i);

                if      (datum == IMAGE) {
                    coordSystem.imageCoordToCropBuffer(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += editBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr2->move_to(currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr2->line_to(currPos.x + 0.5, currPos.y + 0.5);
                }
            }

            if (filled) {
                if (innerLineWidth > 0.) {
                    cr2->fill_preserve();
                    cr2->stroke();
                } else {
                    cr2->fill();
                }
            } else {
                cr2->stroke();
            }
        }
    }
}

void Rectangle::setXYWH(int left, int top, int width, int height)
{
    topLeft.set(left, top);
    bottomRight.set(left + width, top + height);
}

void Rectangle::setXYXY(int left, int top, int right, int bottom)
{
    topLeft.set(left, top);
    bottomRight.set(right, bottom);
}

void Rectangle::setXYWH(rtengine::Coord topLeft, rtengine::Coord widthHeight)
{
    this->topLeft = topLeft;
    this->bottomRight = topLeft + widthHeight;
}

void Rectangle::setXYXY(rtengine::Coord topLeft, rtengine::Coord bottomRight)
{
    this->topLeft = topLeft;
    this->bottomRight = bottomRight;
}

void Rectangle::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width( getOuterLineWidth() );

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen(topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen(bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
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

void Rectangle::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        }

        cr->set_line_width( innerLineWidth );

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen(topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToScreen(bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        if (filled && state != INSENSITIVE) {
            cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

            if (state == INSENSITIVE) {
                std::valarray<double> ds(1);
                ds[0] = 4;
                cr->set_source_rgba(1.0, 1.0, 1.0, 0.618);
                cr->stroke_preserve();
                cr->set_source_rgba(0.0, 0.0, 0.0, 0.618);
                cr->set_dash(ds, 0);
                cr->stroke();
                ds.resize(0);
                cr->set_dash(ds, 0);
            } else {
                cr->stroke();
            }
        }
    }
}

void Rectangle::drawToMOChannel(Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width( getMouseOverLineWidth() );

        rtengine::Coord tl, br;

        if      (datum == IMAGE) {
            coordSystem.imageCoordToCropBuffer(topLeft.x, topLeft.y, tl.x, tl.y);
        } else if (datum == CLICKED_POINT) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            tl = topLeft + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        if      (datum == IMAGE) {
            coordSystem.imageCoordToCropBuffer(bottomRight.x, bottomRight.y, br.x, br.y);
        } else if (datum == CLICKED_POINT) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            br = bottomRight + editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaScreen;
        }

        // drawing the lower byte's value
        unsigned short a = (id + 1) & 0xFF;
        cr->set_source_rgba (0., 0., 0., double(a) / 255.);
        cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

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

        // drawing the higher byte's value
        if (editBuffer->getObjectMode() == OM_65535) {
            a = (id + 1) >> 8;
            cr2->set_source_rgba (0., 0., 0., double(a) / 255.);
            cr->rectangle(tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

            if (filled) {
                if (innerLineWidth > 0.) {
                    cr2->fill_preserve();
                    cr2->stroke();
                } else {
                    cr2->fill();
                }
            } else {
                cr2->stroke();
            }
        }
    }
}

EditSubscriber::EditSubscriber (EditType editType) : ID(EUID_None), editingType(editType), bufferType(BT_SINGLEPLANE_FLOAT), provider(NULL), dragging(false) {}

void EditSubscriber::setEditProvider(EditDataProvider *provider)
{
    this->provider = provider;
}

void EditSubscriber::setEditID(EditUniqueID ID, BufferType buffType)
{
    this->ID = ID;
    bufferType = buffType;
}

bool EditSubscriber::isCurrentSubscriber()
{
    //if (provider && provider->getCurrSubscriber())
    //  return provider->getCurrSubscriber()->getEditID() == ID;

    if (provider) {
        return provider->getCurrSubscriber() == this;
    }

    return false;
}

void EditSubscriber::subscribe()
{
    if (provider) {
        provider->subscribe(this);
    }
}

void EditSubscriber::unsubscribe()
{
    if (provider) {
        provider->unsubscribe();
    }
}

void EditSubscriber::switchOffEditMode()
{
    unsubscribe();
}

EditUniqueID EditSubscriber::getEditID()
{
    return ID;
}

EditType EditSubscriber::getEditingType()
{
    return editingType;
}

BufferType EditSubscriber::getEditBufferType()
{
    return bufferType;
}

bool EditSubscriber::isDragging()
{
    return dragging;
}

//--------------------------------------------------------------------------------------------------


EditDataProvider::EditDataProvider() : currSubscriber(NULL), object(0), posScreen(-1, -1), posImage(-1, -1),
    deltaScreen(0, 0), deltaImage(0, 0), deltaPrevScreen(0, 0), deltaPrevImage(0, 0)
{
    pipetteVal[0] = pipetteVal[1] = pipetteVal[2] = 0.f;
}

void EditDataProvider::subscribe(EditSubscriber *subscriber)
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode();
    }

    currSubscriber = subscriber;
}

void EditDataProvider::unsubscribe()
{
    currSubscriber = NULL;
}

void EditDataProvider::switchOffEditMode()
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode ();
    }
}

CursorShape EditDataProvider::getCursor(int objectID)
{
    if (currSubscriber) {
        currSubscriber->getCursor(objectID);
    }

    return CSOpenHand;
}

EditSubscriber* EditDataProvider::getCurrSubscriber()
{
    return currSubscriber;
}


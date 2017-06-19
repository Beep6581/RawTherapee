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
#include "../rtengine/icons.h"

ObjectMOBuffer::ObjectMOBuffer (EditDataProvider *dataProvider) : objectMap (nullptr), objectMode (OM_255), dataProvider (dataProvider) {}

ObjectMOBuffer::~ObjectMOBuffer()
{
    flush();
}


/* Upgrade or downgrade the objectModeType */
void ObjectMOBuffer::setObjectMode (ObjectMode newType)
{
    if (!objectMap) {
        objectMode = newType;
        return;
    }

    int w = objectMap->get_width ();
    int h = objectMap->get_height ();

    if (w && h) {
        switch (newType) {
            case (OM_255):
                if (objectMode == OM_65535) {
                    objectMap->unreference();
                    objectMap = Cairo::ImageSurface::create (Cairo::FORMAT_A8, w, h);
                }

                break;

            case (OM_65535):
                if (objectMode == OM_255) {
                    objectMap->unreference();
                    objectMap = Cairo::ImageSurface::create (Cairo::FORMAT_RGB16_565, w, h);
                }

                break;
        }
    }

    objectMode = newType;
}

void ObjectMOBuffer::flush()
{
    if (objectMap ) {
        objectMap.clear();
    }
}

EditSubscriber *ObjectMOBuffer::getEditSubscriber ()
{
    if (dataProvider) {
        return dataProvider->getCurrSubscriber();
    } else {
        return nullptr;
    }
}


// Resize buffers if they already exist
void ObjectMOBuffer::resize (int newWidth, int newHeight)
{
    if (!dataProvider) {
        return;
    }

    if (const auto currSubscriber = dataProvider->getCurrSubscriber ()) {
        if (currSubscriber->getEditingType() == ET_OBJECTS) {
            if (objectMap && (objectMap->get_width() != newWidth || objectMap->get_height() != newHeight)) {
                objectMap.clear();
            }

            if (!objectMap && newWidth > 0 && newHeight > 0) {
                objectMap = Cairo::ImageSurface::create (objectMode == OM_255 ? Cairo::FORMAT_A8 : Cairo::FORMAT_RGB16_565, newWidth, newHeight);
            }

        } else {
            flush();
        }
    } else {
        flush();
    }
}

int ObjectMOBuffer::getObjectID (const rtengine::Coord& location)
{
    int id = 0;

    if (!objectMap || location.x < 0 || location.y < 0 || location.x >= objectMap->get_width() || location.y >= objectMap->get_height()) {
        return -1;
    }

    if (objectMode == OM_255) {
        id = (unsigned char) (* ( objectMap->get_data() + location.y * objectMap->get_stride() + location.x ));
    } else {
        id = (unsigned short) (* ( objectMap->get_data() + location.y * objectMap->get_stride() + location.x ));
    }

    return id - 1;
}

bool ObjectMOBuffer::bufferCreated()
{
    EditSubscriber* subscriber;

    if (dataProvider && (subscriber = dataProvider->getCurrSubscriber())) {
        return subscriber->getEditingType() == ET_OBJECTS ? bool (objectMap) : false;
    }

    return false;
}

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

void Circle::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width ( getOuterLineWidth() );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);
        cr->stroke();
    }
}

void Circle::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
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

        cr->set_line_width ( innerLineWidth );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (filled && state != INSENSITIVE) {
            cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);

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

void Circle::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );
        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);

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

        cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, 0., 2.*rtengine::RT_PI);


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

void Arcellipse::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width ( getOuterLineWidth() );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);
//        double radius2_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius2)) : double (radius2);
        double scalx_ = scalx; //radius2_ / radius_;

        double scaly_ = scaly;
        double begang_ = begang;
        double endang_ = endang;
        double translax_ = translax;
        double translay_ = translay;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->save();
        cr->translate (translax_, translay_);

        cr->scale (scalx_, scaly_);
        cr->translate (- translax_ , - translay_);

        cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, begang_, endang_);

        cr->restore();
        cr->stroke();
    }
}



void Arcellipse::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
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

        cr->set_line_width ( innerLineWidth );

        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);
//        double radius2_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius2)) : double (radius2);

        double scalx_ = scalx; //radius2_ / radius_;

        double scaly_ = scaly;
        double begang_ = begang;
        double endang_ = endang;
        double translax_ = translax;
        double translay_ = translay;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        if (filled && state != INSENSITIVE) {
            cr->save();

            cr->translate (translax_, translay_);

            cr->scale (scalx_, scaly_);
            cr->translate (- translax_ , - translay_);

            cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, begang_, endang_);

            cr->restore();

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            cr->save();

            cr->translate (translax_ , translay_);

            cr->scale (scalx_, scaly_);
            cr->translate (- translax_ , - translay_);

            cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, begang_, endang_);

            cr->restore();

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




void Arcellipse::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );
        rtengine::Coord center_ = center;
        double radius_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius)) : double (radius);
//        double radius2_ = radiusInImageSpace ? coordSystem.scaleValueToCanvas (double (radius2)) : double (radius2);

        double scalx_ = scalx ; //radius2_ / radius_;
        double scaly_ = scaly;
        double begang_ = begang;
        double endang_ = endang;
        double translax_ = translax;
        double translay_ = translay;

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (center.x, center.y, center_.x, center_.y);
        } else if (datum == CLICKED_POINT) {
            center_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            center_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->save();

        cr->translate (translax_ , translay_);

        cr->scale (scalx_, scaly_);
        cr->translate (- translax_ , - translay_);
        cr->arc (center_.x + 0.5, center_.y + 0.5, radius_, begang_, endang_);

        cr->restore();

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

void Beziers::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;
        /*
                if (flags & F_AUTO_COLOR) {
                    color = getOuterLineColor();
                } else {
                    color = outerLineColor;
                }
        */
        //      cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        //   cr->set_line_width ( getOuterLineWidth() );
        cr->set_line_width ( innerLineWidth );

        rtengine::Coord begin_ = begin;
        rtengine::Coord inter_ = inter;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen (inter.x, inter.y, inter_.x, inter_.y);
            coordSystem.imageCoordToScreen (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->save();
        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->curve_to (begin_.x + 0.5, begin_.y + 0.5, inter_.x + 0.5, inter_.y + 0.5, end_.x + 0.5, end_.y + 0.5);
        cr->restore();
        cr->stroke();
    }
}

void Beziers::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && innerLineWidth > 0.) {
        /*
        if (state != INSENSITIVE) {
            RGBColor color;

            if (flags & F_AUTO_COLOR) {
                color = getInnerLineColor();
            } else {
                color = innerLineColor;
            }

            cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        }
        */
        cr->set_line_width (innerLineWidth);

        rtengine::Coord begin_ = begin;
        rtengine::Coord inter_ = inter;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToScreen (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToScreen (inter.x, inter.y, inter_.x, inter_.y);
            coordSystem.imageCoordToScreen (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        cr->save();

        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->curve_to (begin_.x + 0.5, begin_.y + 0.5, inter_.x + 0.5, inter_.y + 0.5, end_.x + 0.5, end_.y + 0.5);
        cr->restore();
        cr->stroke();
        /*
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
                */
    }
}

void Beziers::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );
        rtengine::Coord begin_ = begin;
        rtengine::Coord inter_ = inter;
        rtengine::Coord end_ = end;

        if (datum == IMAGE) {
            coordSystem.imageCoordToCropCanvas (begin.x, begin.y, begin_.x, begin_.y);
            coordSystem.imageCoordToCropCanvas (inter.x, inter.y, inter_.x, inter_.y);
            coordSystem.imageCoordToCropCanvas (end.x, end.y, end_.x, end_.y);
        } else if (datum == CLICKED_POINT) {
            begin_ += objectBuffer->getDataProvider()->posScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen;
            end_ += objectBuffer->getDataProvider()->posScreen;
        } else if (datum == CURSOR) {
            begin_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            inter_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            end_ += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
        }

        // setting the color to the objet's ID
        /*
        if (objectBuffer->getObjectMode() == OM_255) {
            cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
        } else {
            cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
        }
        */
        cr->save();

        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->curve_to (begin_.x + 0.5, begin_.y + 0.5, inter_.x + 0.5, inter_.y + 0.5, end_.x + 0.5, end_.y + 0.5);
        cr->restore();

        cr->stroke();
    }
}

void Line::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width ( getOuterLineWidth() );

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

        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to (end_.x + 0.5, end_.y + 0.5);
        cr->stroke();
    }
}

void Line::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
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

        cr->set_line_width (innerLineWidth);

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

        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to (end_.x + 0.5, end_.y + 0.5);

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

void Line::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );
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

        cr->move_to (begin_.x + 0.5, begin_.y + 0.5);
        cr->line_to (end_.x + 0.5, end_.y + 0.5);
        cr->stroke();
    }
}

void Polyline::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE && points.size() > 1) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width ( getOuterLineWidth() );

        rtengine::Coord currPos;

        for (unsigned int i = 0; i < points.size(); ++i) {
            currPos  = points.at (i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToScreen (points.at (i).x, points.at (i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += objectBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            }

            if (!i) {
                cr->move_to (currPos.x + 0.5, currPos.y + 0.5);
            } else {
                cr->line_to (currPos.x + 0.5, currPos.y + 0.5);
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

void Polyline::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
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

        cr->set_line_width ( innerLineWidth );

        if (filled && state != INSENSITIVE) {
            rtengine::Coord currPos;

            for (unsigned int i = 0; i < points.size(); ++i) {
                currPos  = points.at (i);

                if      (datum == IMAGE) {
                    coordSystem.imageCoordToScreen (points.at (i).x, points.at (i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += objectBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to (currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to (currPos.x + 0.5, currPos.y + 0.5);
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
                currPos  = points.at (i);

                if (datum == IMAGE) {
                    coordSystem.imageCoordToScreen (points.at (i).x, points.at (i).y, currPos.x, currPos.y);
                } else if (datum == CLICKED_POINT) {
                    currPos += objectBuffer->getDataProvider()->posScreen;
                } else if (datum == CURSOR) {
                    currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
                }

                if (!i) {
                    cr->move_to (currPos.x + 0.5, currPos.y + 0.5);
                } else {
                    cr->line_to (currPos.x + 0.5, currPos.y + 0.5);
                }
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

        for (unsigned int i = 0; i < points.size(); ++i) {
            cr->set_line_width ( getMouseOverLineWidth() );
            currPos  = points.at (i);

            if      (datum == IMAGE) {
                coordSystem.imageCoordToCropCanvas (points.at (i).x, points.at (i).y, currPos.x, currPos.y);
            } else if (datum == CLICKED_POINT) {
                currPos += objectBuffer->getDataProvider()->posScreen;
            } else if (datum == CURSOR) {
                currPos += objectBuffer->getDataProvider()->posScreen + objectBuffer->getDataProvider()->deltaScreen;
            }

            if (!i) {
                cr->move_to (currPos.x + 0.5, currPos.y + 0.5);
            } else {
                cr->line_to (currPos.x + 0.5, currPos.y + 0.5);
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
    }
}

void Rectangle::setXYWH (int left, int top, int width, int height)
{
    topLeft.set (left, top);
    bottomRight.set (left + width, top + height);
}

void Rectangle::setXYXY (int left, int top, int right, int bottom)
{
    topLeft.set (left, top);
    bottomRight.set (right, bottom);
}

void Rectangle::setXYWH (rtengine::Coord topLeft, rtengine::Coord widthHeight)
{
    this->topLeft = topLeft;
    this->bottomRight = topLeft + widthHeight;
}

void Rectangle::setXYXY (rtengine::Coord topLeft, rtengine::Coord bottomRight)
{
    this->topLeft = topLeft;
    this->bottomRight = bottomRight;
}

void Rectangle::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if ((flags & F_VISIBLE) && state != INSENSITIVE) {
        RGBColor color;

        if (flags & F_AUTO_COLOR) {
            color = getOuterLineColor();
        } else {
            color = outerLineColor;
        }

        cr->set_source_rgb (color.getR(), color.getG(), color.getB());
        cr->set_line_width ( getOuterLineWidth() );

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

        cr->rectangle (tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

        if (filled) {
            cr->fill_preserve();
            cr->stroke();
        } else {
            cr->stroke();
        }
    }
}

void Rectangle::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
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

        cr->set_line_width ( innerLineWidth );

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

        if (filled && state != INSENSITIVE) {
            cr->rectangle (tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

            if (innerLineWidth > 0.) {
                cr->fill_preserve();
                cr->stroke();
            } else {
                cr->fill();
            }
        } else if (innerLineWidth > 0.) {
            cr->rectangle (tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

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

void Rectangle::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        cr->set_line_width ( getMouseOverLineWidth() );

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

        cr->rectangle (tl.x + 0.5, tl.y + 0.5, br.x - tl.x, br.y - tl.y);

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

void OPIcon::drivenPointToRectangle (const rtengine::Coord &pos,
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

OPIcon::OPIcon (const Cairo::RefPtr<Cairo::ImageSurface> &normal,
                const Cairo::RefPtr<Cairo::ImageSurface> &active,
                const Cairo::RefPtr<Cairo::ImageSurface> &prelight,
                const Cairo::RefPtr<Cairo::ImageSurface> &dragged,
                const Cairo::RefPtr<Cairo::ImageSurface> &insensitive,
                DrivenPoint drivenPoint) :
    drivenPoint (drivenPoint)
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

OPIcon::OPIcon (Glib::ustring normalImage, Glib::ustring activeImage, Glib::ustring prelightImage,
                Glib::ustring  draggedImage, Glib::ustring insensitiveImage, DrivenPoint drivenPoint) : drivenPoint (drivenPoint)
{
    if (!normalImage.empty()) {
        normalImg = Cairo::ImageSurface::create_from_png ( rtengine::findIconAbsolutePath (normalImage) );
    }

    if (!prelightImage.empty()) {
        prelightImg = Cairo::ImageSurface::create_from_png ( rtengine::findIconAbsolutePath (prelightImage) );
    }

    if (!activeImage.empty()) {
        activeImg = Cairo::ImageSurface::create_from_png ( rtengine::findIconAbsolutePath (activeImage) );
    }

    if (!draggedImage.empty()) {
        draggedImg = Cairo::ImageSurface::create_from_png ( rtengine::findIconAbsolutePath (draggedImage) );
    }

    if (!insensitiveImage.empty()) {
        insensitiveImg = Cairo::ImageSurface::create_from_png ( rtengine::findIconAbsolutePath (insensitiveImage) );
    }
}

const Cairo::RefPtr<Cairo::ImageSurface> OPIcon::getNormalImg()
{
    return normalImg;
}
const Cairo::RefPtr<Cairo::ImageSurface> OPIcon::getPrelightImg()
{
    return prelightImg;
}
const Cairo::RefPtr<Cairo::ImageSurface> OPIcon::getActiveImg()
{
    return activeImg;
}
const Cairo::RefPtr<Cairo::ImageSurface> OPIcon::getDraggedImg()
{
    return draggedImg;
}
const Cairo::RefPtr<Cairo::ImageSurface> OPIcon::getInsensitiveImg()
{
    return insensitiveImg;
}

void OPIcon::drawImage (const Cairo::RefPtr<Cairo::ImageSurface> &img,
                        Cairo::RefPtr<Cairo::Context> &cr, ObjectMOBuffer *objectBuffer,
                        EditCoordSystem &coordSystem)
{
    int imgW = img->get_width();
    int imgH = img->get_height();

    rtengine::Coord pos;

    if (datum == IMAGE) {
        coordSystem.imageCoordToScreen (position.x, position.y, pos.x, pos.y);
    } else if (datum == CLICKED_POINT) {
        pos = position + objectBuffer->getDataProvider()->posScreen;
    } else if (datum == CURSOR)
        pos = position + objectBuffer->getDataProvider()->posScreen
              + objectBuffer->getDataProvider()->deltaScreen;

    rtengine::Coord tl, br; // Coordinate of the rectangle in the CropBuffer coordinate system
    drivenPointToRectangle (pos, tl, br, imgW, imgH);

    cr->set_source (img, tl.x, tl.y);
    cr->set_line_width (0.);
    cr->rectangle (tl.x, tl.y, imgW, imgH);
    cr->fill();
}

void OPIcon::drawMOImage (const Cairo::RefPtr<Cairo::ImageSurface> &img, Cairo::RefPtr<Cairo::Context> &cr,
                          unsigned short id, ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    // test of F_HOVERABLE has already been done

    int imgW = img->get_width();
    int imgH = img->get_height();

    rtengine::Coord pos;

    if (datum == IMAGE) {
        coordSystem.imageCoordToCropCanvas (position.x, position.y, pos.x, pos.y);
    } else if (datum == CLICKED_POINT) {
        pos = position + objectBuffer->getDataProvider()->posScreen;
    } else if (datum == CURSOR)
        pos = position + objectBuffer->getDataProvider()->posScreen
              + objectBuffer->getDataProvider()->deltaScreen;

    rtengine::Coord tl, br; // Coordinate of the rectangle in the CropBuffer coordinate system
    drivenPointToRectangle (pos, tl, br, imgW, imgH);

    // drawing the lower byte's value
    if (objectBuffer->getObjectMode() == OM_255) {
        cr->set_source_rgba (0., 0., 0., ((id + 1) & 0xFF) / 255.);
    } else {
        cr->set_source_rgba (0., 0., 0., (id + 1) / 65535.);
    }

    cr->set_line_width (0.);
    cr->rectangle (tl.x, tl.y, imgW, imgH);
    cr->fill();
}

void OPIcon::drawOuterGeometry (Cairo::RefPtr<Cairo::Context> &cr,
                                ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem) {}

void OPIcon::drawInnerGeometry (Cairo::RefPtr<Cairo::Context> &cr,
                                ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_VISIBLE) {
        // Here we will handle fall-back solutions

        State tmpState = state;  // can be updated through the successive test

        if (tmpState == INSENSITIVE) {
            if (!insensitiveImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawImage (insensitiveImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == DRAGGED) {
            if (!draggedImg) {
                tmpState = ACTIVE;
            } else {
                OPIcon::drawImage (draggedImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == ACTIVE) {
            if (!activeImg) {
                tmpState = PRELIGHT;
            } else {
                OPIcon::drawImage (activeImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == PRELIGHT) {
            if (!prelightImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawImage (prelightImg, cr, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == NORMAL && normalImg) {
            OPIcon::drawImage (normalImg, cr, objectBuffer, coordSystem);
        }
    }
}

void OPIcon::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, unsigned short id,
                              ObjectMOBuffer *objectBuffer, EditCoordSystem &coordSystem)
{
    if (flags & F_HOVERABLE) {
        // Here we will handle fallback solutions
        State tmpState = state;

        if (tmpState == INSENSITIVE) {
            if (!insensitiveImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawMOImage (insensitiveImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == DRAGGED) {
            if (!draggedImg) {
                tmpState = ACTIVE;
            } else {
                OPIcon::drawMOImage (draggedImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == ACTIVE) {
            if (!activeImg) {
                tmpState = PRELIGHT;
            } else {
                OPIcon::drawMOImage (activeImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == PRELIGHT) {
            if (!prelightImg) {
                tmpState = NORMAL;
            } else {
                OPIcon::drawMOImage (prelightImg, cr, id, objectBuffer, coordSystem);
                return;
            }
        }

        if (tmpState == NORMAL && normalImg) {
            OPIcon::drawMOImage (normalImg, cr, id, objectBuffer, coordSystem);
        }
    }
}

EditSubscriber::EditSubscriber (EditType editType) : ID (EUID_None), editingType (editType), bufferType (BT_SINGLEPLANE_FLOAT), provider (nullptr), action (ES_ACTION_NONE) {}

void EditSubscriber::setEditProvider (EditDataProvider *provider)
{
    this->provider = provider;
}

void EditSubscriber::setEditID (EditUniqueID ID, BufferType buffType)
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
        provider->subscribe (this);
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

BufferType EditSubscriber::getPipetteBufferType()
{
    return bufferType;
}

bool EditSubscriber::isDragging()
{
    return action == ES_ACTION_DRAGGING;
}

bool EditSubscriber::isPicking()
{
    return action == ES_ACTION_PICKING;
}

//--------------------------------------------------------------------------------------------------


EditDataProvider::EditDataProvider() : currSubscriber (nullptr), object (0), posScreen (-1, -1), posImage (-1, -1),
    deltaScreen (0, 0), deltaImage (0, 0), deltaPrevScreen (0, 0), deltaPrevImage (0, 0)
{
    pipetteVal[0] = pipetteVal[1] = pipetteVal[2] = 0.f;
}

void EditDataProvider::subscribe (EditSubscriber *subscriber)
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode();
    }

    currSubscriber = subscriber;
}

void EditDataProvider::unsubscribe()
{
    currSubscriber = nullptr;
}

void EditDataProvider::switchOffEditMode()
{
    if (currSubscriber) {
        currSubscriber->switchOffEditMode ();
    }
}

CursorShape EditDataProvider::getCursor (int objectID)
{
    if (currSubscriber) {
        currSubscriber->getCursor (objectID);
    }

    return CSOpenHand;
}

EditSubscriber* EditDataProvider::getCurrSubscriber()
{
    return currSubscriber;
}


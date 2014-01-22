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

RGBColor Geometry::getInnerLineColor () {
	RGBColor color;
	if (flags & AUTO_COLOR) {
		if      (state == NORMAL)   { color.setColor (1., 1., 1.); }  // White
		else if (state == PRELIGHT) { color.setColor (1., 1., 0.); }  // Orange
		else if (state == DRAGGED)  { color.setColor (1., 0., 0.); }  // Red
	}
	else { color = innerLineColor; }
	return color;
}

RGBColor Geometry::getOuterLineColor () {
	RGBColor color;
	if (flags & AUTO_COLOR) {
		/*
		if      (state == NORMAL)   { color.setColor (0., 0., 0.); }  // Black
		else if (state == PRELIGHT) { color.setColor (0., 0., 0.); }  // Black
		else if (state == DRAGGED)  { color.setColor (1., 0., 0.); }  // Black
		*/
		color.setColor (0., 0., 0.);  // Black
	}
	else { color = outerLineColor; }
	return color;
}

void Circle::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if (flags & ACTIVE) {
		cr->set_line_width( getOuterLineWidth() );
		Coord center_ = center;
		double radius_ = radiusInImageSpace ? coordSystem.scaleValueToImage(double(radius)) : double(radius);
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(center.x, center.y, center_.x, center_.y);
		}
		else if (datum == CLICKED_POINT) {
			center_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}
		cr->set_source_rgb (outerLineColor.getR(), outerLineColor.getG(), outerLineColor.getB());
		cr->arc(center_.x, center_.y, radius_, 0., 2.*M_PI);
	}
}

void Circle::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if (flags & ACTIVE) {
		Coord center_ = center;
		double radius_ = radiusInImageSpace ? coordSystem.scaleValueToImage(double(radius)) : double(radius);
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(center.x, center.y, center_.x, center_.y);
		}
		else if (datum == CLICKED_POINT) {
			center_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}
		cr->set_source_rgb (innerLineColor.getR(), innerLineColor.getG(), innerLineColor.getB());
		if (filled) {
			cr->arc(center_.x, center_.y, radius_, 0., 2.*M_PI);
			cr->set_line_width( innerLineWidth );
			if (innerLineWidth > 0.) {
				cr->fill_preserve();
				cr->stroke();
			}
			else
				cr->fill();
		}
		else if (innerLineWidth > 0.) {
			cr->arc(center_.x, center_.y, radius_, 0., 2.*M_PI);
			cr->stroke();
		}
	}
}

void Circle::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if (flags & ACTIVE) {
		cr->set_line_width( getMouseOverLineWidth() );
		Coord center_ = center;
		double radius_ = radiusInImageSpace ? coordSystem.scaleValueToImage(double(radius)) : double(radius);
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(center.x, center.y, center_.x, center_.y);
		}
		else if (datum == CLICKED_POINT) {
			center_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			center_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}

		// drawing the lower byte's value
		unsigned short a = (id+1) & 0xFF;
		cr->set_source_rgba (0.,0., 0., double(a)/255.);
		cr->arc(center_.x, center_.y, radius_, 0, 2.*M_PI);
		if (filled) {
			if (innerLineWidth > 0.) {
				cr->fill_preserve();
				cr->stroke();
			}
			else
				cr->fill();
		}
		else
			cr->stroke();

		// drawing the higher byte's value
		if (editBuffer->getObjectMode() == OM_65535) {
			a = (id+1)>>8;
			cr2->set_source_rgba (0.,0., 0., double(a)/255.);
			cr2->arc(center_.x, center_.y, radius_, 0, 2.*M_PI);
			if (filled) {
				if (innerLineWidth > 0.) {
					cr2->fill_preserve();
					cr2->stroke();
				}
				else
					cr2->fill();
			}
			else
				cr2->stroke();
		}
	}
}

void Line::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if (flags & ACTIVE) {
		cr->set_line_width( getOuterLineWidth() );
		Coord begin_ = begin;
		Coord end_ = end;
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(begin.x, begin.y, begin_.x, begin_.y);
			coordSystem.imageCoordToScreen(end.x, end.y, end_.x, end_.y);
		}
		else if (datum == CLICKED_POINT) {
			begin_ += editBuffer->getDataProvider()->posScreen;
			end_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
			end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}
		cr->set_source_rgb (outerLineColor.getR(), outerLineColor.getG(), outerLineColor.getB());
		cr->move_to(begin_.x, begin_.y);
		cr->line_to(end_.x, end_.y);
		cr->stroke();
	}
}

void Line::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if ((flags & ACTIVE) && innerLineWidth > 0.) {
		Coord begin_ = begin;
		Coord end_ = end;
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(begin.x, begin.y, begin_.x, begin_.y);
			coordSystem.imageCoordToScreen(end.x, end.y, end_.x, end_.y);
		}
		else if (datum == CLICKED_POINT) {
			begin_ += editBuffer->getDataProvider()->posScreen;
			end_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
			end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}
		cr->set_line_width(innerLineWidth);
		cr->set_source_rgb (innerLineColor.getR(), innerLineColor.getG(), innerLineColor.getB());
		cr->move_to(begin_.x, begin_.y);
		cr->line_to(end_.x, end_.y);
		cr->stroke();
	}
}

void Line::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if (flags & ACTIVE) {
		cr->set_line_width( getMouseOverLineWidth() );
		Coord begin_ = begin;
		Coord end_ = end;
		if (datum == IMAGE) {
			coordSystem.imageCoordToScreen(begin.x, begin.y, begin_.x, begin_.y);
			coordSystem.imageCoordToScreen(end.x, end.y, end_.x, end_.y);
		}
		else if (datum == CLICKED_POINT) {
			begin_ += editBuffer->getDataProvider()->posScreen;
			end_ += editBuffer->getDataProvider()->posScreen;
		}
		else if (datum == CURSOR) {
			begin_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
			end_ += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;
		}

		// drawing the lower byte's value
		unsigned short a = (id+1) & 0xFF;
		cr->set_source_rgba (0.,0., 0., double(a)/255.);
		cr->move_to(begin_.x, begin_.y);
		cr->line_to(end_.x, end_.y);
		cr->stroke();

		// drawing the higher byte's value
		if (editBuffer->getObjectMode() == OM_65535) {
			a = (id+1)>>8;
			cr2->set_source_rgba (0.,0., 0., double(a)/255.);
			cr2->move_to(begin_.x, begin_.y);
			cr2->line_to(end_.x, end_.y);
			cr2->stroke();
		}
	}
}

void Polyline::drawOuterGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if ((flags & ACTIVE) && points.size()>1) {
		cr->set_source_rgb (outerLineColor.getR(), outerLineColor.getG(), outerLineColor.getB());

		Coord currPos;
		for (unsigned int i=0; i<points.size(); ++i) {
			cr->set_line_width( getOuterLineWidth() );
			currPos  = points.at(i);

			if      (datum == IMAGE)         coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
			else if (datum == CLICKED_POINT) currPos += editBuffer->getDataProvider()->posScreen;
			else if (datum == CURSOR)        currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;

			if (!i) cr->move_to(currPos.x, currPos.y);
			else    cr->line_to(currPos.x, currPos.y);
		}
		if (filled) {
			cr->fill_preserve();
			cr->stroke();
		}
		else
			cr->fill();
	}
}

void Polyline::drawInnerGeometry(Cairo::RefPtr<Cairo::Context> &cr, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if ((flags & ACTIVE) && points.size()>1) {
		cr->set_source_rgb (innerLineColor.getR(), innerLineColor.getG(), innerLineColor.getB());

		if (filled) {
			Coord currPos;
			for (unsigned int i=0; i<points.size(); ++i) {
				cr->set_line_width( getOuterLineWidth() );
				currPos  = points.at(i);

				if      (datum == IMAGE)         coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
				else if (datum == CLICKED_POINT) currPos += editBuffer->getDataProvider()->posScreen;
				else if (datum == CURSOR)        currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;

				if (!i) cr->move_to(currPos.x, currPos.y);
				else    cr->line_to(currPos.x, currPos.y);
			}

			if (innerLineWidth > 0.) {
				cr->fill_preserve();
				cr->stroke();
			}
			else
				cr->fill();
		}
		else if (innerLineWidth > 0.) {
			Coord currPos;
			for (unsigned int i=0; i<points.size(); ++i) {
				cr->set_line_width( getOuterLineWidth() );
				currPos  = points.at(i);

				if      (datum == IMAGE)         coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
				else if (datum == CLICKED_POINT) currPos += editBuffer->getDataProvider()->posScreen;
				else if (datum == CURSOR)        currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;

				if (!i) cr->move_to(currPos.x, currPos.y);
				else    cr->line_to(currPos.x, currPos.y);
			}
			cr->fill();
		}
	}
}

void Polyline::drawToMOChannel (Cairo::RefPtr<Cairo::Context> &cr, Cairo::RefPtr<Cairo::Context> &cr2, unsigned short id, rtengine::EditBuffer *editBuffer, EditCoordSystem &coordSystem) {
	if ((flags & ACTIVE) && points.size()>1) {
		Coord currPos;

		// drawing the lower byte's value
		unsigned short a = (id+1) & 0xFF;
		cr->set_source_rgba (0.,0., 0., double(a)/255.);
		for (unsigned int i=0; i<points.size(); ++i) {
			cr->set_line_width( getOuterLineWidth() );
			currPos  = points.at(i);

			if      (datum == IMAGE)         coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
			else if (datum == CLICKED_POINT) currPos += editBuffer->getDataProvider()->posScreen;
			else if (datum == CURSOR)        currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;

			if (!i) cr->move_to(currPos.x, currPos.y);
			else    cr->line_to(currPos.x, currPos.y);
		}
		if (filled) {
			if (innerLineWidth > 0.) {
				cr->fill_preserve();
				cr->stroke();
			}
			else
				cr->fill();
		}
		else
			cr->stroke();

		// drawing the higher byte's value
		if (editBuffer->getObjectMode() == OM_65535) {
			a = (id+1)>>8;
			cr2->set_source_rgba (0.,0., 0., double(a)/255.);
			for (unsigned int i=0; i<points.size(); ++i) {
				cr2->set_line_width( getOuterLineWidth() );
				currPos  = points.at(i);

				if      (datum == IMAGE)         coordSystem.imageCoordToScreen(points.at(i).x, points.at(i).y, currPos.x, currPos.y);
				else if (datum == CLICKED_POINT) currPos += editBuffer->getDataProvider()->posScreen;
				else if (datum == CURSOR)        currPos += editBuffer->getDataProvider()->posScreen + editBuffer->getDataProvider()->deltaPrevScreen;

				if (!i) cr2->move_to(currPos.x, currPos.y);
				else    cr2->line_to(currPos.x, currPos.y);
			}
			if (filled) {
				if (innerLineWidth > 0.) {
					cr2->fill_preserve();
					cr2->stroke();
				}
				else
					cr2->fill();
			}
			else
				cr2->stroke();
		}
	}
}

EditSubscriber::EditSubscriber () : ID(EUID_None), editingType(ET_PIPETTE), bufferType(BT_SINGLEPLANE_FLOAT), provider(NULL) {}

void EditSubscriber::setEditProvider(EditDataProvider *provider) {
	this->provider = provider;
}

void EditSubscriber::setEditID(EditUniqueID ID, BufferType buffType) {
	this->ID = ID;
	bufferType = buffType;
}

bool EditSubscriber::isCurrentSubscriber() {
	//if (provider && provider->getCurrSubscriber())
	//	return provider->getCurrSubscriber()->getEditID() == ID;

	if (provider)
		return provider->getCurrSubscriber() == this;

	return false;
}

void EditSubscriber::subscribe() {
	if (provider)
		provider->subscribe(this);
}

void EditSubscriber::unsubscribe() {
	if (provider)
		provider->unsubscribe();
}

void EditSubscriber::switchOffEditMode() {
	unsubscribe();
}

EditUniqueID EditSubscriber::getEditID() {
	return ID;
}

EditType EditSubscriber::getEditingType() {
	return editingType;
}

BufferType EditSubscriber::getEditBufferType() {
	return bufferType;
}


//--------------------------------------------------------------------------------------------------


EditDataProvider::EditDataProvider() : currSubscriber(NULL), object(0), posScreen(-1,-1), posImage(-1,-1),
									   deltaScreen(0,0), deltaImage(0,0), deltaPrevScreen(0,0), deltaPrevImage(0,0) {
	pipetteVal[0] = pipetteVal[1] = pipetteVal[2] = 0.f;
}

void EditDataProvider::subscribe(EditSubscriber *subscriber) {
	currSubscriber = subscriber;
}

void EditDataProvider::unsubscribe() {
	currSubscriber = NULL;
}

void EditDataProvider::switchOffEditMode() {
	if (currSubscriber)
		currSubscriber->switchOffEditMode ();
}

CursorShape EditDataProvider::getCursor(int objectID) {
	return CSOpenHand;
}

EditSubscriber* EditDataProvider::getCurrSubscriber() {
	return currSubscriber;
}


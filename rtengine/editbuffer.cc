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

#include "editbuffer.h"

namespace rtengine {

EditBuffer::EditBuffer(::EditDataProvider *dataProvider) :
	objectMap(NULL), objectMap2(NULL), objectMode(OM_255), dataProvider(dataProvider),
	imgFloatBuffer(NULL), LabBuffer(NULL), singlePlaneBuffer() {}

EditBuffer::~EditBuffer() {
	flush();
	#ifndef NDEBUG
	imgFloatBuffer = (Imagefloat*)(0xbaadf00d);
	#endif
	#ifndef NDEBUG
	LabBuffer = (LabImage*)(0xbaadf00d);
	#endif
}

void EditBuffer::createBuffer(int width, int height) {
	resize (width, height);
}

void EditBuffer::flush() {
	if (imgFloatBuffer) {
		delete imgFloatBuffer;
		imgFloatBuffer = NULL;
	}
	if (LabBuffer) {
		delete LabBuffer;
		LabBuffer = NULL;
	}
	singlePlaneBuffer.flushData();
}

/* Upgrade or downgrade the objectModeType; we're assuming that objectMap has been allocated */
void EditBuffer::setObjectMode(ObjectMode newType) {
	switch (newType) {
	case (OM_255):
		if (objectMap2) {
			objectMap2->unreference();
		}
		objectMode = OM_255;
		break;
	case (OM_65535):
		if (!objectMap2) {
			objectMap2 = Cairo::ImageSurface::create(Cairo::FORMAT_A8, objectMap->get_width(), objectMap->get_height());
		}
		objectMode = OM_65535;
		break;
	}
}

EditUniqueID EditBuffer::getEditID() {
	if (dataProvider && dataProvider->getCurrSubscriber()) {
		return dataProvider->getCurrSubscriber()->getEditID();
	}
	else return EUID_None;
}

// Resize buffers if they already exist
void EditBuffer::resize(int newWidth, int newHeight) {
	EditSubscriber* subscriber = NULL;
	if (dataProvider && (subscriber = dataProvider->getCurrSubscriber())) {
		if (subscriber->getEditingType() == ET_OBJECTS) {
			if (objectMap && (objectMap->get_width() != newWidth || objectMap->get_height() != newHeight))
				objectMap->unreference();

			if (!objectMap) {
				objectMap = Cairo::ImageSurface::create(Cairo::FORMAT_A8, newWidth, newHeight);
			}
			if (objectMode==OM_65535) {
				if (objectMap2) {
					if (objectMap2->get_width() != newWidth || objectMap2->get_height() != newHeight) {
						objectMap2->unreference();
					}
				}
				if (!objectMap2) {
					objectMap2 = Cairo::ImageSurface::create(Cairo::FORMAT_A8, newWidth, newHeight);
				}
			}
			// OM_255 -> deleting objectMap2, if any
			else if (objectMap2)
				objectMap2->unreference();

			// Should never happen!
			if (imgFloatBuffer) {
				delete imgFloatBuffer;
				imgFloatBuffer = NULL;
			}
			if (LabBuffer) {
				delete LabBuffer;
				LabBuffer = NULL;
			}
			if (singlePlaneBuffer.data) {
				singlePlaneBuffer.allocate(0,0);
			}
		}

		if (subscriber->getEditingType() == ET_PIPETTE) {
			if (subscriber->getEditBufferType() == BT_IMAGEFLOAT) {
				if (!imgFloatBuffer)
					imgFloatBuffer = new Imagefloat(newWidth, newHeight);
				else
					imgFloatBuffer->allocate(newWidth, newHeight);
			}
			else if (imgFloatBuffer) {
				delete imgFloatBuffer;
				imgFloatBuffer = NULL;
			}

			if (subscriber->getEditBufferType() == BT_LABIMAGE) {
				if (LabBuffer && (LabBuffer->W != newWidth && LabBuffer->H != newHeight)) {
					delete LabBuffer;
					LabBuffer = NULL;
				}
				if (!LabBuffer)
					LabBuffer = new LabImage(newWidth, newHeight);
			}
			else if (LabBuffer) {
				delete LabBuffer;
				LabBuffer = NULL;
			}

			if (subscriber->getEditBufferType() == BT_SINGLEPLANE_FLOAT) {
				singlePlaneBuffer.allocate(newWidth, newHeight);
			}
			else if (singlePlaneBuffer.data)
				singlePlaneBuffer.allocate(0,0);

			// Should never happen!
			if (objectMap ) objectMap->unreference();
			if (objectMap2) objectMap2->unreference();
		}
	}
}

bool EditBuffer::bufferCreated() {
	if (dataProvider && dataProvider->getCurrSubscriber()) {
		switch (dataProvider->getCurrSubscriber()->getEditBufferType()) {
		case (BT_IMAGEFLOAT):
			return imgFloatBuffer != NULL;
		case (BT_LABIMAGE):
			return LabBuffer != NULL;
		case (BT_SINGLEPLANE_FLOAT):
			return singlePlaneBuffer.data != NULL;
		}
	}
	return false;
}

unsigned short EditBuffer::getObjectID(const Coord& location) {
	unsigned short id = 0;

	if (objectMap)
		id = (unsigned short)(*( objectMap->get_data() + location.y * objectMap->get_stride() + location.x ));

	if (objectMap2)
		id |= (unsigned short)(*( objectMap->get_data() + location.y * objectMap->get_stride() + location.x )) << 8;

	return id;
}

void EditBuffer::getPipetteData(float* v, int x, int y, int squareSize) {
	if (dataProvider && dataProvider->getCurrSubscriber()) {
		switch (dataProvider->getCurrSubscriber()->getEditBufferType()) {
		case (BT_IMAGEFLOAT):
			if (imgFloatBuffer)
				imgFloatBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize, 0);
			else
				v[0] = v[1] = v[2] = -1.f;
			break;
		case (BT_LABIMAGE):
			if (LabBuffer)
				LabBuffer->getPipetteData(v[0], v[1], v[2], x, y, squareSize);
			else
				v[0] = v[1] = v[2] = -1.f;
			break;
		case (BT_SINGLEPLANE_FLOAT):
			if (singlePlaneBuffer.data != NULL) {
				singlePlaneBuffer.getPipetteData(v[0], x, y, squareSize, 0);
				v[1] = v[2] = -1.f;
			}
			else
				v[0] = v[1] = v[2] = -1.f;
		}
	}
}

}

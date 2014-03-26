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
#ifndef _MYDIAGONALCURVE_
#define _MYDIAGONALCURVE_

#include <gtkmm.h>
#include <vector>
#include "curvelistener.h"
#include "cursormanager.h"
#include "mycurve.h"
#include "../rtengine/LUT.h"


// For compatibility and simplicity reason, order shouldn't change, and must be identical to the order specified in the curveType widget
enum DiagonalCurveType {
	DCT_Empty = -1,		// Also used for identity curves
	DCT_Linear,			// 0
	DCT_Spline,			// 1
	DCT_Parametric,		// 2
	DCT_NURBS,			// 3
	// Insert new curve type above this line
	DCT_Unchanged		// Must remain the last of the enum
};

class DiagonalCurveDescr {

	public:
		DiagonalCurveType type;
		std::vector<double> x, y;   // in case of parametric curves the curve parameters are stored in vector x. In other cases these vectors store the coordinates of the bullets.
};

class MyDiagonalCurve : public MyCurve {

	protected:
		DiagonalCurveDescr curve;
		int grab_point;		// the point that the user is moving
		int closest_point;	// the point that is the closest from the cursor
		int lit_point;		// the point that is lit when the cursor is near it
		double clampedX;	// clamped grabbed point X coordinates in the [0;1] range
		double clampedY;	// clamped grabbed point Y coordinates in the [0;1] range
		double deltaX;		// signed X distance of the cursor between two consecutive MOTION_NOTIFY
		double deltaY;		// signed Y distance of the cursor between two consecutive MOTION_NOTIFY
		double distanceX;	// X distance from the cursor to the closest point
		double distanceY;	// Y distance from the cursor to the closest point
		double ugpX;		// unclamped grabbed point X coordinate in the graph
		double ugpY;		// unclamped grabbed point Y coordinate in the graph
		int activeParam;
		unsigned int* bghist;	// histogram values
		bool bghistvalid;
		void draw (int handle);
		void interpolate ();
		void getCursorPosition(Gdk::EventType evType, bool isHint, int evX, int evY, Gdk::ModifierType modifierKey);
		void findClosestPoint();
		CursorShape motionNotify(CursorShape type, double minDistanceX, double minDistanceY, int num);
		std::vector<double> get_vector (int veclen);

	public:
		MyDiagonalCurve ();
		~MyDiagonalCurve ();
		std::vector<double> getPoints ();
		void setPoints (const std::vector<double>& p);
		void setType (DiagonalCurveType t);
		bool handleEvents (GdkEvent* event);
		void setActiveParam (int ac);
		void reset (double identityValue=0.5);
		void updateBackgroundHistogram (LUTu & hist);

		void pipetteMouseOver (EditDataProvider *provider, int modifierKey);
		void pipetteButton1Pressed(EditDataProvider *provider, int modifierKey);
		void pipetteButton1Released(EditDataProvider *provider);
		void pipetteDrag(EditDataProvider *provider, int modifierKey);
};

#endif

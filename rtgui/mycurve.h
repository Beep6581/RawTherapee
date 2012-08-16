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
#ifndef _MYCURVE_
#define _MYCURVE_

#include <gtkmm.h>
#include <vector>
#include "curvelistener.h"
#include "cursormanager.h"
#include "coloredbar.h"
#include "../rtengine/LUT.h"
#include "guiutils.h"
#include "options.h"

#define RADIUS			3	/* radius of the control points */
#define CBAR_WIDTH_STD	13	/* width of the colored bar (border included) for standard themes */
#define CBAR_WIDTH_SLIM	10	/* width of the colored bar (border included) for slim themes */
#define CBAR_MARGIN		2	/* spacing between the colored bar and the graph */
#define SQUARE			2	/* half length of the square shape of the tangent handles */
#define MIN_DISTANCE	5	/* min distance between control points */
#define GRAPH_SIZE		200 /* size of the curve editor graphic */

// For compatibility and simplicity reason, order shouldn't change, and must be identical to the order specified in the curveType widget
enum CurveType {
	CT_Flat,
	CT_Diagonal
};

enum SnapToType {
	ST_None,
	ST_Identity,	// Point snapped to the identity curve
	ST_Neighbors	// Point snapped to the neighbor points
};

enum ResizeState {
	RS_Pending = 1,	// Resize has to occure
	RS_Done    = 2,	// Resize has been done
	RS_Force   = 4	// Resize has to occure even without CONFIGURE event
};

class MyCurveIdleHelper;

class MyCurve : public Gtk::DrawingArea, public BackBuffer, public ColorCaller {

	friend class MyCurveIdleHelper;

	protected:
		CurveListener* listener;
		ColoredBar *leftBar;
		ColoredBar *bottomBar;
		CursorShape cursor_type;
		int graphX, graphY, graphW, graphH; // dimensions of the graphic area, excluding surrounding space for the points of for the colored bar
		int prevGraphW, prevGraphH;         // previous inner width and height of the editor
		Gdk::ModifierType mod_type;
		int cursorX;		// X coordinate in the graph of the cursor
		int cursorY;		// Y coordinate in the graph of the cursor
		std::vector< Point<float> > point;
		std::vector< Point<float> > upoint;
		std::vector< Point<float> > lpoint;
		bool buttonPressed;
		/*
		 * snapToElmt must be interpreted like this:
		 * -100     : no element (default)
		 * -3       : maximum value
		 * -2       : identity value
		 * -1       : minimum value
		 * [0;1000[ : control point that it's snapped to
		 * >=1000   : moved control point which is snapped to to the line made by it neighbors
		 */
		int snapToElmt;
		bool snapTo;
		double snapToMinDist;
		double snapToVal;
		MyCurveIdleHelper* mcih;
		enum ResizeState sized;

	    virtual std::vector<double> get_vector (int veclen) = 0;
		int getGraphMinSize() { return GRAPH_SIZE + RADIUS + 1; }
		bool snapCoordinate(double testedVal, double realVal);

		// return value = new requested height
		int calcDimensions ();

	public:
		MyCurve ();
		~MyCurve ();

		void setCurveListener (CurveListener* cl) { listener = cl; }
		void setColoredBar (ColoredBar *left, ColoredBar *bottom);
		void notifyListener ();
		void updateBackgroundHistogram (LUTu & hist) {return;} ;
		void forceResize() { sized = RS_Force; }
	    void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
		virtual std::vector<double> getPoints () = 0;
		virtual void setPoints (const std::vector<double>& p) = 0;
		virtual bool handleEvents (GdkEvent* event) = 0;
		virtual void reset () = 0;

		static int getBarWidth() { return options.slimUI ? CBAR_WIDTH_SLIM : CBAR_WIDTH_STD; }
};

class MyCurveIdleHelper {
	public:
		MyCurve* myCurve;
		bool destroyed;
		int pending;

		void clearPixmap () { myCurve->setDirty(true); }
};

#endif

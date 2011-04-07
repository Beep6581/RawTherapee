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
#include <curvelistener.h>
#include <cursormanager.h>
#include <colorprovider.h>
#include <LUT.h>

#define RADIUS			3	/* radius of the control points */
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

class MyCurve : public Gtk::DrawingArea {

	friend class MyCurveIdleHelper;

	protected:
		CurveListener* listener;
		ColorProvider* colorProvider;
		CursorShape cursor_type;
		Glib::RefPtr<Gdk::Pixmap> pixmap;
		int innerWidth;		// inner width of the editor, allocated by the system
		int innerHeight;	// inner height of the editor, allocated by the system
		int prevInnerHeight;// previous inner height of the editor
		Gdk::ModifierType mod_type;
		int cursorX;		// X coordinate in the graph of the cursor
		int cursorY;		// Y coordinate in the graph of the cursor
		std::vector<Gdk::Point> point;
		std::vector<Gdk::Point> upoint;
		std::vector<Gdk::Point> lpoint;
		bool buttonPressed;
		enum SnapToType snapTo;
		MyCurveIdleHelper* mcih;
		enum ResizeState sized;

	    virtual std::vector<double> get_vector (int veclen) = 0;
		int getGraphMinSize() { return GRAPH_SIZE + RADIUS + 1; }

	public:
		MyCurve ();
		~MyCurve ();

		void setCurveListener (CurveListener* cl) { listener = cl; }
		void setColorProvider (ColorProvider* cp) { colorProvider = cp; }
		void notifyListener ();
		void updateBackgroundHistogram (LUTu & hist) {return;} ;
		void forceResize() { sized = RS_Force; }
		virtual std::vector<double> getPoints () = 0;
		virtual void setPoints (const std::vector<double>& p) = 0;
		virtual bool handleEvents (GdkEvent* event) = 0;
		virtual void reset () = 0;
};

class MyCurveIdleHelper {
	public:
		MyCurve* myCurve;
		bool destroyed;
		int pending;

		void clearPixmap () { myCurve->pixmap.clear (); }

};

#endif

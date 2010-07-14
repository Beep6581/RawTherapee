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

#define RADIUS			3	/* radius of the control points. Assuming that the center of the spot is in the center of the pixel, the real RADIUS will be this value +0.5 */
#define MIN_DISTANCE	8	/* min distance between control points */
#define GRAPH_SIZE		200 /* size of the curve editor graphic */

// For compatibility and simplicity reason, order shouldn't change, and must be identical to the order specified in the curveType widget
enum CurveType {
	Empty = -1,
	Linear,			// 0
	Spline,			// 1
	Parametric,		// 2
	NURBS,			// 3
	// Insert new curve type above this line
	Unchanged		// Must remain the last of the enum
};

class CurveDescr {

    public:
        CurveType type;
        std::vector<double> x, y;   // in case of parametric curves the curve parameters are stored in vector x. In other cases these vectors store the coordinates of the bullets.
};

class MyCurve;
struct MyCurveIdleHelper {
    MyCurve* myCurve;
    bool destroyed;
    int pending;
};

class MyCurve : public Gtk::DrawingArea {

        friend int mchistupdate (void* data);
    
    protected:
        CurveListener* listener;
        CurveDescr curve;
        CursorShape cursor_type;
        Glib::RefPtr<Gdk::Pixmap> pixmap;
        int height;
        int grab_point;
        int lit_point;
        int last;
        std::vector<Gdk::Point> point;
        std::vector<Gdk::Point> upoint;
        std::vector<Gdk::Point> lpoint;
        int activeParam;
        unsigned int bghist[256];
        bool bghistvalid;
        MyCurveIdleHelper* mcih;
    	int cursor_x, cursor_y;
    	double ugp_x, ugp_y;  // unclamped grabbed point coordinates
        
        void draw (int width, int height, int handle);
        void interpolate (int width, int height);
        std::vector<double> get_vector (int veclen);

    public:
        MyCurve ();
        ~MyCurve ();
        
        void setCurveListener (CurveListener* cl) { listener = cl; }
        std::vector<double> getPoints ();
        void setPoints (const std::vector<double>& p);
        void setType (CurveType t);
        bool handleEvents (GdkEvent* event);
        void notifyListener ();
        void setActiveParam (int ac);
        void updateBackgroundHistogram (unsigned int* hist);
        void reset ();
};

#endif

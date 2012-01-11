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
#include <glib.h>
#include <glib/gstdio.h>
#include "curves.h"
#include <cmath>
#include <vector>
#include "mytime.h"
#include <cstring>

#include <gtkmm.h>
#include <fstream>

namespace rtengine {

FlatCurve::FlatCurve (const std::vector<double>& p, bool isPeriodic, int poly_pn) : leftTangent(NULL), rightTangent(NULL) {

    ppn = poly_pn > 65500 ? 65500 : poly_pn;
    poly_x.clear();
    poly_y.clear();

    kind = FCT_Empty;
    periodic = isPeriodic;

    bool identity = true;

    if (p.size()>4) {
        kind = (FlatCurveType)p[0];
        if (kind==FCT_MinMaxCPoints) {
            int oneMorePoint = periodic ? 1:0;
            N = (p.size()-1)/4;
            x = new double[N+oneMorePoint];
            y = new double[N+oneMorePoint];
            leftTangent = new double[N+oneMorePoint];
            rightTangent = new double[N+oneMorePoint];
            int ix = 1;
            for (int i=0; i<N; i++) {
                x[i] = p[ix++];
                y[i] = p[ix++];
                leftTangent[i] = p[ix++];
                rightTangent[i] = p[ix++];
                if (y[i] != 0.5f)
                	identity = false;
            }
            // The first point is copied to the end of the point list, to handle the curve periodicity
            if (periodic) {
                x[N] = p[1]+1.0;
                y[N] = p[2];
                leftTangent[N] = p[3];
                rightTangent[N] = p[4];
            }
            else {
                N--;
            }
            if (!identity && N > 0+(periodic?1:0) ) {
                CtrlPoints_set ();
                fillHash();
            }
        }
        /*else if (kind==FCT_Parametric) {
        }*/
        if (identity)
            kind = FCT_Empty;
    }
}

FlatCurve::~FlatCurve () {

    delete [] x;
    delete [] y;
    delete [] leftTangent;
    delete [] rightTangent;
    delete [] ypp;
    poly_x.clear();
    poly_y.clear();
}

void FlatCurve::CtrlPoints_set () {

	int nbSubCurvesPoints = N*6;

    std::vector<double> sc_x(nbSubCurvesPoints);  // X sub-curve points (  XP0,XP1,XP2,  XP2,XP3,XP4,  ...)
    std::vector<double> sc_y(nbSubCurvesPoints);  // Y sub-curve points (  YP0,YP1,YP2,  YP2,YP3,YP4,  ...)
    std::vector<double> sc_length(N*2);           // Length of the subcurves
    std::vector<bool> sc_isLinear(N*2);           // true if the subcurve is linear
    double total_length=0.;

	// Create the list of Bezier sub-curves
    // CtrlPoints_set is called if N > 1

    unsigned int j = 0;
    unsigned int k = 0;

    for (int i = 0; i < N;) {
        double length;
        double dx;
        double dy;
        double xp1, xp2, yp2, xp3;
        bool startLinear, endLinear;

        startLinear = (rightTangent[i]   == 0.) || (y[i] == y[i+1]);
        endLinear   = (leftTangent [i+1] == 0.) || (y[i] == y[i+1]);

        if (startLinear && endLinear) {
        	// line shape
        	sc_x[j]   = x[i];
        	sc_y[j++] = y[i];
        	sc_x[j]   = x[i+1];
        	sc_y[j] = y[i+1];
        	sc_isLinear[k] = true;
        	i++;

    		dx = sc_x[j] - sc_x[j-1];
    		dy = sc_y[j] - sc_y[j-1];
    		length = sqrt(dx*dx + dy*dy);
    		j++;

    		// Storing the length of all sub-curves and the total length (to have a better distribution
    		// of the points along the curve)
    	    sc_length[k++] = length;
    	    total_length += length;
        }
        else {
        	if (startLinear) {
        		xp1 = x[i];
        	}
        	else {
        		//xp1 = (xp4 - xp0) * rightTangent0 + xp0;
        		xp1 = (x[i+1] - x[i]) * rightTangent[i] + x[i];
        	}
        	if (endLinear) {
        		xp3 = x[i+1];
        	}
        	else {
				//xp3 = (xp0 - xp4]) * leftTangent4 + xp4;
				xp3 = (x[i] - x[i+1]) * leftTangent[i+1] + x[i+1];
        	}

            xp2 = (xp1 + xp3) / 2.0;
            yp2 = (y[i] + y[i+1]) / 2.0;

            if (rightTangent[i]+leftTangent[i+1] > 1.0) {	// also means that start and end are not linear
            	xp1 = xp3 = xp2;
            }

        	if (startLinear) {
        		// Point 0, 2
            	sc_x[j]   = x[i];
            	sc_y[j++] = y[i];
            	sc_x[j]   = xp2;
            	sc_y[j]   = yp2;
            	sc_isLinear[k] = true;

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length = sqrt(dx*dx + dy*dy);
        		j++;

        		// Storing the length of all sub-curves and the total length (to have a better distribution
        		// of the points along the curve)
        	    sc_length[k++] = length;
        	    total_length += length;
        	}
        	else {
        		// Point 0, 1, 2
            	sc_x[j]   = x[i];
            	sc_y[j++] = y[i];
            	sc_x[j]   = xp1;
            	sc_y[j]   = y[i];

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length = sqrt(dx*dx + dy*dy);
        		j++;

            	sc_x[j] = xp2;
            	sc_y[j] = yp2;
            	sc_isLinear[k] = false;

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length += sqrt(dx*dx + dy*dy);
        		j++;

        		// Storing the length of all sub-curves and the total length (to have a better distribution
        		// of the points along the curve)
        	    sc_length[k++] = length;
        	    total_length += length;
        	}
        	if (endLinear) {
        		// Point 2, 4
            	sc_x[j]   = xp2;
            	sc_y[j++] = yp2;
            	sc_x[j]   = x[i+1];
            	sc_y[j]   = y[i+1];
            	sc_isLinear[k] = true;

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length = sqrt(dx*dx + dy*dy);
        		j++;

        		// Storing the length of all sub-curves and the total length (to have a better distribution
        		// of the points along the curve)
        	    sc_length[k++] = length;
        	    total_length += length;
        	}
        	else {
        		// Point 2, 3, 4
            	sc_x[j]   = xp2;
            	sc_y[j++] = yp2;
            	sc_x[j]   = xp3;
            	sc_y[j]   = y[i+1];

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length = sqrt(dx*dx + dy*dy);
        		j++;

            	sc_x[j] = x[i+1];
            	sc_y[j] = y[i+1];
            	sc_isLinear[k] = false;

        		dx = sc_x[j] - sc_x[j-1];
        		dy = sc_y[j] - sc_y[j-1];
        		length += sqrt(dx*dx + dy*dy);
        		j++;

        		// Storing the length of all sub-curves and the total length (to have a better distribution
        		// of the points along the curve)
        	    sc_length[k++] = length;
        	    total_length += length;
        	}
        	i++;
        }
    }

    poly_x.clear();
   	poly_y.clear();
    j = 0;

    // adding an initial horizontal line if necessary
    if (!periodic && sc_x[j] != 0.) {
		poly_x.push_back(0.);
		poly_y.push_back(sc_y[j]);
	}

	// the first point of the curves
	poly_x.push_back(sc_x[j]);
	poly_y.push_back(sc_y[j]);

	firstPointIncluded = false;

	// create the polyline with the number of points adapted to the X range of the sub-curve
	for (unsigned int i=0; i < k; i++) {
		if (sc_isLinear[i]) {
			j++; // skip the first point
			poly_x.push_back(sc_x[j]);
			poly_y.push_back(sc_y[j++]);
		}
	else {
			nbr_points = (int)(((double)(ppn) * sc_length[i] )/ total_length);
			if (nbr_points<0){
				for(unsigned int it=0;it < sc_x.size(); it+=3) printf("sc_length[%d/3]=%f \n",it,sc_length[it/3]);
				printf("Flat curve: error detected!\n i=%d periodic=%d nbr_points=%d ppn=%d N=%d sc_length[i/3]=%f total_length=%f",i,periodic,nbr_points,ppn,N,sc_length[i/3],total_length);
				exit(0);
			}
			// increment along the curve, not along the X axis
			increment = 1.0 / (double)(nbr_points-1);
			x1 = sc_x[j]; y1 = sc_y[j++];
			x2 = sc_x[j]; y2 = sc_y[j++];
			x3 = sc_x[j]; y3 = sc_y[j++];
			AddPolygons ();
	}
    }

    // adding the final horizontal segment, always (see under)
   	poly_x.push_back(3.0);		// 3.0 is a hack for optimization purpose of the getVal method (the last value has to be beyond the normal range)
	poly_y.push_back(sc_y[j-1]);

    /*
    // Checking the values
    Glib::ustring fname = "Curve.xyz"; // TopSolid'Design "plot" file format
	std::ofstream f (fname.c_str());
	f << "$" << std::endl;;
	for (unsigned int iter = 0; iter < poly_x.size(); iter++) {
		f << poly_x[iter] << ", " << poly_y[iter] << ", 0." << std::endl;;
	}
	f << "$" << std::endl;;
	f.close ();
	*/
}

double FlatCurve::getVal (double t) {

    switch (kind) {

    case FCT_MinMaxCPoints : {
        /* To be updated if we have to handle non periodic flat curves
        // values under and over the first and last point
        if (t>x[N-1])
            return y[N-1];
        else if (t<x[0])
            return y[0];
        */

        // magic to handle curve periodicity : we look above the 1.0 bound for the value
        if (t < poly_x[0]) t += 1.0;

        // do a binary search for the right interval:
        int k_lo = 0, k_hi = poly_x.size() - 1;
        while (k_hi - k_lo > 1){
            int k = (k_hi + k_lo) / 2;
            if (poly_x[k] > t)
                k_hi = k;
            else
                k_lo = k;
        }

        double dx = poly_x[k_hi] - poly_x[k_lo];
        double dy = poly_y[k_hi] - poly_y[k_lo];
        return poly_y[k_lo] + (t - poly_x[k_lo]) * dy / dx;
		break;
    }
    /*case Parametric : {
        break;
    }*/
    case FCT_Empty :
    case FCT_Linear :  // Linear doesn't exist yet and is then considered as identity
    default:
		return 0.5;
    }
}

void FlatCurve::getVal (const std::vector<double>& t, std::vector<double>& res) {

    res.resize (t.size());
    for (unsigned int i=0; i<t.size(); i++)
        res[i] = getVal(t[i]);
}

}

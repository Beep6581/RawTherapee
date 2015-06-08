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
#ifndef _COORD_H_
#define _COORD_H_

class PolarCoord;

// Do not confuse with rtengine::Coord2D, Coord is for the GUI
class Coord {
public:
	int x;
	int y;

	Coord() : x(-1), y(-1) {}
	Coord(int x, int y) : x(x), y(y) {}

	void set (int x, int y) {
		this->x = x;
		this->y = y;
	}

	void setFromPolar(PolarCoord polar);

	/// @brief Clip the coord to stay in the width x height bounds
	/// @return true if the x or y coordinate has changed
	bool clip(int width, int height) {
		int trimmedX = rtengine::LIM<int>(x, 0, width);
		int trimmedY = rtengine::LIM<int>(y, 0, height);
		bool retval = trimmedX!=x || trimmedY!=y;
		x = trimmedX;
		y = trimmedY;
		return retval;
	}

	void operator+=(const Coord & rhs) {
		x += rhs.x;
		y += rhs.y;
	}
	void operator-=(const Coord & rhs) {
		x -= rhs.x;
		y -= rhs.y;
	}
	void operator*=(double scale) {
		x *= scale;
		y *= scale;
	}
	Coord operator+(Coord & rhs) {
		Coord result(x+rhs.x, y+rhs.y);
		return result;
	}
	Coord operator-(Coord & rhs) {
		Coord result(x-rhs.x, y-rhs.y);
		return result;
	}
	Coord operator*(double scale) {
		Coord result(x*scale, y*scale);
		return result;
	}
};

class PolarCoord {
public:
	double radius;
	double angle; // degree

	PolarCoord() : radius(1.), angle(0.) {}
	PolarCoord(double radius, double angle) : radius(radius), angle(angle) {}

	void set (double radius, double angle) {
		this->radius = radius;
		this->angle = angle;
	}

	void setFromCartesian(Coord start, Coord end) {
		Coord delta(end.x-start.x, end.y-start.y);
		setFromCartesian(delta);
	}

	void setFromCartesian(Coord delta) {
		if (!delta.x && !delta.y) {
			// null vector, we set to a default value
			radius = 1.;
			angle = 0.;
			return;
		}
		double x_ = double(delta.x);
		double y_ = double(delta.y);
		radius = sqrt(x_*x_+y_*y_);
		if (delta.x>0.) {
			if (delta.y>=0.)
				angle = atan(y_/x_)/(2*M_PI)*360.;
			else if (delta.y<0.)
				angle = (atan(y_/x_)+2*M_PI)/(2*M_PI)*360.;
		}
		else if (delta.x<0.)
			angle = (atan(y_/x_)+M_PI)/(2*M_PI)*360.;
		else if (delta.x==0.) {
			if (delta.y>0.)
				angle = 90.;
			else
				angle = 270.;
		}
	}

	void operator+=(const PolarCoord & rhs) {
		Coord thisCoord, rhsCoord;
		thisCoord.setFromPolar(*this);
		rhsCoord.setFromPolar(rhs);
		thisCoord += rhsCoord;
		setFromCartesian(thisCoord);
	}
	void operator-=(const PolarCoord & rhs) {
		Coord thisCoord, rhsCoord;
		thisCoord.setFromPolar(*this);
		rhsCoord.setFromPolar(rhs);
		thisCoord -= rhsCoord;
		setFromCartesian(thisCoord);
	}
	void operator*=(double scale) {
		radius *= scale;
	}
	PolarCoord operator+(PolarCoord & rhs) {
		Coord thisCoord, rhsCoord;
		thisCoord.setFromPolar(*this);
		rhsCoord.setFromPolar(rhs);
		thisCoord += rhsCoord;
		PolarCoord result;
		result.setFromCartesian(thisCoord);
		return result;
	}
	PolarCoord operator-(PolarCoord & rhs) {
		Coord thisCoord, rhsCoord;
		thisCoord.setFromPolar(*this);
		rhsCoord.setFromPolar(rhs);
		thisCoord -= rhsCoord;
		PolarCoord result;
		result.setFromCartesian(thisCoord);
		return result;
	}
	Coord operator*(double scale) {
		Coord result(radius*scale, angle);
		return result;
	}

};

#endif

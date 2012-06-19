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
#ifndef _THRESHOLDSELECTOR_
#define _THRESHOLDSELECTOR_

#include "guiutils.h"
#include "../rtengine/procparams.h"

/*
 * This widget will let you select a linearly variable threshold, creating a ramp up
 * if you want to go from a null effect to a full effect
 *    0-0-ramp up-1-1
 * or a ramp down if you want the contrary
 *    1-1-ramp down-0-0
 *
 * You can optionally create a double threshold
 *    0-0-ramp up-1-1-ramp down-0-0
 * or
 *    1-1-ramp down-0-0-ramp up-1-1
 *
 * Please note that the values are related to the cursors, depending on their position
 * on the graph. E.g. the "bottomLeft" value is related to the bottom left cursor.
 */
class ThresholdSelector : public Gtk::DrawingArea {

	public:

		enum ThreshCursorId {
			TS_UNDEFINED=-1,
			TS_BOTTOMLEFT,
			TS_TOPLEFT,
			TS_BOTTOMRIGHT,
			TS_TOPRIGHT
		};


	protected:

		sigc::signal<void> sig_val_changed;

		Glib::RefPtr<Gdk::GC> gc_;
		Glib::RefPtr<Gdk::Pixmap> backBuffer;
		std::vector<GradientMilestone> bgGradient;

		bool doubleThresh;  // If true: there curve is a double threshold (0 to 1 to 0, or 1 to 0 to 1).
		bool initalEq1;     // If true: the curve start at 1 (top); if false: the curve start at 0 (bottom)
		unsigned int precision;  // Decimal number if this object has to handle "double" values
		ThreshCursorId litCursor;
		ThreshCursorId oldLitCursor;
		double boundary1[2], boundary2[2];
		double tmpX, tmpPos;

		ThreshCursorId movedCursor, secondaryMovedCursor;
		double minVal, maxVal;
		double defPos[4];
		double positions[4];
		unsigned short wslider;

		const static int hb = 3;  // horizontal border
		const static int vb = 2;  // vertical border

		void initValues (double minValue, double maxValue, bool startAtOne);
		void findLitCursor(int posX, int posY);
		void findSecondaryMovedCursor(guint state);
		void findBoundaries(double &min, double &max);
		double to01(double value);
		void updateTooltip();

	public:

		sigc::signal<void> signal_value_changed();

		ThresholdSelector(double minValue, double maxValue, double defBottom, double defTop, unsigned int precision, bool startAtOne);
		ThresholdSelector(double minValue, double maxValue, double defBottomLeft, double defTopLeft, double defBottomRight, double defTopRight, unsigned int precision, bool startAtOne);

		double shapeValue (double value) { return round(value*pow(double(10), precision)) / pow(double(10), precision); }

		template <typename T>
		void setDefaults (const rtengine::procparams::Threshold<T> &t) {
			defPos[TS_BOTTOMLEFT] = double(t.value[0]);  // should we use shapeValue() ?
			defPos[TS_TOPLEFT]    = double(t.value[1]);
			if (doubleThresh) {
				defPos[TS_BOTTOMRIGHT] = double(t.value[2]);
				defPos[TS_TOPRIGHT]    = double(t.value[3]);
			}
		}

		void setDefaults (double bottom, double top);
		void setDefaults (double bottomLeft, double topLeft, double bottomRight, double topRight);

		template <typename T>
		void setPositions (const rtengine::procparams::Threshold<T> &tValues) {
			positions[TS_BOTTOMLEFT]  = static_cast<double>(tValues.value[TS_BOTTOMLEFT]);
			positions[TS_TOPLEFT]     = static_cast<double>(tValues.value[TS_TOPLEFT]);
			if (tValues.isDouble()) {
				positions[TS_BOTTOMRIGHT] = static_cast<double>(tValues.value[TS_BOTTOMRIGHT]);
				positions[TS_TOPRIGHT]    = static_cast<double>(tValues.value[TS_TOPRIGHT]);
			}
			updateTooltip();
			queue_draw();
		}
		void setPositions (double bottom, double top);
		void setPositions (double bottomLeft, double topLeft, double bottomRight, double topRight);

		template <typename T>
		rtengine::procparams::Threshold<T> getPositions () {
			if (doubleThresh) {
				rtengine::procparams::Threshold<T> rThresh(
						static_cast<T>(shapeValue(positions[TS_BOTTOMLEFT])),
						static_cast<T>(shapeValue(positions[TS_TOPLEFT])),
						static_cast<T>(shapeValue(positions[TS_BOTTOMRIGHT])),
						static_cast<T>(shapeValue(positions[TS_TOPRIGHT])),
						initalEq1
				);
				return rThresh;
			}
			else {
				rtengine::procparams::Threshold<T> rThresh(
						static_cast<T>(shapeValue(positions[TS_BOTTOMLEFT])),
						static_cast<T>(shapeValue(positions[TS_TOPLEFT])),
						initalEq1
				);
				return rThresh;
			}
		}

		void getPositions (Glib::ustring& bottom, Glib::ustring& top);
		void getPositions (Glib::ustring& bottomLeft, Glib::ustring& topLeft, Glib::ustring& bottomRight, Glib::ustring& topRight);

		void setBgGradient (const std::vector<GradientMilestone> &milestones);
		bool isStartAtOne() { return initalEq1; }
		bool isDouble() { return doubleThresh; }
		void on_realize ();
		bool on_expose_event(GdkEventExpose* event);
		bool on_button_press_event (GdkEventButton* event);
		bool on_button_release_event (GdkEventButton* event);
		bool on_motion_notify_event (GdkEventMotion* event);
		bool on_leave_notify_event (GdkEventCrossing* event);
		void styleChanged (const Glib::RefPtr<Gtk::Style>& style);
		unsigned int getPrecision () { return precision; }
		void reset ();
};

#endif


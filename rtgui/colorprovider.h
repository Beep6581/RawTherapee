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
#ifndef _COLORPROVIDER_
#define _COLORPROVIDER_

#include <gtkmm.h>

class ColorProvider;

/*
 * The ColorCaller is the class that will query the ColorProvider
 */
class ColorCaller {
	protected:
		// a class can handle several ColorCaller;
		// colorCallerId will let the provider identify the caller
		int colorCallerId;
		ColorProvider* colorProvider;

	public:
		double ccRed;
		double ccGreen;
		double ccBlue;

		ColorCaller() : colorCallerId(-1), colorProvider(NULL), ccRed(0.), ccGreen(0.), ccBlue(0.) {}
		void setColorProvider (ColorProvider* p, int id) { colorProvider = p; colorCallerId = id; }
};

/*
 * Use it to let your widget feed a colored bar or graph lines with the wanted colors
 * If you doesn't need to dynamically feed a widget with colors (e.g. curve's graph),
 * you don't need to declare the instanciator class as BEING a ColorProvider, you'll
 * still be able to set gradients for e.g. ColoredBar(s)
 */
class ColorProvider {

	public:
		virtual ~ColorProvider() {};
		virtual void colorForValue (double valX, double valY, int callerId, ColorCaller* caller) {};
};

#endif

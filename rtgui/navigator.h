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
#ifndef _NAVIGATOR_
#define _NAVIGATOR_

#include <gtkmm.h>
#include <previewwindow.h>
#include <pointermotionlistener.h>

class Navigator : public Gtk::Frame, public PointerMotionListener {

	protected:
		Gtk::Label* position;
		Gtk::Label *R, *G, *B;
		Gtk::Label *H, *S, *V;

		void rgb2hsv 	(int r, int g, int b, int &h, int &s, int &v);
		void setInvalid ();
	public:
		PreviewWindow* previewWindow;

		Navigator ();

		// pointermotionlistener interface
		void pointerMoved (bool validPos, int x, int y, int r, int g, int b);

};

#endif

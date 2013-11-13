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

#ifndef _IMAGEDIMENSIONS_
#define _IMAGEDIMENSIONS_

class PreviewProps {
	public:
		int x, y, w, h, skip;
		PreviewProps (int _x, int _y, int _w, int _h, int _skip)
					 : x(_x), y(_y), w(_w), h(_h), skip(_skip) {}
};

/*
 * Description of an image dimension, with getter and setter
 */
class ImageDimensions {

	public:
		int width;
		int height;

	public:
		ImageDimensions() : width(-1), height(-1) {}
		int getW       () { return width;  }
		int getH       () { return height; }
		int getWidth   () { return width;  }
		int getHeight  () { return height; }
		void transform (PreviewProps pp, int tran, int &sx1, int &sy1, int &sx2, int &sy2);
};


#endif

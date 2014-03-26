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
#ifndef __CROPGUILISTENER__
#define __CROPGUILISTENER__

class CropGUIListener {

  public:
    virtual ~CropGUIListener() {}
    virtual void cropMoved          (int &x, int &y, int &w, int &h) =0;
    virtual void cropWidth1Resized  (int &x, int &y, int &w, int &h) =0;
    virtual void cropWidth2Resized  (int &x, int &y, int &w, int &h) =0;
    virtual void cropHeight1Resized (int &x, int &y, int &w, int &h) =0;
    virtual void cropHeight2Resized (int &x, int &y, int &w, int &h) =0;
    virtual void cropTopLeftResized     (int &x, int &y, int &w, int &h) =0;
    virtual void cropTopRightResized    (int &x, int &y, int &w, int &h) =0;
    virtual void cropBottomLeftResized  (int &x, int &y, int &w, int &h) =0;
    virtual void cropBottomRightResized (int &x, int &y, int &w, int &h) =0;
    virtual void cropInit           (int &x, int &y, int &w, int &h) =0;
    virtual void cropResized        (int &x, int &y, int& x2, int& y2) =0;
    virtual void cropManipReady     () =0;
    virtual bool inImageArea        (int x, int y) =0;
    virtual double getRatio         () =0;
};

#endif

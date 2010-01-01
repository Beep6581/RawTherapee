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
#include <image8.h>
#include <rtengine.h>

using namespace rtengine;

const unsigned char* IImage8::getData () {}

Image8::Image8 () 
  : width(-1), height(-1), data(NULL) {
}

Image8::Image8 (int w, int h) 
  : width(w), height (h), data(NULL) {
  
    allocate (w, h);
}

void Image8::allocate (int width, int height) {

    if (data!=NULL) 
        delete [] data;
  
    data = new unsigned char [width * height * 3];
    this->width = width;
    this->height = height;
}

Image8::~Image8 () {
  
    if (data!=NULL) 
        delete [] data;
}

void Image8::getScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;

    if (bps==8)
        memcpy (buffer, data + row*width*3, width*3);
    else if (bps==16) {
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0, ix = row*width*3; i<width*3; i++, ix++)
         sbuffer[i] =  data[ix] << 8;
    }
}

void Image8::setScanline (int row, unsigned char* buffer, int bps) {

    if (data==NULL)
        return;

    if (bps==8)
        memcpy (data + row*width*3, buffer, width*3);
    else if (bps==16) {
        unsigned short* sbuffer = (unsigned short*) buffer;
        for (int i=0, ix = row*width*3; i<width*3; i++, ix++)
          data[ix] = sbuffer[i] >> 8;
    }
}

unsigned char Image8::r (int row, int col) {

    return data[3*(row*width+col)];
}

unsigned char Image8::g (int row, int col) {

    return data[3*(row*width+col)+1];
}

unsigned char Image8::b (int row, int col) {

    return data[3*(row*width+col)+2];
}

void Image8::r (int row, int col, unsigned char val) {

    data[3*(row*width+col)] = val;
}

void Image8::g (int row, int col, unsigned char val) {

    data[3*(row*width+col)+1] = val;
}

void Image8::b (int row, int col, unsigned char val) {

    data[3*(row*width+col)+2] = val;
}


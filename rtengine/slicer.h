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
#ifndef _SLICER_
#define _SLICER_

//The image is divided in blocks even on single processor machine, mainly to decrease memory consumption
//maximum number of pixel per block
#define PIXELS_PER_BLOCK    250000

/*
 * Used to specify a subregion of an image and to specify a cell in this subregion
 */
class Block
{
public:
    unsigned int posX;
    unsigned int posY;
    unsigned int width;     // If 0, use the full width of the image
    unsigned int height;    // If 0, use the full height of the image
    Block();
    Block(unsigned int x, unsigned int y, unsigned int w, unsigned int h);
};

/*
 * This class handle the best slicing of the image with a given number of pixels per block and the number of
 * processor, and tries to create blocks as square as possible. There can be a different number of block on
 * each line, and the pixel per block requested may be oversized by very few percents.
 */
class Slicer
{
protected:
    bool portrait;                  // Orientation of the sub-region
    unsigned int imWidth;           // Image width
    unsigned int imHeight;          // Image height
    Block region;                   // Sub-region to process
    double hBlockNumber;            // Horizontal number of block for the sub-region
    unsigned int vBlockNumber;      // Vertical number of block for the sub-region
    double blockWidth;

public:
    unsigned int blockNumber;       // number of block for the sub-region
    unsigned int maxPixelNumber;    // number of pixel of the biggest block (for memory allocation purpose)
    Slicer(unsigned int imageWidth, unsigned int imageHeight, Block *subRegion, unsigned int pixels);
    void get_block(unsigned int blockId, Block *block);
};

#endif

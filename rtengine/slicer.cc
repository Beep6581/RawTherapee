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

#include <cmath>
#include <gtkmm.h>
#include "rt_math.h"

#include "slicer.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

// If no parameter set, everything = 0 -> process all the image
Block::Block() {
	posX = 0;
	posY = 0;
	width = 0;
	height = 0;
}

Block::Block(unsigned int x, unsigned int y, unsigned int w, unsigned int h) {
	posX = x;
	posY = y;
	width = w;
	height = h;
}

/*
 * Slice a sub-region to process in blocks who's size is given by the number of processor
 * and the number of pixel per block (and hence the memory footprint)
 */
Slicer::Slicer(unsigned int imageWidth, unsigned int imageHeight, Block *subRegion, unsigned int pixels ) {
	// If the sub-region has a portrait shape, X and Y coordinates are swapped for better result
	// It will be swapped back when sending back the block coordinates
	region.width = !(subRegion->width) ? imageWidth : subRegion->width;
	region.height = !(subRegion->height) ? imageHeight : subRegion->height;	// Assuming that the sub-region is under posY
	if (region.width < region.height) {
		region.width = !(subRegion->height) ? imageHeight : subRegion->height;
		region.height = !(subRegion->width) ? imageWidth : subRegion->width;	// Assuming that the sub-region is under posY
		portrait = true;
		imWidth = imageHeight;
		imHeight = imageWidth;
		region.posX = subRegion->posY;
		region.posY = subRegion->posX;
	}
	else {
		portrait = false;
		imWidth = imageWidth;
		imHeight = imageHeight;
		region.posX = subRegion->posX;
		region.posY = subRegion->posY;
	}
	double subRegionRatio = (double)(region.width) / (double)(region.height);

	//total number of core/processor
#ifdef _OPENMP
	unsigned int procNumber = omp_get_num_procs();
#else
	unsigned int procNumber = 1;
#endif

	//calculate the number of block
	blockNumber = (double(region.width*region.height) / (double)pixels);
	blockNumber = int((rtengine::max(blockNumber, 1U) + (double)procNumber/2.)/procNumber)*procNumber;
	vBlockNumber = (unsigned int)(sqrt((double)blockNumber / subRegionRatio)+0.5);
	vBlockNumber = CLAMP(vBlockNumber, 1, blockNumber);
	hBlockNumber = (double)blockNumber / (double)vBlockNumber;
	blockWidth = 1.0 / hBlockNumber;

	double maxPixelNumberX = (double)region.height / (double)vBlockNumber;
	double maxPixelNumberY = (double)region.width / (double)((unsigned int)hBlockNumber);
	if (maxPixelNumberX - (double)((unsigned int)maxPixelNumberX) != 0.) maxPixelNumberX += 1.;
	if (maxPixelNumberY - (double)((unsigned int)maxPixelNumberY) != 0.) maxPixelNumberY += 1.;
	maxPixelNumber = (unsigned int)maxPixelNumberX * (unsigned int)maxPixelNumberY;

}

// return the absolute position and size of the requested block
void Slicer::get_block(unsigned int numBlock, Block *block) {
	double roundingTradeOff = (hBlockNumber - (double)((int)hBlockNumber)) == 0.5 ? 2.1 : 2.0;
	unsigned int alreadyCompletedLineNbr = (unsigned int)((double)(numBlock) * blockWidth + (blockWidth/roundingTradeOff));

	unsigned int prevLineEnd = (unsigned int)((double)alreadyCompletedLineNbr * hBlockNumber + 0.5);
	unsigned int myLineEnd = (unsigned int)((double)(alreadyCompletedLineNbr+1) * hBlockNumber + 0.5);

	unsigned int nbrCellsOnMyLine = myLineEnd - prevLineEnd;
	unsigned int cellOnMyLine = numBlock - prevLineEnd;

	unsigned int blockStart = (unsigned int)(((double)region.width / (double)nbrCellsOnMyLine)*(double)(cellOnMyLine));
	unsigned int blockEnd = (unsigned int)(((double)region.width / (double)nbrCellsOnMyLine)*(double)(cellOnMyLine+1));
	block->width = blockEnd - blockStart;
	block->posX = region.posX + blockStart;
	if (cellOnMyLine == (nbrCellsOnMyLine-1)) {
		// We make sure that the last block of the row take the rest of the remaining X space
		block->width = region.posX + region.width - block->posX;
	}

	blockStart = (unsigned int)(((double)region.height / (double)vBlockNumber)*(double)(alreadyCompletedLineNbr));
	blockEnd = (unsigned int)(((double)region.height / (double)vBlockNumber)*(double)(alreadyCompletedLineNbr+1));
	block->height = blockEnd - blockStart;
	block->posY = region.posY + blockStart;
	if (alreadyCompletedLineNbr == (vBlockNumber-1)) {
		block->height = region.posY + region.height - block->posY;
	}

	if (portrait) {
		// we swap back the X/Y coordinates
		unsigned int temp;

		temp = block->posX;
		block->posX = block->posY;
		block->posY = temp;

		temp = block->width;
		block->width = block->height;
		block->height = temp;

	}
}

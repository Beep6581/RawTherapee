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

#include <math.h>
#include <omp.h>
#include <slicer.h>

using namespace rtengine;

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
Slicer::Slicer(unsigned int imageWidth, unsigned int imageHeight, Block *subRegion, unsigned int pixels, const char* nomFichier) {
	maxPixelNumberR = 0;
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
	procNumber = omp_get_num_procs();
#else
	procNumber = 1;
#endif

	//calculate the number of block
	blockNumber = (double(region.width*region.height) / (double)pixels);
	blockNumber = int((MAX(blockNumber, 1) + (double)procNumber/2.)/procNumber)*procNumber;
	vBlockNumber = (unsigned int)(sqrt((double)blockNumber / subRegionRatio)+0.5);
	vBlockNumber = CLAMP(vBlockNumber, 1, blockNumber);
	hBlockNumber = (double)blockNumber / (double)vBlockNumber;
	blockWidth = 1.0 / hBlockNumber;
	blockHeight = 1.0 / (double)vBlockNumber;

	double maxPixelNumberX = (double)region.height / (double)vBlockNumber;
	double maxPixelNumberY = (double)region.width / (double)((unsigned int)hBlockNumber);
	if (maxPixelNumberX - (double)((unsigned int)maxPixelNumberX) != 0) maxPixelNumberX += 1.;
	if (maxPixelNumberY - (double)((unsigned int)maxPixelNumberY) != 0) maxPixelNumberY += 1.;
	maxPixelNumber = (unsigned int)maxPixelNumberX * (unsigned int)maxPixelNumberY;

	printf ("\n*****************\nimageWidth : %d\nimageHeight : %d\nregion.width : %d\nregion.height : %d\npixels : %d\n\n", imageWidth, imageHeight, region.width, region.height, pixels);
	printf ("subRegionRatio=%.5f\nPortrait = %d\n", (float)subRegionRatio, portrait);
	printf ("Total block number : %d\nHorizontal block number : %.3f\nVertical block number : %d\nmaxPixelNumber  : %d\n", blockNumber, hBlockNumber, vBlockNumber, maxPixelNumber);
	fichierSortie = fopen(nomFichier, "wt");
	fprintf(fichierSortie,"<?xml version='1.0' encoding='utf-8' standalone='no'?>\n");
	fprintf(fichierSortie,"<svg xmlns:svg='http://www.w3.org/2000/svg' xmlns='http://www.w3.org/2000/svg' version='1.0' ");
	fprintf(fichierSortie,"width='%d' height='%d' id='svg2'><defs id='defs4' />\n", (portrait?imageHeight:imageWidth), (portrait?imageWidth:imageHeight));

	fprintf(fichierSortie,"<rect width='%d' height='%d' x='%d' y='%d'", ((portrait?region.height:region.width)-2), ((portrait?region.width:region.height)-2), ((portrait?region.posY:region.posX)+1), ((portrait?region.posX:region.posY)+1));
	fprintf(fichierSortie," id='rectSR' style='opacity:0.3;fill:none;");
	fprintf(fichierSortie,"fill-opacity:1;fill-rule:evenodd;stroke:#FF0000;stroke-width:2;stroke-linecap:square;stroke-linejoin:miter;");
	fprintf(fichierSortie,"stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:1' />\n");

}

Slicer::~Slicer() {
	printf("maxPixelNumberR = %d (block #%d)\n", maxPixelNumberR, maxPixBlockNumber);
	fprintf(fichierSortie,"</svg>\n");
	fclose(fichierSortie);
}

// return the absolute position and size of the requested block
void Slicer::get_block(unsigned int numBlock, Block *block) {
	volatile double roundingTradeOff = (hBlockNumber - (double)((int)hBlockNumber)) == 0.5 ? 2.1 : 2.0;
	volatile unsigned int nbrLigneDejaComplete = (unsigned int)((double)(numBlock) * blockWidth + (blockWidth/roundingTradeOff));

	volatile unsigned int prevLigneFin = (unsigned int)((double)nbrLigneDejaComplete * hBlockNumber + 0.5);
	volatile unsigned int maLigneFin = (unsigned int)((double)(nbrLigneDejaComplete+1) * hBlockNumber + 0.5);

	volatile unsigned int nbCelluleSurMaLigne = maLigneFin - prevLigneFin;
	volatile unsigned int celluleSurMaLigne = numBlock - prevLigneFin;

	volatile unsigned int blockStart = (unsigned int)(((double)region.width / (double)nbCelluleSurMaLigne)*(double)(celluleSurMaLigne));
	volatile unsigned int blockEnd = (unsigned int)(((double)region.width / (double)nbCelluleSurMaLigne)*(double)(celluleSurMaLigne+1));
	block->width = blockEnd - blockStart;
	block->posX = region.posX + blockStart;
	if (celluleSurMaLigne == (nbCelluleSurMaLigne-1)) {
		// We make sure that the last block of the row take the rest of the remaining X space
		block->width = region.posX + region.width - block->posX;
	}

	blockStart = (unsigned int)(((double)region.height / (double)vBlockNumber)*(double)(nbrLigneDejaComplete));
	blockEnd = (unsigned int)(((double)region.height / (double)vBlockNumber)*(double)(nbrLigneDejaComplete+1));
	block->height = blockEnd - blockStart;
	block->posY = region.posY + blockStart;
	if (nbrLigneDejaComplete == (vBlockNumber-1)) {
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
	unsigned int currBlockSize = block->width * block->height;
	if (currBlockSize > maxPixelNumberR) {
		maxPixelNumberR = currBlockSize;
		maxPixBlockNumber = numBlock+1;
	}

	char *couleur[4] = {"70d26e", "3ba6c3", "d0d26e", "3b62c3"};
	char *coul = couleur[(nbrLigneDejaComplete%2)*2 + celluleSurMaLigne%2];
	fprintf(fichierSortie,"<rect width='%d' height='%d' x='%d' y='%d'", (block->width-2), (block->height-2), (block->posX+1), (block->posY+1));
	fprintf(fichierSortie," id='rect%d' style='opacity:0.3;fill:#%s;fill-opacity:1;fill-rule:evenodd;stroke:#%s;stroke-width:2;stroke-linecap:square;", numBlock, coul, coul);
	fprintf(fichierSortie,"stroke-linejoin:miter;stroke-miterlimit:4;stroke-dasharray:none;stroke-dashoffset:0;stroke-opacity:1' />\n", numBlock);
	fprintf(fichierSortie,"<text x='%d' y='%d'", (block->posX + block->width/2), (block->posY + block->height/2));
	fprintf(fichierSortie," id='text%d' xml:space='preserve' style='font-size:15px;font-style:normal;", numBlock);
	fprintf(fichierSortie,"font-variant:normal;font-weight:normal;font-stretch:normal;text-align:start;text-anchor:start;fill:#000000;");
	fprintf(fichierSortie,"fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;");
	fprintf(fichierSortie,"font-family:Arial;-inkscape-font-specification:Arial'>\n");
	fprintf(fichierSortie,"<tspan x='%d' y='%d'", (block->posX + block->width/2), (block->posY + block->height/2));
	fprintf(fichierSortie," id='tspan%d'>%d</tspan>\n</text>\n", numBlock, (numBlock+1));
}

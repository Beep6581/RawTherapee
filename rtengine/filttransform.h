/*
 * filttransform.h
 *
 *  Created on: Oct 4, 2010
 *      Author: gabor
 */

#ifndef FILTTRANSFORM_H_
#define FILTTRANSFORM_H_

#include "filter.h"
#include <vector>

//
// T r a n s f o r m   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

namespace rtengine {

class TransformFilterDescriptor : public FilterDescriptor {

	public:
        TransformFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern TransformFilterDescriptor transformFilterDescriptor;

class Coord2D {

    public:
        double x, y;
        Coord2D (double x_, double y_) : x(x_), y(y_) {}
        Coord2D () {}
        void set (double x_, double y_) { x = x_; y = y_; }
};

class TransformFilter : public Filter {

        bool needsCA            ();
        bool needsDistortion    ();
        bool needsRotation      ();
        bool needsPerspective   ();
        bool needsVignetting    ();
        bool needsTransform     ();
        void simpltransform     (MultiImage* sourceImage, MultiImage* targetImage);
        void vignetting         (MultiImage* sourceImage, MultiImage* targetImage);
        void transformNonSep    (MultiImage* sourceImage, MultiImage* targetImage);
        void transformSep       (MultiImage* sourceImage, MultiImage* targetImage);
        inline void cubintch    (unsigned short** src, int xs, int ys, double Dx, double Dy, unsigned short *r, double mul);
        inline void cubint      (MultiImage* src, int xs, int ys, double Dx, double Dy, unsigned short *r, unsigned short *g, unsigned short *b, double mul);
        bool transCoord         (Dim fullSize, std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1);
        bool transCoord         (Dim fullSize, ImageView target, ImageView& source, double ascaleDef = -1);
        double getTransformAutoFill ();

    public:
        TransformFilter ();

    	void      process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
        ImageView calculateSourceImageView (const ImageView& requestedImView);
        void      reverseTransPoint (int x, int y, int& xv, int& yv);
};

}
#endif /* FILTTRANSFORM_H_ */

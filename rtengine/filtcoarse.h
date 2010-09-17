/*
 * filtcoarse.h
 *
 *  Created on: Sep 17, 2010
 *      Author: gabor
 */

#ifndef FILTCOARSE_H_
#define FILTCOARSE_H_

#include "filter.h"

//
// C o a r s e T r a n s f o r m   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Must be applied after demosaicing!

namespace rtengine {

class CoarseTransformFilterDescriptor : public FilterDescriptor {

	public:
        CoarseTransformFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern CoarseTransformFilterDescriptor coarseTransformFilterDescriptor;

class CoarseTransformFilter : public Filter {

        void vflip (MultiImage* image);
        void hflip (MultiImage* image);
        void rotate90  (unsigned short** si, unsigned short* ti, int sW, int sH, Buffer<int>* buffer);
        void rotate180 (unsigned short** si, unsigned short* ti, int sW, int sH, Buffer<int>* buffer);
        void rotate270 (unsigned short** si, unsigned short* ti, int sW, int sH, Buffer<int>* buffer);

    public:
        CoarseTransformFilter ();

    	void      process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
        ImageView calculateTargetImageView (const ImageView& requestedImView);
        ImageView calculateSourceImageView (const ImageView& requestedImView);
        void      getFullImageSize (int& w, int& h);
        void      getReqiredBufferSize (int& w, int& h);
        void      reverseTransPoint (int x, int y, int& xv, int& yv);
};

}
#endif /* FILTCOARSE_H_ */

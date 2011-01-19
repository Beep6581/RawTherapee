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
		void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern CoarseTransformFilterDescriptor coarseTransformFilterDescriptor;

class CoarseTransformFilter : public Filter {

        void vflip (MultiImage* image);
        void hflip (MultiImage* image);
        void rotate90  (float** si, float** ti, int sW, int sH, Buffer<float>* buffer);
        void rotate180 (float** si, float** ti, int sW, int sH, Buffer<float>* buffer);
        void rotate270 (float** si, float** ti, int sW, int sH, Buffer<float>* buffer);

    public:
        CoarseTransformFilter ();

    	void      process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
        ImageView calculateTargetImageView (const ImageView& requestedImView);
        ImageView calculateSourceImageView (const ImageView& requestedImView);
        Dim       getFullImageSize ();
        Dim       getReqiredBufferSize ();
        void      reverseTransPoint (int x, int y, int& xv, int& yv);
};

}
#endif /* FILTCOARSE_H_ */

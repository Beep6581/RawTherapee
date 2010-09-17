/*
 * filtsharpener.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTSHARPENER_H_
#define FILTSHARPENER_H_

#include "filter.h"
#include "multiimage.h"

//
// S h a r p e n   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//

namespace rtengine {

class SharpenFilterDescriptor : public FilterDescriptor {

	public:
    SharpenFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern SharpenFilterDescriptor sharpenFilterDescriptor;

class SharpenFilter : public Filter {

        void dcdamping (Buffer<float>* aI, MultiImage* aO, float damping);
        void deconvsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* b2);
        void sharpenHaloCtrl (MultiImage* sourceImage, MultiImage* targetImage, Buffer<unsigned short>* blurmap, Buffer<unsigned short>* base);
        void usmsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<unsigned short>* b2);

	public:
        SharpenFilter ();
        void getReqiredBufferSize (int& w, int& h);
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTSHARPENER_H_ */

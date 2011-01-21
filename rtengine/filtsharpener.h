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
    	void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern SharpenFilterDescriptor sharpenFilterDescriptor;

class SharpenFilter : public Filter {

        void dcdamping (Buffer<float>* aI, MultiImage* aO, float damping);
        void deconvsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* b2);
        void sharpenHaloCtrl (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* blurmap, Buffer<float>* base);
        void usmsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* b2);

	public:
        SharpenFilter ();
        Dim  getReqiredBufferSize ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTSHARPENER_H_ */

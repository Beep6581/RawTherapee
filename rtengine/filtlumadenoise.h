/*
 * filtlumadenoise.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTLUMADENOISE_H_
#define FILTLUMADENOISE_H_

#include "filter.h"
#include "multiimage.h"

//
// L u m i n a n c e   D e n o i s i n g   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//

namespace rtengine {

class LumaDenoiseFilterDescriptor : public FilterDescriptor {

	public:
    LumaDenoiseFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern LumaDenoiseFilterDescriptor lumaDenoiseFilterDescriptor;

class LumaDenoiseFilter : public Filter {

	public:
        LumaDenoiseFilter ();
        Dim  getReqiredBufferSize ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTLUMADENOISE_H_ */

/*
 * filtcsconv.h
 *
 *  Created on: Sep 23, 2010
 *      Author: gabor
 */

#ifndef FILTCSCONV_H_
#define FILTCSCONV_H_

#include "filter.h"
#include "multiimage.h"
#include <lcms2.h>

//
// C o l o r   S p a c e   C o n v e r s i o n   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Assumptions: it must be applied after demosaicing, after highlight recovery
//

namespace rtengine {

class ColorSpaceConvFilterDescriptor : public FilterDescriptor {

	public:
        ColorSpaceConvFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern ColorSpaceConvFilterDescriptor colorSpaceConvFilterDescriptor;

class ColorSpaceConvFilter : public Filter {

        cmsHTRANSFORM hTransform;
        cmsHPROFILE trIn;
        cmsHPROFILE trOut;

	public:
        ColorSpaceConvFilter ();
        ~ColorSpaceConvFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTCSCONV_H_ */

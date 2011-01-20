/*
 * filtlumacurve.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTLUMACURVE_H_
#define FILTLUMACURVE_H_

#include "filter.h"
#include "multiimage.h"

//
// L u m i n a n c e   C u r v e    f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

namespace rtengine {

class LumaCurveFilterDescriptor : public FilterDescriptor {

	public:
    	LumaCurveFilterDescriptor ();
    	void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern LumaCurveFilterDescriptor lumaCurveFilterDescriptor;

class LumaCurveFilter : public Filter {

        float* curve;
        unsigned int* histogram;

	public:
        LumaCurveFilter ();
        ~LumaCurveFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTLUMACURVE_H_ */

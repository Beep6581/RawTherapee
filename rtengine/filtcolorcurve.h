/*
 * filtcolorcurve.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTCOLORCURVE_H_
#define FILTCOLORCURVE_H_

#include "filter.h"
#include "multiimage.h"

//
// C o l o r   C u r v e   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//

namespace rtengine {

class ColorCurveFilterDescriptor : public FilterDescriptor {

	public:
        ColorCurveFilterDescriptor ();
		void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern ColorCurveFilterDescriptor colorCurveFilterDescriptor;

class ColorCurveFilter : public Filter {

        float* curve;

        void generateCurve (float boost, float limit);

	public:
        ColorCurveFilter ();
        ~ColorCurveFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTCOLORCURVE_H_ */

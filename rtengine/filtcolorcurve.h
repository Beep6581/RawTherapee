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
		void createAndAddToList (Filter* tail) const;
};

extern ColorCurveFilterDescriptor colorCurveFilterDescriptor;

class ColorCurveFilter : public Filter {

        double* curve;

        void generateCurve ();

	public:
        ColorCurveFilter ();
        ~ColorCurveFilter ();
    	void process (MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTCOLORCURVE_H_ */

/*
 * filthlrec.h
 *
 *  Created on: Sep 23, 2010
 *      Author: gabor
 */

#ifndef FILTHLREC_H_
#define FILTHLREC_H_

#include "filter.h"
#include "multiimage.h"
#include "matrix33.h"

//
// H i g h l i g h t   R e c o v e r y   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Assumptions: it must be applied after demosaicing, before color space conversion
//

namespace rtengine {

class HighlightRecoveryFilterDescriptor : public FilterDescriptor {

	public:
        HighlightRecoveryFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern HighlightRecoveryFilterDescriptor highlightRecoveryFilterDescriptor;

class HighlightRecoveryFilter : public Filter {

        void luminance (MultiImage* sourceImage, MultiImage* targetImage, int maxval);
        void cieblend (MultiImage* sourceImage, MultiImage* targetImage, int maxval, Matrix33 cam);

	public:
        HighlightRecoveryFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTHLREC_H_ */

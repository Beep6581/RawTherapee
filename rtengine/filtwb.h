/*
 * filtwb.h
 *
 *  Created on: Aug 27, 2010
 *      Author: gabor
 */

#ifndef FILTWB_H_
#define FILTWB_H_

#include "filter.h"
#include "multiimage.h"

//
// W h i t e   B a l a n c e   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Assumptions: it must be applied before color space conversion
//

namespace rtengine {

class WhiteBalanceFilterDescriptor : public FilterDescriptor {

	public:
		WhiteBalanceFilterDescriptor ();
    	void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern WhiteBalanceFilterDescriptor whiteBalanceFilterDescriptor;

class WhiteBalanceFilter : public Filter {

	public:
		WhiteBalanceFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTWB_H_ */

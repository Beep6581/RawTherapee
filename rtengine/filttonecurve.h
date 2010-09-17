/*
 * filttonecurve.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTTONECURVE_H_
#define FILTTONECURVE_H_

#include "filter.h"
#include "multiimage.h"

//
// T o n e   C u r v e    f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// This is a two-phase filter. The pre-filter is responsible to handle auto-exposure calculations
// and maintains the histogram that is used when calculating the tone curve.
// This way the histogram and the auto exp values need not to be recalculated when curve parameters
// are changed.

namespace rtengine {

class ToneCurveFilterDescriptor : public FilterDescriptor {

	public:
        ToneCurveFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern ToneCurveFilterDescriptor toneCurveFilterDescriptor;

class PreToneCurveFilter : public Filter {

        int*   histogram;
    public:
        PreToneCurveFilter ();
        ~PreToneCurveFilter ();
        void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
        int* getHistogram ();
};

class ToneCurveFilter : public Filter {

        int* curve;
        int* bchistogram;
        PreToneCurveFilter* ptcFilter;

	public:
        ToneCurveFilter ();
        ~ToneCurveFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTTONECURVE_H_ */

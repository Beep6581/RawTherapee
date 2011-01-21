/*
 * filtsh.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTSH_H_
#define FILTSH_H_

#include "filter.h"

//
// S h a d o w s / H i g h l i g h t s   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// This is a two-phase filter. The pre-filter is responsible to calculate shadow map
// and its statistics (max, min, avg).
// This way the shadow map does not need not to be recalculated when shadows/highlights parameters
// are changed.

namespace rtengine {

class ShadowsHighlightsFilterDescriptor : public FilterDescriptor {

	public:
        ShadowsHighlightsFilterDescriptor ();
		void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern ShadowsHighlightsFilterDescriptor shadowsHighlightsFilterDescriptor;

class PreShadowsHighlightsFilter : public Filter {

        Buffer<float>* map;
        float   max, min, avg;

    public:
        PreShadowsHighlightsFilter ();
        ~PreShadowsHighlightsFilter ();
        void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
        float** getSHMap ();
        float   getMapMax ();
        float   getMapMin ();
        float   getMapAvg ();
        Dim getReqiredBufferSize ();
};


class ShadowsHighlightsFilter : public Filter {

        PreShadowsHighlightsFilter* pshFilter;

	public:
        ShadowsHighlightsFilter (PreShadowsHighlightsFilter* pshf);
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTSH_H_ */

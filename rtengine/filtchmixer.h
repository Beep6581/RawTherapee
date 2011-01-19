/*
 * filtchmixer.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTCHMIXER_H_
#define FILTCHMIXER_H_

#include "filter.h"
#include "multiimage.h"

//
// C o l o r   M i x e r   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//

namespace rtengine {

class ColorMixerFilterDescriptor : public FilterDescriptor {

	public:
        ColorMixerFilterDescriptor ();
		void getDefaultParameters (ProcParams& defProcParams) const;
		void createAndAddToList (Filter* tail) const;
};

extern ColorMixerFilterDescriptor colorMixerFilterDescriptor;

class ColorMixerFilter : public Filter {

	public:
        ColorMixerFilter ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer);
};

}
#endif /* FILTCHMIXER_H_ */

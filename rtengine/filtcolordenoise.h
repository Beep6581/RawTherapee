/*
 * filtcolordenoise.h
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#ifndef FILTCOLORDENOISE_H_
#define FILTCOLORDENOISE_H_

#include "filter.h"
#include "multiimage.h"

//
// C o l o r   D e n o i s i n g   f i l t e r
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//

namespace rtengine {

class ColorDenoiseFilterDescriptor : public FilterDescriptor {

	public:
        ColorDenoiseFilterDescriptor ();
		void createAndAddToList (Filter* tail) const;
};

extern ColorDenoiseFilterDescriptor colorDenoiseFilterDescriptor;

class ColorDenoiseFilter : public Filter {

	public:
        ColorDenoiseFilter ();
        Dim  getReqiredBufferSize ();
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);
};

}
#endif /* FILTCOLORDENOISE_H_ */

/*
 * filterfactory.cc
 *
 *  Created on: Aug 18, 2010
 *      Author: gabor
 */

#include "filterfactory.h"

namespace rtengine {

extern FilterFactory filterFactory;

FilterFactory::FilterFactory() {

    registerFilterDescriptor (&whiteBalanceFilterDescriptor);
    registerFilterDescriptor (&demosaicFilterDescriptor);
    registerFilterDescriptor (&toneCurveFilterDescriptor);
    registerFilterDescriptor (&lumaDenoiseFilterDescriptor);
    registerFilterDescriptor (&colorDenoiseFilterDescriptor);
    registerFilterDescriptor (&colorMixerFilterDescriptor);
    registerFilterDescriptor (&shadowsHighlightsFilterDescriptor);
    registerFilterDescriptor (&sharpenFilterDescriptor);
    registerFilterDescriptor (&colorCurveFilterDescriptor);
}

void FilterFactory::registerFilterDescriptor (FilterDescriptor* descr) {

	filterDescriptors[descr->getName ()] = descr;
}

void FilterFactory::createFilterAddToList (const std::string& name, Filter* tail) {

	FilterDescriptor* fDescr = filterDescriptors[name];
	if (fDescr)
	    fDescr->createAndAddToList (tail);
}
}

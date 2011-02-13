/*
 * filterfactory.cc
 *
 *  Created on: Aug 18, 2010
 *      Author: gabor
 */
#include "procparams.h"
#include "filterfactory.h"
#include "filtwb.h"
#include "filtdemosaic.h"
#include "filttonecurve.h"
#include "filtlumadenoise.h"
#include "filtlumacurve.h"
#include "filtcolordenoise.h"
#include "filtchmixer.h"
#include "filtcoarse.h"
#include "filtsh.h"
#include "filtsharpener.h"
#include "filtcolorcurve.h"
#include "filthlrec.h"
#include "filtcsconv.h"
#include "filtresize.h"
#include "filttransform.h"

namespace rtengine {

FilterFactory* filterFactory;

FilterFactory::FilterFactory() {

    registerFilterDescriptor (&whiteBalanceFilterDescriptor);
    registerFilterDescriptor (&demosaicFilterDescriptor);
    registerFilterDescriptor (&highlightRecoveryFilterDescriptor);
    registerFilterDescriptor (&colorSpaceConvFilterDescriptor);
    registerFilterDescriptor (&toneCurveFilterDescriptor);
    registerFilterDescriptor (&coarseTransformFilterDescriptor);
    registerFilterDescriptor (&lumaDenoiseFilterDescriptor);
    registerFilterDescriptor (&lumaCurveFilterDescriptor);
    registerFilterDescriptor (&colorDenoiseFilterDescriptor);
    registerFilterDescriptor (&colorMixerFilterDescriptor);
    registerFilterDescriptor (&shadowsHighlightsFilterDescriptor);
    registerFilterDescriptor (&sharpenFilterDescriptor);
    registerFilterDescriptor (&colorCurveFilterDescriptor);
    registerFilterDescriptor (&resizeFilterDescriptor);
    registerFilterDescriptor (&transformFilterDescriptor);
}

void FilterFactory::registerFilterDescriptor (FilterDescriptor* descr) {

	filterDescriptors[descr->getName ()] = descr;
	descr->getDefaultParameters (defaultProcParams);
}

FilterDescriptor* FilterFactory::getFilterDescriptor(const String& name) {

    return filterDescriptors[name];
}
}

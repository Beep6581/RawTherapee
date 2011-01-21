/*
 * filtsh.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtsh.h"
#include "rtengine.h"
#include "macros.h"
#include "improcfuns.h"
#include "buffer.h"
#include "gauss.h"
#include "bilateral2.h"
#include "filterchain.h"
#include "iccstore.h"

namespace rtengine {

class PreShadowsHighlightsFilterDescriptor : public FilterDescriptor {
    public:
		PreShadowsHighlightsFilterDescriptor ()
				: FilterDescriptor ("PreShadowsHighlights", MultiImage::RGB, MultiImage::RGB) {
				addTriggerEvent (EvSHEnabled);
				addTriggerEvent (EvSHRadius);
			}
        void createAndAddToList (Filter* tail) const {}
};

PreShadowsHighlightsFilterDescriptor preShadowsHighlightsFilterDescriptor;
ShadowsHighlightsFilterDescriptor shadowsHighlightsFilterDescriptor;

ShadowsHighlightsFilterDescriptor::ShadowsHighlightsFilterDescriptor ()
	: FilterDescriptor ("ShadowsHighlights", MultiImage::RGB, MultiImage::RGB) {

    addTriggerEvent (EvSHEnabled);
	addTriggerEvent (EvSHHighlights);
    addTriggerEvent (EvSHShadows);
    addTriggerEvent (EvSHHLTonalW);
    addTriggerEvent (EvSHSHTonalW);
    addTriggerEvent (EvSHLContrast);
}

void ShadowsHighlightsFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setBoolean ("ShadowsHighlightsEnabled", false);
	defProcParams.setBoolean ("ShadowsHighlightsHQ", false);
	defProcParams.setFloat   ("ShadowsHighlightsHighlights", 10);
	defProcParams.setFloat   ("ShadowsHighlightsShadows", 10);
	defProcParams.setFloat   ("ShadowsHighlightsHTonalWidth", 80);
	defProcParams.setFloat   ("ShadowsHighlightsSTonalWIdth", 80);
	defProcParams.setFloat   ("ShadowsHighlightsLocalContrast", 0);
	defProcParams.setFloat   ("ShadowsHighlightsRadius", 40);
}

void ShadowsHighlightsFilterDescriptor::createAndAddToList (Filter* tail) const {

    PreShadowsHighlightsFilter* pshf = new PreShadowsHighlightsFilter ();
    tail->addNext (pshf);
    pshf->addNext (new ShadowsHighlightsFilter (pshf));
}

PreShadowsHighlightsFilter::PreShadowsHighlightsFilter ()
    : Filter (&preShadowsHighlightsFilterDescriptor), map (NULL) {
}

PreShadowsHighlightsFilter::~PreShadowsHighlightsFilter () {

    delete map;
}

void PreShadowsHighlightsFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	bool enabled = procParams->getBoolean ("ShadowsHighlightsEnabled");

    if (enabled) {

    	float shradius 	= procParams->getFloat ("ShadowsHighlightsRadius");
    	bool shhq 		= procParams->getFloat ("ShadowsHighlightsHQ");

        // calculate radius (procparams contains relative radius only)
        ImageSource* imgsrc = getFilterChain ()->getImageSource ();
        Dim fsize = imgsrc->getFullImageSize ();
        float radius = sqrt ((float)(fsize.width*fsize.width+fsize.height*fsize.height)) / 2.0 / 1800.0 * shradius;

        // allocate map
        delete map;
        map = new Buffer<float> (sourceImage->width, sourceImage->height);

        // fill with luminance
        Matrix33 wprof = iccStore.workingSpaceMatrix (procParams->icm.working);
        float lumimulr = wprof.data[1][0];
        float lumimulg = wprof.data[1][1];
        float lumimulb = wprof.data[1][2];

        for (int i=0; i<map->height; i++)
            for (int j=0; j<map->width; j++)
            	map->rows[i][j] = lumimulr*sourceImage->r[i][j] + lumimulg*sourceImage->g[i][j] + lumimulb*sourceImage->b[i][j];

		Dim size (sourceImage->width, sourceImage->height);
		gaussHorizontal<float> (map, map, size, (double*)(buffer->data), radius, multiThread);
		gaussVertical<float>   (map, map, size, (double*)(buffer->data), radius, multiThread);

        // update average, minimum, maximum
        Filter* p = getParentFilter ();
        if (!p) {
            avg = 0;
            int n = 1;
            min = FLT_MAX;
            max = 0;
            for (int i=32; i<map->height-32; i++)
                for (int j=32; j<map->width-32; j++) {
                    float val = map->rows[i][j];
                    if (val < min)
                        min = val;
                    if (val > max)
                        max = val;
                    avg = 1.0/n * val + (1.0 - 1.0/n) * avg;
                    n++;
                }
        }
        else {
            Filter* root = p;
            while (root->getParentFilter())
                root = root->getParentFilter();
            min = ((PreShadowsHighlightsFilter*)root)->min;
            max = ((PreShadowsHighlightsFilter*)root)->max;
            avg = ((PreShadowsHighlightsFilter*)root)->avg;
        }
    }
    else {
        delete map;
        map = NULL;
    }

    // we have to copy image data if input and output are not the same
    if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

float** PreShadowsHighlightsFilter::getSHMap () {

    return map->rows;
}

float PreShadowsHighlightsFilter::getMapMax () {

    return max;
}

float PreShadowsHighlightsFilter::getMapMin () {

    return min;
}

float PreShadowsHighlightsFilter::getMapAvg () {

    return avg;
}

Dim PreShadowsHighlightsFilter::getReqiredBufferSize () {

    if (procParams->getBoolean ("ShadowsHighlightsEnabled")) {
        Dim sdim = getScaledTargetImageView().getSize();
    	bool shhq = procParams->getFloat ("ShadowsHighlightsHQ");
        if (!shhq) {
            if (sdim.height > sdim.width)
                return Dim (2, sdim.height*omp_get_max_threads());
            else
                return Dim (sdim.width*omp_get_max_threads(), 2);
        }
        else
            return sdim;
    }
    else
        return Dim ();
}

ShadowsHighlightsFilter::ShadowsHighlightsFilter (PreShadowsHighlightsFilter* pshf)
	: Filter (&shadowsHighlightsFilterDescriptor), pshFilter (pshf) {
}

void ShadowsHighlightsFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

    // apply filter
	bool enabled 	 = procParams->getBoolean ("ShadowsHighlightsEnabled");
	float highlights = procParams->getFloat ("ShadowsHighlightsHighlights");
	float shadows    = procParams->getFloat ("ShadowsHighlightsShadows");
	float lce        = procParams->getFloat ("ShadowsHighlightsLocalContrast");

	bool processSH  = enabled && (highlights>0 || shadows>0);
    bool processLCE = enabled && lce>0;
    double lceamount = lce / 200.0;

    if (processSH || processLCE) {

    	float htonalwidth = procParams->getFloat ("ShadowsHighlightsHTonalWidth");
    	float stonalwidth = procParams->getFloat ("ShadowsHighlightsSTonalWidth");

        float** shmap = pshFilter->getSHMap ();

        float h_th, s_th;
        if (shmap) {
            h_th = pshFilter->getMapMax() - htonalwidth * (pshFilter->getMapMax() - pshFilter->getMapAvg()) / 100.0;
            s_th = stonalwidth * (pshFilter->getMapAvg() - pshFilter->getMapMin()) / 100.0;
        }

        Matrix33 wprof = iccStore.workingSpaceMatrix (procParams->icm.working);
        float lumimulr = wprof.rowsum(0);
        float lumimulg = wprof.rowsum(1);
        float lumimulb = wprof.rowsum(2);

        #pragma omp parallel for if (multiThread)
        for (int i=0; i<sourceImage->height; i++) {
            for (int j=0; j<sourceImage->width; j++) {
                float r = sourceImage->r[i][j];
                float g = sourceImage->g[i][j];
                float b = sourceImage->b[i][j];
                float mapval = shmap[i][j];
                double factor = 1.0;
                if (processSH) {
                    if (mapval > h_th)
                        factor = (h_th + (100.0 - highlights) * (mapval - h_th) / 100.0) / mapval;
                    else if (mapval < s_th)
                        factor = (s_th - (100.0 - shadows) * (s_th - mapval) / 100.0) / mapval;
                }
                if (processLCE) {
                    double sub = lceamount*(mapval-factor*(r*lumimulr + g*lumimulg + b*lumimulb));
                    r = factor*r-sub;
                    g = factor*g-sub;
                    b = factor*b-sub;
                }
                else {
                    r *= factor;
                    g *= factor;
                    b *= factor;
                }
                targetImage->r[i][j] = r;
                targetImage->g[i][j] = g;
                targetImage->b[i][j] = b;
            }
        }
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}

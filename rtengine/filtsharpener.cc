/*
 * filtsharpener.cc
 *
 *  Created on: Sep 16, 2010
 *      Author: gabor
 */

#include "filtsharpener.h"
#include "rtengine.h"
#include "macros.h"
#include "gauss.h"
#include "bilateral2.h"
#include "minmax.h"

namespace rtengine {

SharpenFilterDescriptor sharpenFilterDescriptor;

SharpenFilterDescriptor::SharpenFilterDescriptor ()
	: FilterDescriptor ("Sharpener", MultiImage::Lab, MultiImage::Lab) {

    addTriggerEvent (EvShrEnabled);
    addTriggerEvent (EvShrRadius);
    addTriggerEvent (EvShrAmount);
    addTriggerEvent (EvShrThresh);
    addTriggerEvent (EvShrEdgeOnly);
    addTriggerEvent (EvShrEdgeRadius);
    addTriggerEvent (EvShrEdgeTolerance);
    addTriggerEvent (EvShrHaloControl);
    addTriggerEvent (EvShrHaloAmount);
    addTriggerEvent (EvShrMethod);
    addTriggerEvent (EvShrDRadius);
    addTriggerEvent (EvShrDAmount);
    addTriggerEvent (EvShrDDamping);
    addTriggerEvent (EvShrDIterations);

    applyOnThumbnail = false;
}

void SharpenFilterDescriptor::getDefaultParameters (ProcParams& defProcParams) const {

	defProcParams.setBoolean ("Sharpening", "Enabled", true);
	defProcParams.setFloat   ("Sharpening", "Radius", 1.0);
	defProcParams.setFloat   ("Sharpening", "Amount", 90);
	defProcParams.setFloat   ("Sharpening", "Threshold", 1.5);
	defProcParams.setBoolean ("Sharpening", "EdgesOnly", false);
	defProcParams.setFloat   ("Sharpening", "EdgesRadius", 3);
	defProcParams.setFloat   ("Sharpening", "EdgesTolerance", 1.5);
	defProcParams.setBoolean ("Sharpening", "HaloControl", false);
	defProcParams.setFloat   ("Sharpening", "HaloControlAmount", 85);
	defProcParams.setString  ("Sharpening", "Method", "usm");
	defProcParams.setFloat   ("Sharpening", "DeconvAmount", 75);
	defProcParams.setFloat   ("Sharpening", "DeconvRadius", 0.75);
	defProcParams.setInteger ("Sharpening", "DeconvIter", 30);
	defProcParams.setFloat   ("Sharpening", "DeconvDamping", 20);
}

void SharpenFilterDescriptor::createAndAddToList (Filter* tail) const {

	tail->addNext (new SharpenFilter ());
}

SharpenFilter::SharpenFilter ()
	: Filter (&sharpenFilterDescriptor) {
}

void SharpenFilter::dcdamping (Buffer<float>* aI, MultiImage* aO, float damping) {

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<aI->height; i++)
        for (int j=0; j<aI->width; j++) {
            float I = aI->rows[i][j];
            float O = (float)(aO->cieL[i][j]);
            if (O==0.0 || I==0.0) {
                aI->rows[i][j] = 0.0;
                continue;
            }
            float U = -(O * log(I/O) - I + O) * 2.0 / (damping*damping);
            U = MIN(U,1.0);
            U = U*U*U*U*(5.0-U*4.0);
            aI->rows[i][j] = (O - I) / I * U + 1.0;
        }
}

void SharpenFilter::deconvsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* b2) {

    int W = sourceImage->width, H = sourceImage->height;
    Dim size (W, H);
    double scale = getScale ();

	float amount = procParams->getFloat ("Sharpening", "DeconvAmount") / 100.0;
	float radius = procParams->getFloat ("Sharpening", "DeconvRadius") / scale;
	int   iter   = procParams->getInteger ("Sharpening", "DeconvIter");
	float damping = procParams->getFloat ("Sharpening", "DeconvDamping") / 5.0;



    Buffer<float>* tmpI = new Buffer<float> (W, H);
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            tmpI->rows[i][j] = sourceImage->cieL[i][j];

    Buffer<float>* tmp = b2;

    #ifdef _OPENMP
    double* buffer = new double [std::max(W,H)*omp_get_max_threads()];
	#else
    double* buffer = new double [std::max(W,H)];
	#endif

    bool needdamp = damping > 0;
    for (int k=0; k<iter; k++) {

        // apply blur function (gaussian blur)
        gaussHorizontal<float> (tmpI, tmp, size, buffer, radius, multiThread);
        gaussVertical<float>   (tmp, tmp,  size, buffer, radius, multiThread);

        if (!needdamp) {
            #pragma omp parallel for if (multiThread)
            for (int i=0; i<H; i++)
                for (int j=0; j<W; j++)
                    if (tmp->rows[i][j]>0)
                        tmp->rows[i][j] = (float)sourceImage->cieL[i][j] / tmp->rows[i][j];
        }
        else
            dcdamping (tmp, sourceImage, damping);

        gaussHorizontal<float> (tmp, tmp, size, buffer, radius, multiThread);
        gaussVertical<float>   (tmp, tmp, size, buffer, radius, multiThread);

        #pragma omp parallel for if (multiThread)
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                tmpI->rows[i][j] = tmpI->rows[i][j] * tmp->rows[i][j];
    }
    delete [] buffer;

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            targetImage->cieL[i][j] = sourceImage->cieL[i][j]*(1.0-amount) + tmpI->rows[i][j]*amount;
            targetImage->ciea[i][j] = sourceImage->ciea[i][j];
            targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        }
    delete tmpI;
}

void SharpenFilter::sharpenHaloCtrl (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* blurmap, Buffer<float>* base) {

	float amount    = procParams->getFloat ("Sharpening", "Amount") / 100;
	float threshold = procParams->getFloat ("Sharpening", "Threshold") / 100;
	float scale     = 1.0 - procParams->getFloat ("Sharpening", "HaloControlAmount") / 100;

    float** nL = base->rows;
    int W = blurmap->width;
    int H = blurmap->height;
    #pragma omp parallel for if (multiThread)
    for (int i=2; i<H-2; i++) {
        float max1 = 0, max2 = 0, min1 = 0, min2 = 0, maxn, minn, np1, np2, np3, min, max;
        for (int j=2; j<W-2; j++) {
            if (i>1 && j>1 && i<H-2 && j<W-2) {
                float diff = base->rows[i][j] - blurmap->rows[i][j];
                if (fabs(diff) > threshold) {
                    // compute maximum/minimum in a delta environment
                    np1 = 2*(nL[i-2][j] + nL[i-2][j+1] + nL[i-2][j+2] + nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2]) / 27 + nL[i-1][j+1] / 3;
                    np2 = 2*(nL[i-1][j] + nL[i-1][j+1] + nL[i-1][j+2] + nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2]) / 27 + nL[i][j+1] / 3;
                    np3 = 2*(nL[i][j] + nL[i][j+1] + nL[i][j+2] + nL[i+1][j] + nL[i+1][j+1] + nL[i+1][j+2] + nL[i+2][j] + nL[i+2][j+1] + nL[i+2][j+2]) / 27 + nL[i+1][j+1] / 3;
                    MINMAX3(np1,np2,np3,maxn,minn);
                    MAX3(max1,max2,maxn,max);
                    MIN3(min1,min2,minn,min);
                    max1 = max2; max2 = maxn;
                    min1 = min2; min2 = minn;
                    if (max < sourceImage->cieL[i][j])
                        max = sourceImage->cieL[i][j];
                    if (min > sourceImage->cieL[i][j])
                        min = sourceImage->cieL[i][j];
                    float val = sourceImage->cieL[i][j] + amount * diff;
                    // applying halo control
                    if (val > max)
                        val = max + (val-max) * scale;
                    else if (val < min)
                        val = min - (min-val) * scale;
                    targetImage->cieL[i][j] = val;
                }
            }
            else
                targetImage->cieL[i][j] = sourceImage->cieL[i][j];
        }
    }
}

void SharpenFilter::usmsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* b2) {

    int W = sourceImage->width, H = sourceImage->height;
    Dim size (W, H);
    double scale = getScale ();

	float amount     = procParams->getFloat ("Sharpening", "Amount") / 100;
	float threshold  = procParams->getFloat ("Sharpening", "Threshold") / 100;
	bool edgesonly   = procParams->getBoolean ("Sharpening", "EdgesOnly");
	float radius     = procParams->getFloat ("Sharpening", "Radius") / scale;
	bool halocontrol = procParams->getBoolean ("Sharpening", "HaloControl");
	float eradius    = procParams->getFloat ("Sharpening", "EdgesRadius") / scale;
	float etolerance = procParams->getFloat   ("Sharpening", "EdgesTolerance") / 100;

    Buffer<float>* b3;
    Buffer<float> sourceView = sourceImage->getBufferView (sourceImage->cieL);

    #ifdef _OPENMP
    double* buffer = new double [std::max(W,H)*omp_get_max_threads()];
	#else
    double* buffer = new double [std::max(W,H)];
	#endif

    if (!edgesonly) {
        gaussHorizontal<float> (&sourceView, b2, size, buffer, radius, multiThread);
        gaussVertical<float>   (b2,          b2, size, buffer, radius, multiThread);
    }
    else {
        b3 = new Buffer<float> (W, H);
        bilateral (sourceImage->cieL, b3->rows, b2->rows, W, H, eradius, etolerance, multiThread);
        gaussHorizontal<float> (b3, b2, size, buffer, radius / scale, multiThread);
        gaussVertical<float>   (b2, b2, size, buffer, radius / scale, multiThread);
    }
    delete [] buffer;

    Buffer<float>* base = &sourceView;
    if (edgesonly)
        base = b3;

    if (!halocontrol) {
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++) {
                float diff = base->rows[i][j] - b2->rows[i][j];
                if (fabs(diff) > threshold)
                	targetImage->cieL[i][j] = sourceImage->cieL[i][j] + amount * diff;
            }
    }
    else
        sharpenHaloCtrl (sourceImage, targetImage, b2, base);

    if (edgesonly)
        delete b3;
}

Dim SharpenFilter::getReqiredBufferSize () {

    return getScaledTargetImageView().getSize();
}

void SharpenFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<float>* buffer) {

	bool enabled  = procParams->getBoolean ("Sharpening", "Enabled");
	String method = procParams->getString ("Sharpening", "Method");
	float uamount = procParams->getFloat ("Sharpening", "Amount");
	float damount = procParams->getFloat ("Sharpening", "DeconvAmount");

    if (getTargetImageView().skip==1 && enabled && method=="rld" && damount>0 && sourceImage->width>=8 && sourceImage->height>=8) {

        deconvsharpening (sourceImage, targetImage, (Buffer<float>*)buffer);

        if (sourceImage != targetImage)
        	for (int i=0; i<targetImage->height; i++)
        		for (int j=0; j<targetImage->width; j++) {
        			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
        			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        		}
    }
    else if (getTargetImageView().skip==1 && enabled && method=="usm" && uamount>0 && sourceImage->width>=8 && sourceImage->height>=8) {

    	usmsharpening (sourceImage, targetImage, (Buffer<float>*)buffer);

        if (sourceImage != targetImage)
        	for (int i=0; i<targetImage->height; i++)
        		for (int j=0; j<targetImage->width; j++) {
        			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
        			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        		}
    }
    else if (sourceImage != targetImage)
        targetImage->copyFrom (sourceImage);
}

}

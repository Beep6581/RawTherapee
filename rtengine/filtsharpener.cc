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

    Buffer<float>* tmpI = new Buffer<float> (W, H);
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++)
            tmpI->rows[i][j] = (float)sourceImage->cieL[i][j];

    Buffer<float>* tmp = b2;

    double* buffer = new double [std::max(W,H)*omp_get_max_threads()];

    float damping = procParams->sharpening.deconvdamping / 5.0;
    bool needdamp = procParams->sharpening.deconvdamping > 0;
    for (int k=0; k<procParams->sharpening.deconviter; k++) {

        // apply blur function (gaussian blur)
        gaussHorizontal<float> (tmpI, tmp, size, buffer, procParams->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp,  size, buffer, procParams->sharpening.deconvradius / scale, multiThread);

        if (!needdamp) {
            #pragma omp parallel for if (multiThread)
            for (int i=0; i<H; i++)
                for (int j=0; j<W; j++)
                    if (tmp->rows[i][j]>0)
                        tmp->rows[i][j] = (float)sourceImage->cieL[i][j] / tmp->rows[i][j];
        }
        else
            dcdamping (tmp, sourceImage, damping);

        gaussHorizontal<float> (tmp, tmp, size, buffer, procParams->sharpening.deconvradius / scale, multiThread);
        gaussVertical<float>   (tmp, tmp, size, buffer, procParams->sharpening.deconvradius / scale, multiThread);

        #pragma omp parallel for if (multiThread)
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++)
                tmpI->rows[i][j] = tmpI->rows[i][j] * tmp->rows[i][j];
    }
    delete [] buffer;

    #pragma omp parallel for if (multiThread)
    for (int i=0; i<H; i++)
        for (int j=0; j<W; j++) {
            targetImage->cieL[i][j] = sourceImage->cieL[i][j]*(100-procParams->sharpening.deconvamount) / 100 + (int)CLIP(tmpI->rows[i][j])*procParams->sharpening.deconvamount / 100;
            targetImage->ciea[i][j] = sourceImage->ciea[i][j];
            targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        }
    delete tmpI;
}

void SharpenFilter::sharpenHaloCtrl (MultiImage* sourceImage, MultiImage* targetImage, Buffer<unsigned short>* blurmap, Buffer<unsigned short>* base) {

    int scale = 100 * (100 - procParams->sharpening.halocontrol_amount);
    unsigned short** nL = base->rows;
    int W = blurmap->width;
    int H = blurmap->height;
    #pragma omp parallel for if (multiThread)
    for (int i=2; i<H-2; i++) {
        int max1 = 0, max2 = 0, min1 = 0, min2 = 0, maxn, minn, np1, np2, np3, min, max;
        for (int j=2; j<W-2; j++) {
            if (i>1 && j>1 && i<H-2 && j<W-2) {
                int diff = base->rows[i][j] - blurmap->rows[i][j];
                if (ABS(diff) > procParams->sharpening.threshold) {
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
                    int val = sourceImage->cieL[i][j] + procParams->sharpening.amount * diff / 100;
                    int newL = CLIP(val);
                    // applying halo control
                    if (newL > max)
                        newL = max + (newL-max) * scale / 10000;
                    else if (newL<min)
                        newL = min - (min-newL) * scale / 10000;
                    targetImage->cieL[i][j] = newL;
                }
            }
            else
                targetImage->cieL[i][j] = sourceImage->cieL[i][j];
        }
    }
}

void SharpenFilter::usmsharpening (MultiImage* sourceImage, MultiImage* targetImage, Buffer<unsigned short>* b2) {

    int W = sourceImage->width, H = sourceImage->height;
    Dim size (W, H);
    double scale = getScale ();

    Buffer<unsigned short>* b3;
    Buffer<unsigned short> sourceView = sourceImage->getBufferView (sourceImage->cieL);

    double* buffer = new double [std::max(W, H) * omp_get_max_threads()];
    if (procParams->sharpening.edgesonly==false) {
        gaussHorizontal<unsigned short> (&sourceView, b2, size, buffer, procParams->sharpening.radius / scale, multiThread);
        gaussVertical<unsigned short>   (b2,          b2, size, buffer, procParams->sharpening.radius / scale, multiThread);
    }
    else {
        b3 = new Buffer<unsigned short> (W, H);
        bilateral<unsigned short, unsigned int> (sourceImage->cieL, b3->rows, b2->rows, W, H, procParams->sharpening.edges_radius / scale, procParams->sharpening.edges_tolerance, multiThread);
        gaussHorizontal<unsigned short> (b3, b2, size, buffer, procParams->sharpening.radius / scale, multiThread);
        gaussVertical<unsigned short>   (b2, b2, size, buffer, procParams->sharpening.radius / scale, multiThread);
    }
    delete [] buffer;

    Buffer<unsigned short>* base = &sourceView;
    if (procParams->sharpening.edgesonly)
        base = b3;

    if (procParams->sharpening.halocontrol==false) {
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<H; i++)
            for (int j=0; j<W; j++) {
                int diff = base->rows[i][j] - b2->rows[i][j];
                if (ABS(diff) > procParams->sharpening.threshold) {
                    int val = sourceImage->cieL[i][j] + procParams->sharpening.amount * diff / 100;
                    targetImage->cieL[i][j] = CLIP(val);
                }
            }
    }
    else
        sharpenHaloCtrl (sourceImage, targetImage, b2, base);

    if (procParams->sharpening.edgesonly)
        delete b3;
}

Dim SharpenFilter::getReqiredBufferSize () {

    return getScaledTargetImageView().getSize();
}

void SharpenFilter::process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer) {

    if (getTargetImageView().skip==1 && procParams->sharpening.enabled && procParams->sharpening.method=="rld" && procParams->sharpening.deconvamount>0 && sourceImage->width>=8 && sourceImage->height>=8) {

        deconvsharpening (sourceImage, targetImage, (Buffer<float>*)buffer);

        if (sourceImage != targetImage)
        	for (int i=0; i<targetImage->height; i++)
        		for (int j=0; j<targetImage->width; j++) {
        			targetImage->ciea[i][j] = sourceImage->ciea[i][j];
        			targetImage->cieb[i][j] = sourceImage->cieb[i][j];
        		}
    }
    else if (getTargetImageView().skip==1 && procParams->sharpening.enabled && procParams->sharpening.method=="usm" && procParams->sharpening.amount>0 && sourceImage->width>=8 && sourceImage->height>=8) {

        usmsharpening (sourceImage, targetImage, (Buffer<unsigned short>*)buffer);

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

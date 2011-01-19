#include "filterchain.h"
#include "imageview.h"
#include "settings.h"
#include "filterfactory.h"
#include "filtinitial.h"
#include "image16.h"
#include "image8.h"
#include "iccstore.h"
#include "curves.h"
#include <lcms2.h>
#include <iostream>

namespace rtengine {

FilterChain::FilterChain (ImProcListener* listener, ImageSource* imgSource, ProcParams* params, bool multiThread)
	: listener(listener), imgSource(imgSource), procParams(params), multiThread(multiThread), invalidated(true) {

    last = first = new InitialFilter (imgSource);

	if (params->getBoolean ("FilterOrderCustom") && !params->getStringList ("FilterOrderFilterList").empty())
		filterOrder = params->getStringList ("FilterOrderFilterList");
	else
		filterOrder = Settings::settings->filterList;
	setupChain (NULL);
}

FilterChain::FilterChain (ImProcListener* listener, FilterChain* previous)
	: listener(listener), imgSource(previous->imgSource),
	  procParams(previous->procParams), multiThread(previous->multiThread), invalidated(true) {

    last = first = new InitialFilter (imgSource);

    filterOrder = previous->filterOrder;
	setupChain (previous);
}

void FilterChain::setupChain (FilterChain* previous) {

    // First filter already there, it is the initial filter. We go on from the second one.
	Filter* parent = previous ? previous->first->next : NULL;
	for (int i=0; i<filterOrder.size(); i++)
		if (filterOrder[i] == "--Cache--")
			last->forceOutputCache = true;
		else {
		    FilterDescriptor* fDescr = filterFactory.getFilterDescriptor (filterOrder[i]);
		    if (fDescr && (
		            (fDescr->isAppliedOnThumbnail() && imgSource->isThumbnail()) ||
		            (fDescr->isAppliedOnRawImage() && imgSource->isRaw()) ||
		            (fDescr->isAppliedOnStdImage() && !imgSource->isRaw()))) {
                fDescr->createAndAddToList (last);
                while (last->next) {
                    last = last->next;
                    if (parent) {
                        parent = parent->next;
                        last->parent = parent;
                    }
                    last->multiThread = multiThread;
                    last->myFilterChain = this;
                }
		    }
		}
	for (Filter* curr = first; curr; curr = curr->next) {
	    curr->setProcParams (procParams);
	    curr->myFilterChain = this;
	}
}

void FilterChain::setNextChain (FilterChain* other) {

    if (!other)
        for (Filter* curr = first; curr; curr = curr->next)
            curr->parent = NULL;
    else
        for (Filter *curr = first, *pcurr = other->first; curr; curr = curr->next, pcurr = pcurr->next)
            curr->parent = pcurr;
}

FilterChain::~FilterChain () {

	Filter* curr = first;
	while (curr) {
		Filter* tmp = curr;
		curr = curr->next;
		delete tmp;
	}
}

void FilterChain::invalidate () {

    invalidated = true;
}

void FilterChain::setupProcessing (const std::set<ProcEvent>& events, bool useShortCut) {

	// tell the listener the full size of the image with the current settings and let it decide the portion to refresh
	Dim fullSize = getFullImageSize ();
	ImageView reqView;
	if (listener)
	    reqView = listener->getViewToProcess (fullSize);
	else {
        reqView.w = fullSize.width;
        reqView.h = fullSize.height;
	}
	// walk through the list and find first filter that needs to be refreshed because of the changed view
	Filter* curr = last;
	firstToUpdate = last;
	ImageView view = reqView;
	while (curr) {
	    ImageView tiv = curr->calculateTargetImageView (view);
        if (curr->targetImageView != tiv)
            firstToUpdate = curr;
        curr->targetImageView = tiv;
        view  = curr->calculateSourceImageView (view);
        curr->sourceImageView = view;
        curr = curr->prev;
	}

	// calculate and set source and target image sizes of the filters
	for (curr = first; curr; curr = curr->next) {
        ImageView scSourceIV = curr->sourceImageView.getScaled (curr==first ? imgSource->getScale () : curr->prev->getScale ());
        ImageView scTargetIV = curr->targetImageView.getScaled (curr->getScale ());
        if (scSourceIV != curr->scaledSourceImageView || scTargetIV != curr->scaledTargetImageView) {
            curr->scaledSourceImageView = scSourceIV;
            curr->scaledTargetImageView = scTargetIV;
            curr->valid = false;
        }
	}

	// set mandatory cache points
	for (curr = first->next; curr != last; curr = curr->next)
		if (curr->forceOutputCache || curr->sourceImageView != curr->targetImageView
		        || curr->scaledSourceImageView != curr->scaledTargetImageView
		        || (curr->next && (curr->targetImageView != curr->next->sourceImageView || curr->scaledTargetImageView != curr->next->scaledSourceImageView)))
			curr->hasOutputCache = true;
		else
			curr->hasOutputCache = false;
	last->hasOutputCache = true;

	// walk further to find first filter that needs to be recalculated
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (invalidated || !curr->valid || curr->isTriggerEvent (events))
			firstToUpdate = curr;

	// walk further to find closest cache
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (!curr->prev || curr->prev->hasOutputCache) {
			firstToUpdate = curr;
			break;
		}

	// set every filter from firstToUpdate to be invalid
	for (curr = firstToUpdate; curr; curr = curr->next)
	    curr->valid = false;

	// detect shortcut possibilities in the filter group
	if (useShortCut) {
	    // clear shortcut pointers
	    for (curr = first; curr; curr = curr->next)
	        curr->shortCutPrev = NULL;
	    // find possible shortcut filter
	    for (curr = last; curr && curr != firstToUpdate; curr = curr->prev) {
	        Filter* p = curr->prev->parent;
	        while (p) {
                double skip = curr->getScale() / p->getScale();
	            if (p->hasOutputCache && fabs(skip-round(skip))<1e-12 && curr->sourceImageView.isPartOf (p->targetImageView)) {
	                curr->shortCutPrev = p;
	                firstToUpdate = curr;
	                break;
	            }
	            p = p->parent;
	        }
	        if (curr->shortCutPrev)
	            break;
	    }
	}

	// set up caches
	for (curr = first; curr; curr = curr->next)
		curr->setupCache ();
}

void FilterChain::process (const std::set<ProcEvent>& events, Buffer<float>* buffer, MultiImage* worker) {

    for (Filter* curr = firstToUpdate; curr; curr = curr->next) {

	    if (Settings::settings->verbose)
	        std::cout << "Applying filter " << curr->getDescriptor()->getName() << " "
                << curr->sourceImageView << ", " << curr->targetImageView << " in progress...";
	    std::flush (std::cout);

	    MultiImage* sourceImage = NULL;
		MultiImage* targetImage = NULL;

		// set up source image for the filter
		if (curr != first) {
		    if (curr->shortCutPrev) {
		        sourceImage = worker;
                worker->setDimensions (curr->scaledSourceImageView.w, curr->scaledSourceImageView.h);
                int skip = (int)round (curr->getScale() / curr->shortCutPrev->getScale());
                worker->copyFrom (curr->shortCutPrev->outputCache, curr->scaledSourceImageView.x - curr->shortCutPrev->scaledTargetImageView.x * skip, curr->scaledSourceImageView.y - curr->shortCutPrev->scaledTargetImageView.y * skip, skip);
		    }
		    else if (curr->prev->hasOutputCache && curr->prev->targetImageView == curr->sourceImageView && curr->prev->scaledTargetImageView == curr->scaledSourceImageView)
				sourceImage = curr->prev->outputCache;
			else {
				sourceImage = worker;
				worker->setDimensions (curr->scaledSourceImageView.w, curr->scaledSourceImageView.h);
				if (curr->prev->hasOutputCache && curr->prev->targetImageView != curr->sourceImageView || curr->prev->scaledTargetImageView != curr->scaledSourceImageView) {
				    // There is a cache, but image views do not fit. Assume that sourceImageView is the part of prev->targetImageView.
                    int skip = (int)round (curr->getScale() / curr->prev->getScale());
				    worker->copyFrom (curr->prev->outputCache, curr->scaledSourceImageView.x - curr->prev->scaledTargetImageView.x * skip, curr->scaledSourceImageView.y - curr->prev->scaledTargetImageView.y * skip, skip);
				}
			}
			sourceImage->convertTo (curr->descriptor->getInputColorSpace());
		}

		// set up target image for the filter
		if (curr->hasOutputCache)
			targetImage = curr->outputCache;
		else {
			worker->setDimensions (curr->scaledTargetImageView.w, curr->scaledTargetImageView.h);
			targetImage = worker;
		}
		targetImage->switchTo (curr->descriptor->getOutputColorSpace());

		// apply filter
		curr->process (events, sourceImage, targetImage, buffer);
		curr->valid = true;

		if (Settings::settings->verbose)
            std::cout << "ready." << std::endl;
	}
	invalidated = false;
}

Dim FilterChain::getReqiredBufferSize () {

	Dim bufferSize;
	for (Filter* curr = first; curr; curr = curr->next)
		bufferSize.setMax (curr->getReqiredBufferSize ());

	return bufferSize;
}

Dim FilterChain::getReqiredWorkerSize () {

    Dim workerSize;
    for (Filter* curr = firstToUpdate; curr; curr = curr->next)
        if (!curr->hasOutputCache || (curr->prev && (curr->prev->targetImageView != curr->sourceImageView || curr->prev->scaledTargetImageView != curr->scaledSourceImageView)))
            workerSize.setMax (curr->scaledSourceImageView.getSize ());

    return workerSize;
}

Dim FilterChain::getFullImageSize () {

	if (!last)
		return Dim();
	else
		return last->getFullImageSize ();
}

// this method is called by the gui to determine the image view to request from rtengine such that
// the size of the result is given by a gui window
double FilterChain::getScale (int skip) {

    if (!last)
        return 1.0;
    else
        return last->getTargetScale (skip);
}

ImageView FilterChain::getLastImageView () {

    if (last)
        return last->targetImageView;
    else
        return ImageView();
}

double FilterChain::getLastScale () {

    return last ? last->getScale() : 1.0;
}

Image16* FilterChain::getFinalImage () {

    // TODO! we can do better! instead of converting to rgb, let us convert it to xyz
    last->outputCache->convertTo(MultiImage::RGB, true, procParams->icm.working);

    // calculate crop rectangle
    int cx = procParams->crop.enabled ? procParams->crop.x : 0;
    int cy = procParams->crop.enabled ? procParams->crop.y : 0;
    int cw = procParams->crop.enabled ? procParams->crop.w : last->outputCache->width;
    int ch = procParams->crop.enabled ? procParams->crop.h : last->outputCache->height;

    // adjust to valid values
    if (cx<0) cx = 0;
    if (cy<0) cy = 0;
    if (cx+cw > last->outputCache->width)  cw = last->outputCache->width - cx;
    if (cy+ch > last->outputCache->height) ch = last->outputCache->height - cy;

    // obtain cropped image
    Image16* final = new Image16 (cw, ch);
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<ch; i++)
        for (int j=0; j<cw; j++) {
            final->r[i][j] = last->outputCache->r[i+cy][j+cx];
            final->g[i][j] = last->outputCache->g[i+cy][j+cx];
            final->b[i][j] = last->outputCache->b[i+cy][j+cx];
        }

    if (procParams->icm.output!="") {
        // custom output icm profile specified, call lcms
        cmsHPROFILE oprof = iccStore.getProfile (procParams->icm.output);
        cmsHPROFILE iprof = iccStore.getProfile (procParams->icm.working);
        if (oprof && iprof) {
        	// TODO: do not create cmstransform every time, store it and re-create it only when it has been changed
        	cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_16_PLANAR, Settings::settings->colorimetricIntent, 0);
        	cmsDoTransform (hTransform, final->data, final->data, final->planestride/2);
        	cmsDeleteTransform(hTransform);
        	return final;
        }
    }
    else if (procParams->icm.working != "sRGB") {
        // convert it to sRGB
        Matrix33 m;
        m.multiply (iccStore.workingSpaceMatrix ("sRGB"));
        m.multiply (iccStore.workingSpaceInverseMatrix (procParams->icm.working));
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<ch; i++)
            for (int j=0; j<cw; j++)
                m.transform (final->r[i][j], final->g[i][j], final->b[i][j]);
    }

    // apply gamma correction
    #pragma omp parallel for if (multiThread)
    for (int i=0; i<ch; i++)
        for (int j=0; j<cw; j++) {
            final->r[i][j] = CurveFactory::gamma_srgb (final->r[i][j]);
            final->g[i][j] = CurveFactory::gamma_srgb (final->g[i][j]);
            final->b[i][j] = CurveFactory::gamma_srgb (final->b[i][j]);
        }
    return final;
}

Image8* FilterChain::getDisplayImage () {

    // TODO! we can do better! instead of converting to rgb, let us convert it to xyz
    last->outputCache->convertTo(MultiImage::RGB, true, procParams->icm.working);

    // obtain cropped image
    Image8* final = new Image8 (last->outputCache->width, last->outputCache->height);

    if (Settings::settings->monitorProfile != "") {
        // custom output icm profile specified, call lcms
        cmsHPROFILE oprof = iccStore.getProfile (Settings::settings->monitorProfile);
        cmsHPROFILE iprof = iccStore.getProfile (procParams->icm.working);
        if (oprof && iprof) {
			// TODO: do not create cmstransform every time, store it and re-create it only when it has been changed
			cmsHTRANSFORM hTransform = cmsCreateTransform (iprof, TYPE_RGB_16_PLANAR, oprof, TYPE_RGB_8, Settings::settings->colorimetricIntent, 0);
			cmsDoTransform (hTransform, last->outputCache->getData(), final->data, final->width*final->height);
			cmsDeleteTransform(hTransform);
			return final;
        }
    }

    if (procParams->icm.working == "sRGB") {
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<final->height; i++) {
            int ix = i * final->width * 3;
            for (int j=0; j<final->width; j++) {
                final->data[ix++] = CurveFactory::gamma_srgb (last->outputCache->r[i][j]) >> 8;
                final->data[ix++] = CurveFactory::gamma_srgb (last->outputCache->g[i][j]) >> 8;
                final->data[ix++] = CurveFactory::gamma_srgb (last->outputCache->b[i][j]) >> 8;
            }
        }
    }
    else {
        Matrix33 m;
        m.multiply (iccStore.workingSpaceMatrix ("sRGB"));
        m.multiply (iccStore.workingSpaceInverseMatrix (procParams->icm.working));
        #pragma omp parallel for if (multiThread)
        for (int i=0; i<final->height; i++) {
            int ix = i * final->width * 3;
            for (int j=0; j<final->width; j++) {
                unsigned short r, g, b;
                m.transform (last->outputCache->r[i][j], last->outputCache->g[i][j], last->outputCache->b[i][j], r, g, b);
                final->data[ix++] = CurveFactory::gamma_srgb (r) >> 8;
                final->data[ix++] = CurveFactory::gamma_srgb (g) >> 8;
                final->data[ix++] = CurveFactory::gamma_srgb (b) >> 8;
            }
        }
    }
    return final;
}

}

#include "filterchain.h"
#include "imageview.h"
#include "settings.h"
#include "filterfactory.h"
#include "filtinitial.h"

namespace rtengine {

FilterChain::FilterChain (ImProcListener* listener, ImageSource* imgSource, ProcParams* params, bool multiThread)
	: listener(listener), imgSource(imgSource), procParams(params), multiThread(multiThread), invalidated(true) {

    last = first = new InitialFilter (imgSource);

	if (params->filterOrder.custom && !params->filterOrder.filterlist.empty())
		filterOrder = params->filterOrder.filterlist;
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

void FilterChain::setupProcessing (const std::set<ProcEvent>& events, Dim fullSize, Dim& maxWorkerSize, bool useShortCut) {

	if (!listener)
		return;

	// tell the listener the full size of the image with the current settings and let it decide the portion to refresh
	ImageView reqView = listener->getViewToProcess (fullSize);

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
        ImageView scSourceIV = curr->sourceImageView.getScaled (first ? imgSource->getScale () : curr->prev->getScale ());
        ImageView scTargetIV = curr->targetImageView.getScaled (curr->getScale ());
        if (scSourceIV != curr->scaledSourceImageView || scTargetIV != curr->scaledTargetImageView) {
            curr->scaledSourceImageView = scSourceIV;
            curr->scaledTargetImageView = scTargetIV;
            curr->valid = false;
        }
	}

	// set mandatory cache points
	first->hasOutputCache = true;
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

	// find out the dimensions of the largest worker image necessary
	for (curr = firstToUpdate; curr; curr = curr->next)
		if (!curr->hasOutputCache)
		    maxWorkerSize.setMax (curr->scaledSourceImageView.getSize ());

	// set up caches
	for (curr = first; curr; curr = curr->next)
		curr->setupCache ();
}

void FilterChain::process (const std::set<ProcEvent>& events, Buffer<int>* buffer, MultiImage* worker) {

	for (Filter* curr = firstToUpdate; curr; curr = curr->next) {
		MultiImage* sourceImage = NULL;
		MultiImage* targetImage = NULL;
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
                int skip = (int)round (curr->getScale() / curr->shortCutPrev->getScale());
				if (curr->prev->hasOutputCache && curr->prev->targetImageView != curr->sourceImageView || curr->prev->scaledTargetImageView != curr->scaledSourceImageView)
				    // There is a cache, but image views do not fit. Assume that sourceImageView is the part of prev->targetImageView.
				    worker->copyFrom (curr->prev->outputCache, curr->scaledSourceImageView.x - curr->prev->scaledTargetImageView.x * skip, curr->scaledSourceImageView.y - curr->prev->scaledTargetImageView.y * skip, skip);
			}
			sourceImage->convertTo (curr->descriptor->getInputColorSpace());
		}
		if (curr->hasOutputCache)
			targetImage = curr->outputCache;
		else {
			worker->setDimensions (curr->scaledTargetImageView.w, curr->scaledTargetImageView.h);
			targetImage = worker;
		}
		curr->process (events, sourceImage, targetImage, buffer);
		curr->valid = true;
	}
	invalidated = false;
}

Dim FilterChain::getReqiredBufferSize () {

	Dim bufferSize;
	for (Filter* curr = first; curr; curr = curr->next)
		bufferSize.setMax (curr->getReqiredBufferSize ());
}

Dim FilterChain::getFullImageSize () {

	if (!last)
		return Dim();
	else
		return last->getFullImageSize ();
}

double FilterChain::getScale (int skip) {

    return last ? last->getScale() : 1.0;
}

}

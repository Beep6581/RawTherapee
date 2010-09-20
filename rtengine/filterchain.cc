#include "filterchain.h"
#include "imageview.h"
#include "settings.h"
#include "filterfactory.h"

namespace rtengine {
using namespace procparams;

FilterChain::FilterChain (ImProcListener* listener, ImageSource* imgSource, ProcParams* params, bool multiThread)
	: listener(listener), first(imgSource), last(imgSource), procParams(params), multiThread(multiThread) {

	if (params->filterOrder.custom && !params->filterOrder.filterlist.empty())
		filterOrder = params->filterOrder.filterlist;
	else
		filterOrder = imgSource->isRaw() ? Settings::settings->filterListRawImage : Settings::settings->filterListStdImage;
	setupChain (NULL);
}

FilterChain::FilterChain (ImProcListener* listener, FilterChain* previous)
	: listener(listener), first(previous->first), last(previous->first), procParams(previous->procParams), multiThread(previous->multiThread) {

	filterOrder = previous->filterOrder;
	setupChain (previous);
}

void FilterChain::setupChain (FilterChain* previous) {

	Filter* parent = previous ? previous->first : NULL;
	for (int i=0; i<filterOrder.size(); i++)
		if (filterOrder[i] == "--Cache--")
			last->forceOutputCache = true;
		else {
			filterFactory.createFilterAddToList (filterOrder[i], last);
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

FilterChain::~FilterChain () {

	Filter* curr = first;
	while (curr) {
		Filter* tmp = curr;
		curr = curr->next;
		delete tmp;
	}
}

void FilterChain::invalidate () {

	// TODO!!! IMPORTANT!!
}

void FilterChain::setupProcessing (const std::set<ProcEvent>& events, int fullW, int fullH, int& maxWorkerWidth, int& maxWorkerHeight) {

	if (!listener || !first || !last)
		return;

	// tell the listener the full size of the image with the current settings and let it decide the portion to refresh
	ImageView reqView = listener->getViewToProcess (fullW, fullH);

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

	// set mandatory cache points
	first->hasOutputCache = true;
	for (curr = first->next; curr != last; curr = curr->next)
		if (curr->forceOutputCache || !(curr->sourceImageView == curr->targetImageView) || (curr->next && !(curr->targetImageView == curr->next->sourceImageView)))
			curr->hasOutputCache = true;
		else
			curr->hasOutputCache = false;
	last->prev->hasOutputCache = true;	// cache is needed in the last phase since output format is different
	last->hasOutputCache = false;

	// walk further to find first filter that needs to be recalculated
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (curr->isTriggerEvent (events))
			firstToUpdate = curr;

	// walk further to find closest cache
	for (curr = firstToUpdate; curr; curr = curr->prev)
		if (!curr->prev || curr->prev->hasOutputCache) {
			firstToUpdate = curr;
			break;
		}

	// find out the dimensions of the largest worker image necessary
	for (curr = first->next; curr; curr = curr->next)
		if (!curr->hasOutputCache) {
			if (curr->sourceImageView.getPixelWidth() > maxWorkerWidth)
				maxWorkerWidth = curr->sourceImageView.getPixelWidth ();
			if (curr->sourceImageView.getPixelHeight() > maxWorkerHeight)
				maxWorkerHeight = curr->sourceImageView.getPixelHeight ();
		}

	// set up caches
	for (curr = first; curr; curr = curr->next)
		curr->setupCache ();
}

void FilterChain::process (Buffer<int>* buffer, MultiImage* worker) {

	for (Filter* curr = firstToUpdate; curr; curr = curr->next) {
		MultiImage* sourceImage = NULL;
		MultiImage* targetImage = NULL;
		if (curr != first) {
			if (curr->prev->hasOutputCache && curr->prev->targetImageView == curr->sourceImageView)
				sourceImage = curr->prev->outputCache;
			else {
				sourceImage = worker;
				worker->setDimensions (curr->sourceImageView.getPixelWidth(), curr->sourceImageView.getPixelHeight());
				if (curr->prev->hasOutputCache) // There is a cache, but image views do not fit. Assume that sourceImageView is the part of prev->targetImageView.
					worker->copyFrom (curr->prev->outputCache, curr->sourceImageView.x - curr->prev->targetImageView.x, curr->sourceImageView.y - curr->prev->targetImageView.y, curr->sourceImageView.skip / curr->prev->targetImageView.skip);
			}
			sourceImage->convertTo (curr->inputColorSpace);
		}
		if (curr->hasOutputCache)
			targetImage = curr->outputCache;
		else {
			worker->setDimensions (curr->targetImageView.getPixelWidth(), curr->targetImageView.getPixelHeight());
			targetImage = worker;
		}
		curr->process (sourceImage, targetImage, buffer);
	}
}

void FilterChain::getReqiredBufferSize (int& w, int& h) {

	w = 0;
	h = 0;
	for (Filter* curr = first; curr; curr = curr->next) {
		int wc, hc;
		curr->getReqiredBufferSize (wc, hc);
		if (wc > w)
			w = wc;
		if (hc > h)
			h = hc;
	}
}

void FilterChain::getFullImageSize (int& w, int& h) {

	if (!last)
		return;
	else
		last->getFullImageSize (w, h);
}
}

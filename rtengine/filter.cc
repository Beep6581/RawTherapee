#include "filter.h"

namespace rtengine {

FilterDescriptor::FilterDescriptor (const std::string name, MultiImage::ColorSpace ics, MultiImage::ColorSpace ocs, bool forceCache)
	: name(name), inputColorSpace(ics), outputColorSpace(ocs), forceOutputCache(forceCache), plistener(NULL) {
}

void FilterDescriptor::addTriggerEvent (ProcEvent ev) {

	myEvents.insert (ev);
}

bool FilterDescriptor::myTriggerEvent (ProcEvent ev) const {

	 return myEvents.count (ev);
}

Filter::Filter (FilterDescriptor* descr)
	: descriptor(descr), next(NULL), prev(NULL), parent(NULL), cache(NULL),
	  forceOutputCache(descr->forceOutputCache()) {
}

Filter::~Filter () {

	delete cache;
}

void Filter::addNext (Filter* f) {

    next = f;
    f->prev = this;
}


void Filter::getReqiredBufferSize (ProcEventint& w, int& h) {

	// default implementation: filter does not require any buffer
	w = 0;
	h = 0;
}

void Filter::reverseTransPoint (int x, int y, int& xv, int& yv) {

	// default implementation: filter does not change image geometry
	xv = x;
	yv = y;
}

ImageView Filter::calculateTargetImageView (const ImageView& requestedImView) {

    return requestedImView;
}

ImageView Filter::calculateSourceImageView (const ImageView& requestedImView) {

    return requestedImView;
}

void Filter::getFullImageSize (int& w, int& h) {

	// default implementation: filter does not change image size
	if (prev)
		prev->getFullImageSize (w, h);
	else {
		w = -1; h = -1;
	}
}

bool Filter::isTriggerEvent (const set<ProcEvent>& events) {

	// check if we remain valid with 'events' occurred
	if (valid)
		for (set<ProcEvent>::iterator i=events.begin(); i!=events.end(); i++)
			if (descriptor->myTriggerEvent(*i))
				return true;
	return false;
}

double Filter::getScale () {

    if (parent)
        return parent->getScale ();
    else
        return 1.0 / targetImageView.skip;     // to be overridden in thumbnail image source to incorporate thumbnail size
}

void Filter::setupCache () {

	if (outputCache && !hasOutputCache) {
		delete outputCache;
		outputCache = NULL;
		return;
	}

	if (outputCache && outputCache->width==targetImageView.w && outputCache->height==targetImageView.h)
		return;

	delete outputCache;
	outputCache = new MultiImage (targetImageView.w, targetImageView.h, outputColorSpace);
}

void Filter::setProcParams (const ProcParams* pparams) {

	procParams = pparams;
}

}

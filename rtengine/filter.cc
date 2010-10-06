#include "filter.h"
#include <math.h>
#include "filterchain.h"

namespace rtengine {

FilterDescriptor::FilterDescriptor (const std::string name, MultiImage::ColorSpace ics, MultiImage::ColorSpace ocs, bool forceCache)
	: name(name), inputColorSpace(ics), outputColorSpace(ocs), forceOutCache(forceCache),
	  applyOnRawImage(true), applyOnStdImage(true), applyOnThumbnail(true) {
}

void FilterDescriptor::addTriggerEvent (ProcEvent ev) {

	myEvents.insert (ev);
}

bool FilterDescriptor::myTriggerEvent (ProcEvent ev) const {

	 return myEvents.count (ev);
}

Filter::Filter (FilterDescriptor* descr)
	: descriptor(descr), next(NULL), prev(NULL), shortCutPrev(NULL), parent(NULL),
	  outputCache(NULL), forceOutputCache(descr->forceOutputCache()), valid(false) {
}

Filter::~Filter () {

	delete outputCache;
}

void Filter::addNext (Filter* f) {

    next = f;
    f->prev = this;
}


Dim Filter::getReqiredBufferSize () {

	// default implementation: filter does not require any buffer
    return Dim ();
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

Dim Filter::getFullImageSize () {

	// default implementation: filter does not change image size
	if (prev)
		return prev->getFullImageSize ();
	else
		return Dim (-1, -1);
}

bool Filter::isTriggerEvent (const std::set<ProcEvent>& events) {

    for (std::set<ProcEvent>::iterator i=events.begin(); i!=events.end(); i++)
        if (*i<0 || descriptor->myTriggerEvent(*i))
            return true;
	return false;
}

double Filter::getScale () {

    if (parent)
        return parent->getScale () * parent->targetImageView.skip / sourceImageView.skip;
    else
        return getFilterChain()->getImageSource()->getScale() / targetImageView.skip;
}

void Filter::setupCache () {

	if (outputCache && !hasOutputCache) {
		delete outputCache;
		outputCache = NULL;
		return;
	}

	if (outputCache && outputCache->width==scaledTargetImageView.w && outputCache->height==scaledTargetImageView.h)
		return;

	delete outputCache;
	outputCache = new MultiImage (scaledTargetImageView.w, scaledTargetImageView.h, descriptor->getOutputColorSpace());
}

void Filter::setProcParams (ProcParams* pparams) {

	procParams = pparams;
}

}

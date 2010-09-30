#include "filterchaingroup.h"

namespace rtengine {

FilterChainGroup::FilterChainGroup (ImageSource* imgSource, ProcParams* pparams, bool multiThread)
	: imgSource(imgSource), procParams(pparams), buffer(NULL), workerImage(NULL), multiThread(multiThread) {
}

FilterChainGroup::~FilterChainGroup () {

	delete workerImage;
	delete buffer;
	for (int i=0; i<filterChains.size(); i++)
		delete filterChains[i];
}

void FilterChainGroup::addNewFilterChain (ImProcListener* listener) {

	if (filterChains.empty())
		filterChains.push_back (new FilterChain (listener, imgSource, procParams, multiThread));
	else
		filterChains.push_back (new FilterChain (listener, filterChains.back()));
	filterChains.back()->invalidate ();
}

void FilterChainGroup::removeFilterChain (ImProcListener* listener) {

	for (int i = 0; i<filterChains.size(); i++)
		if (filterChains[i]->getListener() == listener) {
		    if (i<filterChains.size()-1)
	            // remove from the chain
		        filterChains[i+1]->setNextChain (i==0 ? NULL : filterChains[i-1]);
			filterChains.erase (filterChains.begin() + i);
			break;
		}
}

double FilterChainGroup::getScale (ImProcListener* listener, int skip) {

    for (int i = 0; i<filterChains.size(); i++)
        if (filterChains[i]->getListener() == listener)
            return filterChains[i]->getScale (skip);
    return 1.0;
}

void FilterChainGroup::update (ImProcListener* listener) {

	for (int i=0; i<filterChains.size(); i++)
		if (!listener || filterChains[i]->getListener() == listener)
			filterChains[i]->invalidate ();

	std::set<ProcEvent> ev;
	ev.insert (EvAll);
	if (!listener)
	    // process all the filter chains
	    process (ev);
	else {
	    // process only the requested filter chain
	    FilterChain* fChain = NULL;
	    for (int i=0; i<filterChains.size(); i++)
	        if (filterChains[i]->getListener() == listener) {
	            fChain = filterChains[i];
	            break;
	        }
	    if (!fChain)
	        return;
	    Dim fullSize = filterChains[0]->getFullImageSize ();
	    Dim maxWorkerSize;
        fChain->setupProcessing (ev, fullSize, maxWorkerSize, true);
        if (!workerImage || workerImage->getAllocWidth()<maxWorkerSize.width || workerImage->getAllocHeight()<maxWorkerSize.height) {
            delete workerImage;
            workerImage = NULL;
            if (maxWorkerSize.width > 0 && maxWorkerSize.height > 0)
                workerImage = new MultiImage (maxWorkerSize.width, maxWorkerSize.height);
        }
        Dim bufferSize = fChain->getReqiredBufferSize ();
        if (!buffer || bufferSize.width > buffer->width || bufferSize.height > buffer->height)
            updateBuffer (bufferSize);
        fChain->process (ev, buffer, workerImage);
	}
}

void FilterChainGroup::process (const std::set<ProcEvent>& events) {

	if (filterChains.size()==0)
		return;

	Dim fullSize = filterChains[0]->getFullImageSize ();	// calculate it here once, it must be the same for the group
    Dim maxWorkerSize;

	// set up filter chains
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->setupProcessing (events, fullSize, maxWorkerSize, true);

	// re-allocate worker image, if necessary
	if (!workerImage || workerImage->getAllocWidth()!=maxWorkerSize.width || workerImage->getAllocHeight()!=maxWorkerSize.height) {
		delete workerImage;
		workerImage = NULL;
        if (maxWorkerSize.width > 0 && maxWorkerSize.height > 0)
            workerImage = new MultiImage (maxWorkerSize.width, maxWorkerSize.height);
	}

	// calculate common buffer of required size
	Dim bufferSize;
	for (int i=0; i<filterChains.size(); i++)
		bufferSize.setMax (filterChains[i]->getReqiredBufferSize ());
	updateBuffer (bufferSize);

	// process all filter chains
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->process (events, buffer, workerImage);
}

void FilterChainGroup::updateBuffer (Dim size) {

	bool deleteNeeded = buffer && (!size.nonZero() || buffer->width!=size.width || buffer->height!=size.height);
	bool createNeeded = size.nonZero() && (buffer->width!=size.width || buffer->height!=size.height);

	if (deleteNeeded) {
	    delete buffer;
	    buffer = NULL;
	}

	if (createNeeded)
	    buffer = new Buffer<int> (size.width, size.height);
}

}

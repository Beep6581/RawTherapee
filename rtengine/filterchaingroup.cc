#include "filterchaingroup.h"
#include "improclistener.h"

namespace rtengine {

FilterChainGroup::FilterChainGroup (ImageSource* imgSource, ProcParams* pparams, bool multiThread)
	: imgSource(imgSource), procParams(pparams), buffer(NULL), worker(NULL), multiThread(multiThread) {
}

FilterChainGroup::~FilterChainGroup () {

	delete worker;
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

	    fChain->setupProcessing (ev, true);
        updateBuffer (fChain->getReqiredBufferSize ());
        updateWorker (fChain->getReqiredWorkerSize ());
        fChain->process (ev, buffer, worker);
        notifyListener (fChain);
	}
}

void FilterChainGroup::process (const std::set<ProcEvent>& events) {

	if (filterChains.size()==0)
		return;

	// set up filter chains
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->setupProcessing (events, true);

	// calculate common buffer of required size
	Dim bufferSize;
	for (int i=0; i<filterChains.size(); i++)
		bufferSize.setMax (filterChains[i]->getReqiredBufferSize ());
	updateBuffer (bufferSize);

	// calculate size of the worker image
    Dim workerSize;
    for (int i=0; i<filterChains.size(); i++)
        workerSize.setMax (filterChains[i]->getReqiredWorkerSize ());
    updateWorker (workerSize);

	// process all filter chains
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->process (events, buffer, worker);

	// notify listeners
    for (int i=0; i<filterChains.size(); i++)
        notifyListener (filterChains[i]);
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

void FilterChainGroup::updateWorker (Dim size) {

    bool deleteNeeded = worker && (!size.nonZero() || worker->width!=size.width || worker->height!=size.height);
    bool createNeeded = size.nonZero() && (worker->width!=size.width || worker->height!=size.height);

    if (deleteNeeded) {
        delete worker;
        worker = NULL;
    }

    if (createNeeded)
        worker = new MultiImage (size.width, size.height);
}

void FilterChainGroup::notifyListener (FilterChain* chain) {

    ImProcListener* iml = chain->getListener ();
    if (iml) {
        Image8* img = chain->getDisplayImage ();
        iml->imageReady (img, chain->getLastScale(), chain->getFullImageSize(), chain->getLastImageView(), *procParams);
        delete img;
    }
}

}

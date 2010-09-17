#include "filterchaingroup.h"

namespace rtengine {
using namespace procparams;
FilterChainGroup::FilterChainGroup (ImageSource* imgSource, ProcParams* pparams, bool multiThread)
	: imgSource(imgSource), procParams(params), buffer(NULL), bufferWidth(0), bufferHeight(0), workerImage(NULL), multiThread(multiThread) {
}

FilterChainGroup::~FilterChainGroup () {

	delete workerImage;
	delete buffer;
	for (int i=0; i<filterChains.size(); i++)
		delete filterChains[i];
}

void FilterChainGroup::addNewFilterChain (ImProcListener* listener) {

	if (filterChains.empty())
		filterChains.push_back (new FilterChain (listener, imgSource, pparams, multiThread));
	else
		filterChains.push_back (new FilterChain (listener, filterChains.back()));
	filterChains.back()->invalidate ();
}

void FilterChainGroup::removeFilterChain (ImProcListener* listener) {

	for (std::vector<FilterChain*>::iterator i = filterChains.begin(); i!=filterChains.end(); i++)
		if (*i->getListener() == listener) {
			filterChains.erase (i);
			break;
		}
}

void FilterChainGroup::update (ImProcListener* listener) {

	for (int i=0; i<filterChains.size(); i++)
		if (!listener || filterChains[i]->getListener() == listener)
			filterChains[i]->invalidate ();
	// TODO: process () with empty set!!
}

void FilterChainGroup::process (const set<ProcEvent>& events) {

	if (filterChains.size()==0)
		return;

	int w, h;
	filterChains[0]->getFullImageSize (w, h);	// calculate it here once, it must be the same for the group

	// set up filter chains
	int maxWorkerWidth = -1, maxWorkerHeight = -1;
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->setupProcessing (events, w, h, maxWorkerWidth, maxWorkerHeight);

	// re-allocate worker image, if necessary
	if (!workerImage || workerImage->getAllocWidth()!=maxWorkerWidth || workerImage->getAllocHeight()!=maxWorkerHeight) {
		delete workerImage;
		workerImage = NULL;
		if (maxWorkerWidth > 0 && maxWorkerHeight > 0)
			workerImage = new MultiImage (maxWorkerWidth, maxWorkerHeight);
	}

	// calculate common buffer of required size
	int bw = 0, bh = 0;
	for (int i=0; i<filterChains.size(); i++) {
		int wc, hc;
		filterChains[i]->getReqiredBufferSize (wc, hc);
		if (wc > bw)
			bw = wc;
		if (hc > bh)
			bh = hc;
	}
	updateBuffer (bw, bh);

	// process all filter chains
	for (int i=0; i<filterChains.size(); i++)
		filterChains[i]->process (buffer);
}

void FilterChainGroup::updateBuffer (int bw, int bh) {

	bool deleteNeeded = buffer && (bw==0 || bh==0 || buffer->width!=bw || buffer->height!=bh);
	bool createNeeded = bw>0 && bh>0 && (buffer->width!=bw || buffer->height!=bh);

	if (deleteNeeded) {
	    delete buffer;
	    buffer = NULL;
	}

	if (createNeeded)
	    buffer = new Buffer<int> (bw, bh);
}

}

/*
 * interactiveimprocimpl.cc
 *
 *  Created on: Aug 23, 2010
 *      Author: gabor
 */

#include "interactiveimprocimpl.h"

namespace rtengine {

PreviewListenerAdapter::PreviewListenerAdapter () : prevListener (NULL) {}

void PreviewListenerAdapter::setPreviewListener (PreviewImageListener* prevListener) {

	this->prevListener = prevListener;
}

ImageView PreviewListenerAdapter::getViewToProcess (int fullW, int fullH) {

	return ImageView (0, 0, fullW, fullH, 16);
}

void PreviewListenerAdapter::imageReady	(IImage8* img, double scale, int fullW, int fullH, ImageView view, ProcParams params) {

	if (prevListener)
	    prevListener->imageReady (img, scale, fullW, fullH, params);
}


InteractiveImageProcessor* InteractiveImageProcessor::create (InitialImage* initialImage, PreviewImageListener* prevListener) {

	return new InteractiveImProcImpl ((ImageSource*)initialImage, prevListener);
}

InteractiveImProcImpl::InteractiveImProcImpl (ImageSource* imageSource, PreviewImageListener* prevListener)
	: imageSource (imageSource), thread (NULL), updaterRunning (false), destroying (false) {

	prevListAdapter.setPreviewListener (prevListener);
	filterChainGroup = new FilterChainGroup (imageSource, &params);
	filterChainGroup->addNewFilterChain (&prevListAdapter);
}

InteractiveImProcImpl::~InteractiveImProcImpl () {

	destroying = true;
    updaterThreadStart.lock ();
    if (updaterRunning && thread)
        thread->join ();
    mProcessing.lock();
    mProcessing.unlock();

    delete filterChainGroup;

    imageSource->decreaseRef ();
    updaterThreadStart.unlock ();
}

InitialImage* InteractiveImProcImpl::getInitialImage () {

	return imageSource;
}

void InteractiveImProcImpl::getParams (ProcParams& dst) {

	dst = params;
}

ProcParams* InteractiveImProcImpl::getParamsForUpdate (ProcEvent change) {

    paramsUpdateMutex.lock ();
    changeSinceLast.insert (change);
    return &nextParams;
}

void InteractiveImProcImpl::paramsUpdateReady () {

    paramsUpdateMutex.unlock ();
    startProcessing ();
}
void InteractiveImProcImpl::stopProcessing () {

	updaterThreadStart.lock ();
    if (updaterRunning && thread) {
        changeSinceLast.clear ();
        thread->join ();
    }
    updaterThreadStart.unlock ();
}

void InteractiveImProcImpl::startProcessing () {

	#undef THREAD_PRIORITY_NORMAL

	if (!destroying) {
		updaterThreadStart.lock ();
		if (!updaterRunning) {
			thread = NULL;
			updaterRunning = true;
			updaterThreadStart.unlock ();
			thread = Glib::Thread::create(sigc::mem_fun(*this, &InteractiveImProcImpl::process), 0, false, true, Glib::THREAD_PRIORITY_NORMAL);
		}
		else
			updaterThreadStart.unlock ();
	}
}

void InteractiveImProcImpl::process () {

    if (progressListener)
    	progressListener->setBusyFlag (true);

    std::set<ProcEvent> events;
    paramsUpdateMutex.lock ();
    while (!changeSinceLast.empty()) {
        params = nextParams;
        events = changeSinceLast;
        changeSinceLast.clear ();
        paramsUpdateMutex.unlock ();
        mProcessing.lock();
        filterChainGroup->process (events);
        mProcessing.unlock();
        paramsUpdateMutex.lock ();
    }
    paramsUpdateMutex.unlock ();
    updaterRunning = false;

    if (progressListener)
    	progressListener->setBusyFlag (false);
}

void InteractiveImProcImpl::createView (ImProcListener* listener) {

    mProcessing.lock();
	filterChainGroup->addNewFilterChain (listener);
    mProcessing.unlock();
}

void InteractiveImProcImpl::removeView (ImProcListener* listener) {

    mProcessing.lock();
	filterChainGroup->removeFilterChain (listener);
    mProcessing.unlock();
}

void InteractiveImProcImpl::fullUpdate (ImProcListener* listener) {

    mProcessing.lock();
	filterChainGroup->update (listener);
    mProcessing.unlock();
}

ColorTemp InteractiveImProcImpl::getAutoWB () {

	return imageSource->getAutoWB ();

}

ColorTemp InteractiveImProcImpl::getCamWB () {

    return imageSource->getCamWB();
}

ColorTemp InteractiveImProcImpl::getSpotWB (int x, int y, int rectSize) {

	// TODO
}

void InteractiveImProcImpl::getAutoCrop (double ratio, int &x, int &y, int &w, int &h) {

	// TODO
}

void InteractiveImProcImpl::saveInputICCReference (const Glib::ustring& fname) {

	// TODO
}

void InteractiveImProcImpl::setProgressListener (ProgressListener* l) {

	progressListener = l;
}
}

/*
 * interactiveimprocimpl.cc
 *
 *  Created on: Aug 23, 2010
 *      Author: gabor
 */

#include "interactiveimprocimpl.h"

#ifdef QTBUILD
	#define wakeUpThread	 hasJob.wakeAll();
	#define waitForJobFinish imProcThread->jobFinished.wait(&nextJobMutex);
	#define waitForJob 		 parent->hasJob.wait(&parent->nextJobMutex);

#else	
	#define wakeUpThread	 hasJob.signal ();
	#define waitForJobFinish imProcThread->jobFinished.wait(nextJobMutex);
	#define waitForJob 		 parent->hasJob.wait (parent->nextJobMutex);
#endif

namespace rtengine {

PreviewListenerAdapter::PreviewListenerAdapter () : prevListener (NULL) {}

void PreviewListenerAdapter::setPreviewListener (PreviewImageListener* prevListener) {

	this->prevListener = prevListener;
}

ImageView PreviewListenerAdapter::getViewToProcess (Dim fullSize) {

	return ImageView (0, 0, fullSize.width, fullSize.height, 16);
}

void PreviewListenerAdapter::imageReady	(const DisplayImage& img, double scale, Dim fullSize, ImageView view, ProcParams params) {

	if (prevListener)
	    prevListener->imageReady (img, scale, fullSize, params);
}

InteractiveImageProcessor* InteractiveImageProcessor::create (InitialImage* initialImage, PreviewImageListener* prevListener) {

	return new InteractiveImProcImpl ((ImageSource*)initialImage, prevListener);
}

InteractiveImProcImpl::InteractiveImProcImpl (ImageSource* imageSource, PreviewImageListener* prevListener)
	: imageSource (imageSource), destroyThread (false) {

	prevListAdapter.setPreviewListener (prevListener);
	filterChainGroup = new FilterChainGroup (imageSource, &params);
	filterChainGroup->addNewFilterChain (&prevListAdapter);
	
	imProcThread = new ImProcThread (this);
#ifdef QTBUILD
	imProcThread->start ();
#else
	#undef THREAD_PRIORITY_NORMAL
	Glib::Thread* thread = Glib::Thread::create(sigc::mem_fun(*imProcThread, &ImProcThread::run), 0, false, true, Glib::THREAD_PRIORITY_NORMAL);
#endif
}

InteractiveImProcImpl::~InteractiveImProcImpl () {

    nextJobMutex.lock ();
	changeSinceLast.clear ();
	destroyThread = true;
    wakeUpThread;
    waitForJobFinish;
    nextJobMutex.unlock ();

#ifdef QTBUILD
	imProcThread->wait ();
	delete imProcThread;
#endif

    delete filterChainGroup;
    imageSource->decreaseRef ();
}

InitialImage* InteractiveImProcImpl::getInitialImage () {

	return imageSource;
}

void InteractiveImProcImpl::getParams (ProcParams& dst) {

	dst = params;
}

ProcParams* InteractiveImProcImpl::getParamsForUpdate (ProcEvent change) {

    nextJobMutex.lock ();
    changeSinceLast.insert (change);
    return &nextParams;
}

void InteractiveImProcImpl::paramsUpdateReady () {

    wakeUpThread;
    nextJobMutex.unlock ();      
}
void InteractiveImProcImpl::stopProcessing () {

    nextJobMutex.lock ();
	changeSinceLast.clear ();
    wakeUpThread;
    waitForJobFinish;
    nextJobMutex.unlock ();
}

void ImProcThread::run () {

	bool killmyself = false;

	while (!killmyself) {

		parent->nextJobMutex.lock ();
		if (parent->changeSinceLast.empty ()) {
			if (parent->progressListener)
				parent->progressListener->setBusyFlag (false);
			waitForJob;
			if (parent->progressListener)
				parent->progressListener->setBusyFlag (true);
		}
		std::set<ProcEvent> events = parent->changeSinceLast;
		parent->changeSinceLast.clear ();
		parent->params = parent->nextParams;
		bool killmyself = parent->destroyThread;
		parent->nextJobMutex.unlock ();
		
		if (!killmyself) {
			if (!events.empty()) {
				parent->mProcessing.lock();
				parent->filterChainGroup->process (events);
				parent->mProcessing.unlock();
			}
			#ifdef QTBUILD
			jobFinished.wakeAll ();
			#else
			jobFinished.signal ();
			#endif
		}
	}
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

double InteractiveImProcImpl::getScale (ImProcListener* listener, int skip) {

    return  filterChainGroup->getScale (listener, skip);
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

void InteractiveImProcImpl::saveInputICCReference (const String& fname) {

	// TODO
}

void InteractiveImProcImpl::setProgressListener (ProgressListener* l) {

	progressListener = l;
}
}

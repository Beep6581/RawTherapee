/*
 * interactiveimprocimpl.h
 *
 *  Created on: Aug 23, 2010
 *      Author: gabor
 */

#ifndef INTERACTIVEIMPROCIMPL_H_
#define INTERACTIVEIMPROCIMPL_H_

#include "improclistener.h"
#include "rtengine.h"
#include "colortemp.h"
#include "filterchaingroup.h"
#include <set>

#ifdef QTBUILD
#include <QThread>
#endif

namespace rtengine {

class PreviewListenerAdapter : public ImProcListener {

        PreviewImageListener* prevListener;

	public:

		PreviewListenerAdapter ();
		void 		setPreviewListener (PreviewImageListener* prevListener);

		// implements ImProcListener interface:
		ImageView 	getViewToProcess 	(Dim fullSize);
		void		imageReady 			(const DisplayImage& img, double scale, Dim fullSize, ImageView view, ProcParams params);
};

class ThreadJob {
	
	public:
		ProcParams params;
		
		bool destroy;
};

class InteractiveImProcImpl;
#ifdef QTBUILD
class ImProcThread : public QThread {
#else
class ImProcThread {
#endif	
		InteractiveImProcImpl* parent;

	public:
		Condition jobFinished;
		ImProcThread (InteractiveImProcImpl* p) : parent(p) {}
		
		void run ();
};

class InteractiveImProcImpl : public InteractiveImageProcessor {

		friend class ImProcThread;
		
		ImageSource* imageSource;
		FilterChainGroup* filterChainGroup;
		ProcParams params;
		Mutex mProcessing;
		ProgressListener* progressListener;
		PreviewListenerAdapter prevListAdapter;

		// members of the updater:
		ImProcThread* imProcThread;
		std::set<ProcEvent> changeSinceLast;
		ProcParams nextParams;
		Condition hasJob;
		bool destroyThread;
		Mutex nextJobMutex;
		
	public:
		InteractiveImProcImpl  (ImageSource* imageSource, PreviewImageListener* prevListener);
		~InteractiveImProcImpl ();

		InitialImage* 			getInitialImage ();
		void          			getParams (ProcParams& dst);

		ProcParams*             getParamsForUpdate (ProcEvent change);
		void        			paramsUpdateReady ();
		void        			stopProcessing ();

		void		createView  (ImProcListener* listener);
		void		removeView  (ImProcListener* listener);
		void        fullUpdate  (ImProcListener* listener);
        double      getScale    (ImProcListener* listener, int skip);

		ColorTemp   getAutoWB   ();
		ColorTemp   getCamWB    ();
		ColorTemp   getSpotWB   (int x, int y, int rectSize);
		void        getAutoCrop (double ratio, int &x, int &y, int &w, int &h);

		void        saveInputICCReference (const String& fname);

		void        setProgressListener     (ProgressListener* l);
};

}

#endif



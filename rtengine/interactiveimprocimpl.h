/*
 * interactiveimprocimpl.h
 *
 *  Created on: Aug 23, 2010
 *      Author: gabor
 */

#ifndef INTERACTIVEIMPROCIMPL_H_
#define INTERACTIVEIMPROCIMPL_H_

#include "improclistener.h"
#include <glibmm.h>
#include "rtengine.h"
#include "colortemp.h"
#include "filterchaingroup.h"
#include <set>

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

class InteractiveImProcImpl : public InteractiveImageProcessor {

		ImageSource* imageSource;
		FilterChainGroup* filterChainGroup;
		Glib::Mutex mProcessing;
		ProcParams params;
		ProgressListener* progressListener;
		PreviewListenerAdapter prevListAdapter;

		// members of the updater:
		Glib::Thread* thread;
		Glib::Mutex updaterThreadStart;
		Glib::Mutex paramsUpdateMutex;
		bool updaterRunning;
		ProcParams nextParams;
		bool destroying;
		std::set<ProcEvent> changeSinceLast;

		void startProcessing ();
		void process ();

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



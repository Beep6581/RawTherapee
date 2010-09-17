/*
 * interactiveimprocimpl.h
 *
 *  Created on: Aug 23, 2010
 *      Author: gabor
 */

#ifndef INTERACTIVEIMPROCIMPL_H_
#define INTERACTIVEIMPROCIMPL_H_

#include <glibmm.h>
#include "rtengine.h"
#include "colortemp.h"
#include "filterchaingroup.h"
#include <set>

namespace rtengine {

class PreviewListenerAdapter : public ImProcListener {

		PreviewListener* prevListener;

	public:

		PreviewListenerAdapter ();
		void 		setPreviewListener (PreviewListener* prevListener);

		// implements ImProcListener interface:
		ImageView 	getViewToProcess 	(int fullW, int fullH);
		void		imageReady 			(IImage8* img, int fullW, int fullH, ImageView view, ProcParams params);
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
		InteractiveImProcImpl  (ImageSource* imageSource);
		~InteractiveImProcImpl ();

		InitialImage* 			getInitialImage ();
		void          			getParams (ProcParams& dst);

		ProcParams*             getParamsForUpdate (ProcEvent change) =0;
		void        			paramsUpdateReady () =0;
		void        			stopProcessing () =0;

		void		createView  (ImProcListener* listener) =0;
		void		removeView  (ImProcListener* listener) =0;
		void        fullUpdate  (ImProcListener* listener) =0;

		ColorTemp   getAutoWB   ();
		ColorTemp   getCamWB    ();
		ColorTemp   getSpotWB   (int x, int y, int rectSize);
		void        getAutoCrop (double ratio, int &x, int &y, int &w, int &h);

		void        saveInputICCReference (const Glib::ustring& fname);

		void        setProgressListener     (ProgressListener* l);
};

}


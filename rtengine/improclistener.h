#ifndef _IMPROCLISTENER_H_
#define _IMPROCLISTENER_H_

#include "rtcommon.h"
#include "imageview.h"
#include "procparams.h"

namespace rtengine {

class ImProcListener {

	public:
		virtual ImageView 	getViewToProcess 	(Dim fullSize) = 0;
		virtual void		imageReady 			(const DisplayImage& img, double scale, Dim fullSize, ImageView view, ProcParams params) = 0;
};

class PreviewImageListener {

	public:
		virtual void imageReady	(const DisplayImage& img, double scale, Dim fullSize, ProcParams params) = 0;
};

}

#endif

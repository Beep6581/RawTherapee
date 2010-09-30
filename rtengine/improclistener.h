#ifndef _IMPROCLISTENER_H_
#define _IMPROCLISTENER_H_

#include "rtengine.h"
#include "imageview.h"

namespace rtengine {

class ImProcListener {

	public:
		virtual ImageView 	getViewToProcess 	(Dim fullSize) = 0;
		virtual void		imageReady 			(IImage8* img, double scale, Dim fullSize, ImageView view, ProcParams params) = 0;
};

class PreviewImageListener {

	public:
		virtual void imageReady	(IImage8* img, double scale, Dim fullSize, ProcParams params) = 0;
};

}

#endif

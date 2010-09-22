#ifndef _IMPROCLISTENER_H_
#define _IMPROCLISTENER_H_

#include "rtengine.h"
#include "imageview.h"

namespace rtengine {

class ImProcListener {

	public:
		virtual ImageView 	getViewToProcess 	(int fullW, int fullH) = 0;
		virtual void		imageReady 			(IImage8* img, int fullW, int fullH, ImageView view, ProcParams params) = 0;
};

class PreviewImageListener {

	public:
		virtual void imageReady	(IImage8* img, int fullW, int fullH, int skip, ProcParams params) = 0;
};

}

#endif

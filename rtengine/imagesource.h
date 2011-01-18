/*
 * imagesource.h
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#ifndef IMAGESOURCE_H_
#define IMAGESOURCE_H_

#include "rtengine.h"
#include "colortemp.h"
#include "matrix33.h"
#include "coord2d.h"
#include "multiimage.h"

namespace rtengine {

class ImageSource : public InitialImage {

		int references;

	public:
		ImageSource ();

        virtual int load (const Glib::ustring& fileName, ProgressListener* listener = NULL) =0;

        virtual ColorTemp   getCamWB    () =0;
        virtual ColorTemp   getAutoWB   () =0;
        virtual ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) =0;
        virtual void        getAEHistogram (unsigned int* histogram, int& histcompr) =0;
        virtual Matrix33	getCamToRGBMatrix ()=0;
        virtual Matrix33	getRGBToCamMatrix ()=0;

        ImageSource* getImageSource () { return this; }

        virtual bool isRaw () =0;
        virtual bool isThumbnail () =0;
        virtual double getScale () =0;

        // inherited from InitialImage
        void increaseRef ();
		void decreaseRef ();

		virtual Dim  getFullImageSize () =0;
        virtual void getImage (const ImageView& targetImageView, MultiImage* targetImage) =0;
};

}

#endif /* IMAGESOURCE_H_ */

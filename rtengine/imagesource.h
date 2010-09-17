/*
 * imagesource.h
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#ifndef IMAGESOURCE_H_
#define IMAGESOURCE_H_

#include "rtengine.h"
#include "filter.h"
#include "colortemp.h"
#include "matrix33.h"

namespace rtengine {

class ImageSource : public InitialImage, public Filter {

		int references;

	public:
		ImageSource (FilterDescriptor* descr);

        virtual int load (const Glib::ustring& fileName, ProgressListener* listener = NULL);

        virtual ColorTemp   getCamWB    () =0;
        virtual ColorTemp   getAutoWB   () =0;
        virtual ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) =0;
        virtual double      getDefGain  () =0;
        virtual void        getAEHistogram (unsigned int* histogram, int& histcompr) =0;
        virtual Matrix33	getCamToRGBMatrix ()=0;
        virtual Matrix33	getRGBToCamMatrix ()=0;

        virtual bool isRaw () =0;

        // inherited from InitialImage
        void increaseRef ();
		void decreaseRef ();
};

}

#endif /* IMAGESOURCE_H_ */

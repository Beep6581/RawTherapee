/*
 * rawimagesource.h
 *
 *  Created on: Aug 22, 2010
 *      Author: gabor
 */

#ifndef RAWIMAGESOURCE_H_
#define RAWIMAGESOURCE_H_

#include "imagesource.h"
#include "rawimage.h"
#include "imagedata.h"
#include "colortemp.h"

namespace rtengine {

class RawImageSource : public ImageSource {

    private:
	 	int border;
        Glib::ustring fileName;
        RawImage* img;
        ImageData* idata;
        ColorTemp autoWB;
        bool autoWBComputed;
        cmsHPROFILE embProfile;

	public:

        RawImageSource ();
        virtual ~RawImageSource ();

		// inherited from InitialImage
		Glib::ustring 		 getFileName () { return fileName; }
		cmsHPROFILE 		 getEmbeddedProfile ();
		const ImageMetaData* getMetaData () { return idata; }

      // inherited from ImageSource
        int 		load (const Glib::ustring& fileName, ProgressListener* listener = NULL);

        ColorTemp   getCamWB    ();
        ColorTemp   getAutoWB   ();
        ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue);
        double      getDefGain  ();
        void        getAEHistogram (unsigned int* histogram, int& histcompr);
        Matrix33	getCamToRGBMatrix ();
        Matrix33	getRGBToCamMatrix ();

        // inherited from Filter
    	void getFullImageSize (int& w, int& h);
    	void process (const std::set<ProcEvent>& events, MultiImage* sourceImage, MultiImage* targetImage, Buffer<int>* buffer);

        bool isRaw () { return false; }
};

}

#endif /* RAWIMAGESOURCE_H_ */

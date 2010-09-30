/*
 * stdimagesource.h
 *
 *  Created on: Aug 19, 2010
 *      Author: gabor
 */

#ifndef STDIMAGESOURCE_H_
#define STDIMAGESOURCE_H_

#include "imagesource.h"
#include "image16.h"
#include "imagedata.h"
#include "colortemp.h"
#include <glibmm.h>
#include "imageview.h"
#include "multiimage.h"

namespace rtengine {

class StdImageSource : public ImageSource {

    private:
        Glib::ustring fileName;
        Image16* img;
        ImageData* idata;
        ColorTemp autoWB;
        bool autoWBComputed;


	public:

        StdImageSource ();
        virtual ~StdImageSource ();

		// inherited from InitialImage
		Glib::ustring 		 getFileName () { return fileName; }
		cmsHPROFILE 		 getEmbeddedProfile ();
		const ImageMetaData* getMetaData () { return idata; }

      // inherited from ImageSource
        int 		load (const Glib::ustring& fileName, ProgressListener* listener = NULL);

        ColorTemp   getCamWB    () =0;
        ColorTemp   getAutoWB   () =0;
        ColorTemp   getSpotWB   (std::vector<Coord2D> red, std::vector<Coord2D> green, std::vector<Coord2D>& blue) =0;
        double      getDefGain  () {return 1.0; }
        void        getAEHistogram (unsigned int* histogram, int& histcompr) =0;
        Matrix33	getCamToRGBMatrix ()=0;
        Matrix33	getRGBToCamMatrix ()=0;

    	Dim  getFullImageSize ();
    	void getImage (const ImageView& targetImageView, MultiImage* targetImage);

        bool isRaw ()       { return false; }
        bool isThumbnail () { return false; }
        double getScale ()  { return 1.0; }
};

}
#endif /* STDIMAGESOURCE_H_ */

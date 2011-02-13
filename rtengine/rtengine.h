/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _RTENGINE_
#define _RTENGINE_

#include "procparams.h"
#include "procevents.h"
#include <lcms2.h>
#include <string>
#include <time.h>
#include "settings.h"
#include "colortemp.h"
#include "imageview.h"
#include "dim.h"
#include "rtcommon.h"
#include "finalimage16.h"
#include <exiv2/exiv2.hpp>

/**
 * @file 
 * This file contains the main functionality of the raw therapee engine.
 *
 */

namespace rtengine {

  /**
    * This class represents provides functions to obtain exif and IPTC metadata information
    * from the image file
    */
    class ImageMetaData {

        public:
			virtual const Exiv2::ExifData& 	getExifData () const = 0;
			virtual const Exiv2::IptcData& 	getIptcData () const = 0;
			virtual const Exiv2::XmpData& 	getXmpData () const = 0;

          /** @return a struct containing the date and time of the image */
            virtual struct tm   getDateTime () const =0;
          /** @return the ISO of the image */
            virtual int         getISO () const =0;
          /** @return the F number of the image */
            virtual float       getFNumber  () const =0;
          /** @return the focal length used at the exposure */
            virtual float       getFocalLen () const =0;
          /** @return the shutter speed */
            virtual Exiv2::Rational getExposureTime () const =0;
          /** @return the maker of the camera */
            virtual String		getMake     () const =0;
          /** @return the model of the camera */
            virtual String      getModel    () const =0;
          /** @return the lens on the camera  */
            virtual String      getLens     () const =0;
            
            virtual int			getDefaultRotation () const=0;
            
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static String fNumberToString (float fNumber);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static String exposureTimeToString (Exiv2::Rational expTime);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static float fNumberFromString (const String& fNumber);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static Exiv2::Rational exposureTimeFromString (const String& expTime);
            
          /** Reads metadata from file.
            * @param fname is the name of the file
            * @param rml is a struct containing information about metadata location. Use it only for raw files. In case
            * of jpgs and tiffs pass a NULL pointer. 
            * @return The metadata */
            static ImageMetaData* fromFile (const String& fname);
    };

  /** This listener interface is used to indicate the progress of time consuming operations */
    class ProgressListener {
    
         public:
        /** This member function is called when the percentage of the progress has been changed.
          * @param p is a number between 0 and 1 */
          virtual void setProgress (double p) {}
        /** This member function is called when a textual information corresponding to the progress has been changed.
          * @param str is the textual information corresponding to the progress */
          virtual void setProgressStr (String str) {}
        /** This member function is called when the state of the processing has been changed.
          * @param busy =true if the processing has been started, =false if it has been finished */
          virtual void setBusyFlag (bool busy) {}
        /** This member function is called when an error occurs during the operation.
          * @param descr is the error message */
          virtual void error (String descr) {}
          virtual void startTimeConsumingOperation () {}
          virtual void progressReady () {}

    };

    class ImageSource;
   /**
    * This class represents an image loaded into the memory. It is the basis of further processing.
    * The embedded icc profile and metadata information can be obtained through this class, too.
    */
    class InitialImage {

        public:
          /** Returns the file name of the image.
            * @return The file name of the image */
            virtual String getFileName () =0;
          /** Returns the embedded icc profile of the image.
            * @return The handle of the embedded profile */
            virtual cmsHPROFILE getEmbeddedProfile () =0;
          /** Returns a class providing access to the exif and iptc metadata tags of the image.
            * @return An instance of the ImageMetaData class */
            virtual const ImageMetaData* getMetaData () =0;
          /** This class has manual reference counting. You have to call this function each time to make a new reference to an instance. */
            virtual void increaseRef () =0;
          /** This class has manual reference counting. You have to call this function each time to remove a reference 
            * (the last one deletes the instance automatically). */
            virtual void decreaseRef () =0;
            /** This is a function used for internal purposes only. */
            virtual ImageSource* getImageSource () =0;

          /** Loads an image into the memory. If it is a raw file, is is partially demosaiced (the time consuming part is done)
            * @param fname the name of the file
            * @param isRaw shall be true if it is a raw file
            * @param errorCode is a pointer to a variable that is set to nonzero if an error happened (output)
            * @param pl is a pointer pointing to an object implementing a progress listener. It can be NULL, in this case progress is not reported.
            * @return an object representing the loaded and pre-processed image */
            static InitialImage* load (const String& fname, bool isRaw, int& errorCode, ProgressListener* pl = NULL);
    };

    #include "improclistener.h"
    class PreviewImageListener;
    class ImProcListener;

    /** This is a staged, cached image processing manager with partial image update support.  */
    class InteractiveImageProcessor {

        public:
            /** Returns the inital image corresponding to the image processor.
              * @return the inital image corresponding to the image processor and increases its reference count */
            virtual InitialImage* getInitialImage () =0;
            /** Returns the current processing parameters.
              * @param dst is the location where the image processing parameters are copied (it is assumed that the memory is allocated by the caller) */
            virtual void          getParams (ProcParams& dst) =0;
            /** An essential member function. Call this when a setting has been changed. This function returns a pointer to the
              * processing parameters, that you have to update to reflect the changed situation. When ready, call the paramsUpdateReady 
              * function to start the image update.
              * @param change is the ID of the changed setting */
            virtual ProcParams* getParamsForUpdate (ProcEvent change) =0;
            /** An essential member function. This indicates that you are ready with the update of the processing parameters you got 
              * with the getParamsForUpdate call, so the image can be updated. This function returns immediately.
              * The image update starts immediately in the background. If it is ready, the result is passed to a PreviewImageListener
              * and to a DetailedCropListener (if enabled). */
            virtual void        paramsUpdateReady () =0;
            /** Stops image processing. When it returns, the image processing is already stopped. */
            virtual void        stopProcessing () =0;
            /** Creates a new view that acts as a window on the image */
            virtual void		createView  (ImProcListener* listener) =0;
            /** Removes the view corresponding to the listener */
            virtual void		removeView  (ImProcListener* listener) =0;
            /** Performs a full update on the image. If listener!=NULL, then only that image is updated that belongs to the given listener */
            virtual void        fullUpdate (ImProcListener* listener) =0;
            /** Returns the ratio of the dimensions of the obtained and the original image by using the given "skip" parameter */
            virtual double      getScale    (ImProcListener* listener, int skip) =0;

            virtual ColorTemp   getAutoWB   ();
            virtual ColorTemp   getCamWB    ();
            virtual ColorTemp   getSpotWB  (int x, int y, int rectSize) =0;
            virtual void        getAutoCrop (double ratio, int &x, int &y, int &w, int &h) =0;

            virtual void        saveInputICCReference (const String& fname) =0;

            virtual void        setProgressListener     (ProgressListener* l) =0;

            virtual ~InteractiveImageProcessor () {}

        /** Returns a staged, cached image processing manager supporting partial updates.
         *  Warning! the ownership of the initialImage passed to the function is taken over, so do not dereference it!
         * @param initialImage is a loaded and pre-processed initial image
         * @return the staged image processing manager */
            static InteractiveImageProcessor* create (InitialImage* initialImage, PreviewImageListener* prevListener);
    };

/** Checks if a raw file is supported
  * @param fname the name of the file
  * @param rml is a struct constaining informations on the location of the metadata in the raw file
  * @param rotation is the default angle of rotation (0, 90, 180, 270)
  * @param thumbWidth is the width of the embedded thumbnail file
  * @param thumbHeight is the height of the embedded thumbnail file
  * @param thumbOffset is the offset of the embedded thumbnail in the raw file
  * @param thumbType is the type of the embedded thumbnail (=0: no thumbnail, =1: jpeg format, =2: simple continuous image data in rgbrgb... order)
  * @return =0 if not supported */
    int getRawFileBasicInfo (const String& fname, int& rotation, int& thumbWidth, int& thumbHeight, int& thumbOffset, int& thumbType);

/** Returns the available output profile names
  * @return a vector of the available output profile names */
    StringList getOutputProfiles ();

/** Returns the available working profile names
  * @return a vector of the available working profile names */
    StringList getWorkingProfiles ();


    /** This class describes an image processing job. */
    class ProcessingJob {

        public:
    /** Creates a processing job with a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately. 
       * @param fname the name of the file
       * @param isRaw shall be true if it is a raw file
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */
        static ProcessingJob* create (const String& fname, bool isRaw, const ProcParams& pparams);
   
    /** Creates a processing job from a previously loaded file, thus an InitialImage instance.
       * This function always succeeds. It only stores the data into the ProcessingJob class, it does not do any
       * specific operation thus it returns immediately.
       * @param initialImage is a loaded and pre-processed initial image
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */   
        static ProcessingJob* create (InitialImage* initialImage, const ProcParams& pparams);

        virtual ~ProcessingJob () {}
    };

    /** This is a listener class that is notified when the SingleImageProcessor finished processing an image.
     * It supports batch processing by allowing the listener to return an other ProcessingJob instance.
     */
    class FinalImageListener {
        public:
         /** This function is called when processing of a ProcessingJob is ready.
          * It has to return with the next job, or with NULL if there are no jobs left.
          * @param img is the result of the last ProcessingJob
          * @return the next ProcessingJob to process */
            virtual ProcessingJob* imageReady (FinalImage16* img) =0;
    };

    /** Through this class it is possible to process a single image without support for interactive editing
     */
    class SingleImageProcessor {

        public:
            /** This function performs the processing of a single ProcessingJob. It has two modes: if the FinalImageListener
             * passed is NULL, then the call will be synchronous, thus it returns only when the processing is accomplished and
             * returns the result as a return value. Alternatively, if a non-NULL FinalImageListener is given, it starts the processing
             * in the background and terminates immediately by returning a NULL.
             * As soon as the job is ready, the given FinalImageListener is notified.
             * @param pJob is the ProcessingJob to process
             * @param pListener a progress listener to follow the progress of the processing of the current image. It can be NULL.
             * @param fiListener can be NULL if a synchronous call is required and non-NULL is a background processing is initiated.
             * @param errorCode is returned for example if the image specified by a file name can not be loaded or if it is not supported.
             */
            static FinalImage16* process (ProcessingJob* pJob, ProgressListener* pListener, FinalImageListener* fiListener, int& errorCode);
    };
}

#endif


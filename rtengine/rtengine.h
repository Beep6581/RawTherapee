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

#include "rawmetadatalocation.h"
#include "procparams.h"
#include "procevents.h"
#include <lcms2.h>
#include <string>
#include <glibmm.h>
#include <time.h>
#include <rtexif.h>
#include "iimage.h"
#include "settings.h"
#include "colortemp.h"
#include "imageview.h"
#include "dim.h"

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
          /** Checks the availability of exif metadata tags.
            * @return Returns true if image contains exif metadata tags */
            virtual bool hasExif () const =0;        
          /** Returns the directory of exif metadata tags.
            * @return The directory of exif metadata tags */
            virtual const rtexif::TagDirectory* getExifData () const =0;
          /** Checks the availability of IPTC tags.
            * @return Returns true if image contains IPTC tags */
            virtual bool hasIPTC () const =0;
          /** Returns the directory of IPTC tags.
            * @return The directory of IPTC tags */
            virtual const std::vector<IPTCPair> getIPTCData () const =0;
          /** @return a struct containing the date and time of the image */
            virtual struct tm   getDateTime () const =0;
          /** @return the ISO of the image */
            virtual int         getISOSpeed () const =0;
          /** @return the F number of the image */
            virtual double      getFNumber  () const =0;
          /** @return the focal length used at the exposure */
            virtual double      getFocalLen () const =0;
          /** @return the shutter speed */
            virtual double      getShutterSpeed () const =0;
          /** @return the maker of the camera */
            virtual std::string getMake     () const =0;
          /** @return the model of the camera */
            virtual std::string getModel    () const =0;
          /** @return the lens on the camera  */
            virtual std::string getLens     () const =0;
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static std::string apertureToString (double aperture);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static std::string shutterToString (double shutter);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static double apertureFromString (std::string shutter);
          /** Functions to convert between floating point and string representation of shutter and aperture */
            static double shutterFromString (std::string shutter);
            
          /** Reads metadata from file.
            * @param fname is the name of the file
            * @param rml is a struct containing information about metadata location. Use it only for raw files. In case
            * of jpgs and tiffs pass a NULL pointer. 
            * @return The metadata */
            static ImageMetaData* fromFile (const Glib::ustring& fname, RawMetaDataLocation* rml);
    };

  /** This listener interface is used to indicate the progress of time consuming operations */
    class ProgressListener {
    
         public:
        /** This member function is called when the percentage of the progress has been changed.
          * @param p is a number between 0 and 1 */
          virtual void setProgress (double p) {}
        /** This member function is called when a textual information corresponding to the progress has been changed.
          * @param str is the textual information corresponding to the progress */
          virtual void setProgressStr (Glib::ustring str) {}
        /** This member function is called when the state of the processing has been changed.
          * @param busy =true if the processing has been started, =false if it has been finished */
          virtual void setBusyFlag (bool busy) {}
        /** This member function is called when an error occurs during the operation.
          * @param descr is the error message */
          virtual void error (Glib::ustring descr) {}
          virtual void startTimeConsumingOperation () {}
          virtual void progressReady () {}

    };
    
   /**
    * This class represents an image loaded into the memory. It is the basis of further processing.
    * The embedded icc profile and metadata information can be obtained through this class, too.
    */
    class InitialImage {

        public:
          /** Returns the file name of the image.
            * @return The file name of the image */
            virtual Glib::ustring getFileName () =0;
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

          /** Loads an image into the memory. If it is a raw file, is is partially demosaiced (the time consuming part is done)
            * @param fname the name of the file
            * @param isRaw shall be true if it is a raw file
            * @param errorCode is a pointer to a variable that is set to nonzero if an error happened (output)
            * @param pl is a pointer pointing to an object implementing a progress listener. It can be NULL, in this case progress is not reported.
            * @return an object representing the loaded and pre-processed image */
            static InitialImage* load (const Glib::ustring& fname, bool isRaw, int* errorCode, ProgressListener* pl = NULL);
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

            virtual void        saveInputICCReference (const Glib::ustring& fname) =0;

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
    int getRawFileBasicInfo (const Glib::ustring& fname, RawMetaDataLocation& rml, int& rotation, int& thumbWidth, int& thumbHeight, int& thumbOffset, int& thumbType);

/** Returns the available output profile names
  * @return a vector of the available output profile names */
    std::vector<std::string> getOutputProfiles ();

/** Returns the available working profile names
  * @return a vector of the available working profile names */
    std::vector<std::string> getWorkingProfiles ();

    /** This class  holds all the necessary informations to accomplish the full processing of the image */
    class ProcessingJob {

        public:

    /** Creates a processing job from a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately. 
       * @param fname the name of the file
       * @param isRaw shall be true if it is a raw file
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */
        static ProcessingJob* create (const Glib::ustring& fname, bool isRaw, const ProcParams& pparams);
   
    /** Creates a processing job from a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately. This function increases the reference count of the initialImage. If you decide not the process the image you
       * have to cancel it by calling the member function void cancel(). If the image is processed the reference count of initialImage is decreased automatically, thus the ProcessingJob
       * instance gets invalid. You can not use a ProcessingJob instance to process an image twice.
       * @param initialImage is a loaded and pre-processed initial image
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */   
        static ProcessingJob* create (InitialImage* initialImage, const ProcParams& pparams);

    /** Cancels and destroys a processing job. The reference count of the corresponding initialImage (if any) is decreased. After the call of this function the ProcessingJob instance
      * gets invalid, you must not use it any more. Dont call this function while the job is being processed. 
      * @param job is the job to destroy */
        static void destroy (ProcessingJob* job);
    };

/** This function performs all the image processinf steps corresponding to the given ProcessingJob. It returns when it is ready, so it can be slow.
   * The ProcessingJob passed becomes invalid, you can not use it any more.
   * @param job the ProcessingJob to cancel. 
   * @param errorCode is the error code if an error occured (e.g. the input image could not be loaded etc.) 
   * @param pl is an optional ProgressListener if you want to keep track of the progress
   * @return the resulting image, with the output profile applied, exif and iptc data set. You have to save it or you can access the pixel data directly.  */  
    IImage16* processImage (ProcessingJob* job, int& errorCode, ProgressListener* pl = NULL);

/** This class is used to control the batch processing. The class implementing this interface will be called when the full processing of an
   * image is ready and the next job to process is needed. */
    class BatchProcessingListener : public ProgressListener {
        public:
         /** This function is called when an image gets ready during the batch processing. It has to return with the next job, or with NULL if
                        * there is no jobs left.
                        * @param img is the result of the last ProcessingJob 
                        * @return the next ProcessingJob to process */  
            virtual ProcessingJob* imageReady (IImage16* img) =0;
    };
/** This function performs all the image processing steps corresponding to the given ProcessingJob. It runs in the background, thus it returns immediately,
   * When it finishes, it calls the BatchProcessingListener with the resulting image and asks for the next job. It the listener gives a new job, it goes on 
   * with processing. If no new job is given, it finishes.
   * The ProcessingJob passed becomes invalid, you can not use it any more.
   * @param job the ProcessingJob to cancel. 
   * @param bpl is the BatchProcessingListener that is called when the image is ready or the next job is needed. It also acts as a ProgressListener. */  
    void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl);
}

#endif


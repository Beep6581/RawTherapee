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

#include <procparams.h>
#include <procevents.h>
#include <lcms.h>
#include <string>
#include <glibmm.h>
#include <time.h>
#include <rtexif.h>
#include <rawmetadatalocation.h>
#include <iimage.h>
#include <utils.h>
#include <settings.h>
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
            virtual bool hasExif () const {}        
          /** Returns the directory of exif metadata tags.
            * @return The directory of exif metadata tags */
            virtual const rtexif::TagDirectory* getExifData () const {}
          /** Checks the availability of IPTC tags.
            * @return Returns true if image contains IPTC tags */
            virtual bool hasIPTC () const {}
          /** Returns the directory of IPTC tags.
            * @return The directory of IPTC tags */
            virtual const std::vector<procparams::IPTCPair> getIPTCData () const {}
          /** @return a struct containing the date and time of the image */
            virtual struct tm   getDateTime () const {}
          /** @return the ISO of the image */
            virtual int         getISOSpeed () const {}
          /** @return the F number of the image */
            virtual double      getFNumber  () const {}
          /** @return the focal length used at the exposure */
            virtual double      getFocalLen () const {}
          /** @return the shutter speed */
            virtual double      getShutterSpeed () const {}
          /** @return the maker of the camera */
            virtual std::string getMake     () const {}
          /** @return the model of the camera */
            virtual std::string getModel    () const {}
          /** @return the lens on the camera  */
            virtual std::string getLens     () const {}
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
          * @param state =1 if the processing has been started, =0 if it has been stopped */
          virtual void setProgressState (int state) {}
        /** This member function is called when an error occurs during the operation.
          * @param descr is the error message */
          virtual void error (Glib::ustring descr) {}
    };
    
    class ImageSource;

  /**
    * This class represents an image loaded into the memory. It is the basis of further processing.
    * In case of raw files the most time consuming operation, the demosaicing is already performed.
    * The embedded icc profile and metadata information can be obtained through this class, too.
    */
    class InitialImage {

        public:
          /** Returns the file name of the image.
            * @return The file name of the image */
            virtual Glib::ustring getFileName () {}
          /** Returns the embedded icc profile of the image.
            * @return The handle of the embedded profile */
            virtual cmsHPROFILE getEmbeddedProfile () {}
          /** Returns a class providing access to the exif and iptc metadata tags of the image.
            * @return An instance of the ImageMetaData class */
            virtual const ImageMetaData* getMetaData () {}
          /** This is a function used for internal purposes only. */
            virtual ImageSource* getImageSource () {}
          /** This class has manual reference counting. You have to call this function each time to make a new reference to an instance. */
            virtual void increaseRef () {}
          /** This class has manual reference counting. You have to call this function each time to remove a reference 
            * (the last one deletes the instance automatically). */
            virtual void decreaseRef () {}


          /** Loads an image into the memory. If it is a raw file, is is partially demosaiced (the time consuming part is done)
            * @param fname the name of the file
            * @param isRaw shall be true if it is a raw file
            * @param errorCode is a pointer to a variable that is set to nonzero if an error happened (output)
            * @param pl is a pointer pointing to an object implementing a progress listener. It can be NULL, in this case progress is not reported.
            * @return an object representing the loaded and pre-processed image */
            static InitialImage* load (const Glib::ustring& fname, bool isRaw, int* errorCode, ProgressListener* pl = NULL);
    };

    /** When the preview image is ready for display during staged processing (thus the changes have been updated),
      * the staged processor notifies the listener class implementing a PreviewImageListener.
      * It is important to note that the file passed to the listener can be used in a shared manner (copying it is not
      * needed) as long as the mutex corresponding to the image is used every time the image is accessed.
      * If the scale of the preview image is >1, no sharpening, no denoising and no cropping is applied, and 
      * the transform operations (rotate, c/a, vignetting correction) are performed using a faster less quality algorithm. 
      * The image you get with this listener is created to display on the monitor (monitor profile has been already applied). */
    class PreviewImageListener {
        public: 
            /** With this member function the staged processor notifies the listener that it allocated a new
             * image to store the end result of the processing. It can be used in a shared manner. 
             * @param img is a pointer to the image
             * @param scale describes the current scaling applied compared to the 100% size (preview scale + resize scale) 
             * @param cp holds the coordinates of the current crop rectangle */
            virtual void setImage   (IImage8* img, double scale, procparams::CropParams cp) {}
            /** With this member function the staged processor notifies the listener that the image passed as parameter
              * will be deleted, and no longer used to store the preview image.
              * @param img the pointer to the image to be destroyed. The listener has to free the image!  */
            virtual void delImage   (IImage8* img) {}
            /** With this member function the staged processor notifies the listener that the preview image has been updated.
              * @param cp holds the coordinates of the current crop rectangle */
            virtual void imageReady (procparams::CropParams cp) {}
    };

    /** When the detailed crop image is ready for display during staged processing (thus the changes have been updated),
      * the staged processor notifies the listener class implementing a DetailedCropListener.
      * It is important to note that the file passed to the listener can <b>not</b> be used in a shared manner, the class
      * implementing this interface has to store a copy of it. */
    class DetailedCropListener {
        public: 
            /** With this member function the staged processor notifies the listener that the detailed crop image has been updated.
              * @param img is a pointer to the detailed crop image */
            virtual void setDetailedCrop (IImage8* img, procparams::CropParams cp, int cx, int cy, int cw, int ch, int skip) {}
            virtual bool getWindow       (int& cx, int& cy, int& cw, int& ch, int& skip) { return false; }
    };

    /** This listener is used when the full size of the final image has been changed (e.g. rotated by 90 deg.) */
    class SizeListener {
        public: 
            /** This member function is called when the size of the final image has been changed
              * @param w is the width of the final image (without cropping)           
              * @param h is the height of the final image (without cropping)
              * @param ow is the width of the final image (without resizing and cropping)           
              * @param oh is the height of the final image (without resizing and cropping) */
            virtual void sizeChanged (int w, int h, int ow, int oh) {}
    };

    /** This listener is used when the histogram of the final image has changed. */
    class HistogramListener {
        public:
            /** This member function is called when the histogram of the final image has changed.
              * @param redh is the array of size 256 containing the histogram of the red channel
              * @param greenh is the array of size 256 containing the histogram of the green channel
              * @param blueh is the array of size 256 containing the histogram of the blue channel
              * @param lumah is the array of size 256 containing the histogram of the luminance channel */
            virtual void histogramChanged (unsigned int* redh, unsigned int* greenh, unsigned int* blueh, unsigned int* lumah) {}
    };

    /** This listener is used when the auto exposure has been recomputed (e.g. when the clipping ratio changed). */
    class AutoExpListener {
        public:
            /** This member function is called when the auto exposure has been recomputed.
              * @param brightness is the new brightness value (in logarithmic scale)
              * @param black is the new black level (measured in absolute pixel data) */
            virtual void autoExpChanged (double brightness, int black) {}
    };

    /** This class represents a detailed part of the image (looking through a kind of window).
      * It can be created and destroyed with the appropriate members of StagedImageProcessor.
      * Several crops can be assigned to the same image.   */
    class DetailedCrop {
        public:
            /** Sets the window defining the crop. */
            virtual void setWindow   (int cx, int cy, int cw, int ch, int skip) {} 
            /** Perform a full recalculation of the part of the image corresponding to the crop. */
            virtual void fullUpdate  () {}
            /** Sets the listener of the crop. */
            virtual void setListener (DetailedCropListener* il) {}       
            /** Destroys the crop. */
            virtual void destroy () {}       
    };

    /** This is a staged, cached image processing manager with partial image update support.  */
    class StagedImageProcessor {

        public:
            /** Returns the inital image corresponding to the image processor.
              * @return the inital image corresponding to the image processor */
            virtual InitialImage* getInitialImage () {}
            /** Returns the current processing parameters.
              * @param dst is the location where the image processing parameters are copied (it is assumed that the memory is allocated by the caller) */
            virtual void        getParams (procparams::ProcParams* dst) {}
            /** An essential member function. Call this when a setting has been changed. This function returns a pointer to the
              * processing parameters, that you have to update to reflect the changed situation. When ready, call the paramsUpdateReady 
              * function to start the image update.
              * @param change is the ID of the changed setting */
            virtual procparams::ProcParams* getParamsForUpdate (ProcEvent change) {}
            /** An essential member function. This indicates that you are ready with the update of the processing parameters you got 
              * with the getParamsForUpdate call, so the image can be updated. This function returns immediately.
              * The image update starts immediately in the background. If it is ready, the result is passed to a PreviewImageListener
              * and to a DetailedCropListener (if enabled). */
            virtual void        paramsUpdateReady () {}
            /** Stops image processing. When it returns, the image processing is already stopped. */
            virtual void        stopProcessing () {}
            /** Sets the scale of the preview image. The larger the number is, the faster the image updates are (typical values are 4-5).
              * @param scale is the scale of the preview image */
            virtual void        setPreviewScale (int scale) {}
            /** Returns the scale of the preview image. 
              * @return the current scale of the preview image */
            virtual int         getPreviewScale () {}
            /** Performs a full update on the preview image. The resulting image is passed to the listener. */
            virtual void        fullUpdatePreviewImage () {}
            /** Performs a full update on the detailed crops corresponding to the image. The resulting images are passed to the listeners of the crops. */
            virtual void        fullUpdateDetailedCrops () {}
            /** Returns the full width of the resulting image (in 1:1 scale). 
              * @return the width of the final image */
            virtual int         getFullWidth () {}
            /** Returns the full height of the resulting image (in 1:1 scale). 
              * @return the height of the final image */
            virtual int         getFullHeight () {}
            /** Returns the width of the preview image. 
              * @return the width of the preview image */
            virtual int         getPreviewWidth () {}
            /** Returns the height of the preview image. 
              * @return the height of the preview image */
            virtual int         getPreviewHeight () {}

            /** Creates and returns a Crop instance that acts as a window on the image */
            virtual DetailedCrop* createCrop  () {}

            virtual void        getAutoWB   (double& temp, double& green) {}
            virtual void        getCamWB    (double& temp, double& green) {}
            virtual void        getSpotWB  (int x, int y, int rectSize, double& temp, double& green) {}
            virtual void        getAutoCrop (double ratio, int &x, int &y, int &w, int &h) {}

            virtual void        saveInputICCReference (const Glib::ustring& fname) {}

            virtual void        setProgressListener     (ProgressListener* l) {}
            virtual void        setSizeListener         (SizeListener* l) {}
            virtual void        delSizeListener         (SizeListener* l) {}
            virtual void        setAutoExpListener      (AutoExpListener* l) {}
            virtual void        setHistogramListener    (HistogramListener *l) {}
            virtual void        setPreviewImageListener (PreviewImageListener* l) {}

            virtual ~StagedImageProcessor () {}

        /** Returns a staged, cached image processing manager supporting partial updates
        * @param initialImage is a loaded and pre-processed initial image
        * @return the staged image processing manager */
            static StagedImageProcessor* create (InitialImage* initialImage);
            static void destroy (StagedImageProcessor* sip);
    };


/** 
  *  Initializes the RT engine
  * @param s is a struct of basic settings */
    int init (const Settings* s);

/** Checks if a raw file is supported
  * @param fname the name of the file
  * @param rml is a struct constaining informations on the location of the metadata in the raw file
  * @param rotation is the default angle of rotation (0, 90, 180, 270)
  * @param thumbWidth is the width of the embedded thumbnail file
  * @param thumbHeight is the height of the embedded thumbnail file
  * @param thumbOffset is the offset of the embedded thumbnail in the raw file
  * @param thumbType is the type of the embedded thumbnail (=0: no thumbnail, =1: jpeg format, =2: simple continuous image data in rgbrgb... order)
  * @return =0 if not supported */
    int getRawFileBasicInfo (const Glib::ustring& fname, rtengine::RawMetaDataLocation& rml, int& rotation, int& thumbWidth, int& thumbHeight, int& thumbOffset, int& thumbType);

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
        static ProcessingJob* create (const Glib::ustring& fname, bool isRaw, const procparams::ProcParams& pparams);
   
    /** Creates a processing job from a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately. This function increases the reference count of the initialImage. If you decide not the process the image you
       * have to cancel it by calling the member function void cancel(). If the image is processed the reference count of initialImage is decreased automatically, thus the ProcessingJob
       * instance gets invalid. You can not use a ProcessingJob instance to process an image twice.
       * @param initialImage is a loaded and pre-processed initial image
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */   
        static ProcessingJob* create (InitialImage* initialImage, const procparams::ProcParams& pparams);

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
            virtual ProcessingJob* imageReady (IImage16* img) {}
    };
/** This function performs all the image processinf steps corresponding to the given ProcessingJob. It runs in the background, thus it returns immediately,
   * When it finishes, it calls the BatchProcessingListener with the resulting image and asks for the next job. It the listener gives a new job, it goes on 
   * with processing. If no new job is given, it finishes.
   * The ProcessingJob passed becomes invalid, you can not use it any more.
   * @param job the ProcessingJob to cancel. 
   * @param bpl is the BatchProcessingListener that is called when the image is ready or the next job is needed. It also acts as a ProgressListener. */  
    void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl);

    
    extern Glib::Mutex* lcmsMutex;
}

#endif


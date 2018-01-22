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

#include "imageformat.h"
#include "rt_math.h"
#include "procparams.h"
#include "procevents.h"
#include <lcms2.h>
#include <string>
#include <glibmm.h>
#include <ctime>
#include "../rtexif/rtexif.h"
#include "rawmetadatalocation.h"
#include "iimage.h"
#include "utils.h"
#include "../rtgui/threadutils.h"
#include "settings.h"
#include "LUT.h"

/**
 * @file
 * This file contains the main functionality of the RawTherapee engine.
 *
 */

class EditDataProvider;

namespace rtengine
{

class IImage8;
class IImage16;
class IImagefloat;
class ImageSource;

/**
  * This class provides functions to obtain exif and IPTC metadata information
  * from any of the sub-frame of an image file
  */
class FramesMetaData
{

public:
    /** @return Returns the number of root Metadata */
    virtual unsigned int getRootCount () const = 0;
    /** @return Returns the number of frame contained in the file based on Metadata */
    virtual unsigned int getFrameCount () const = 0;

    /** Checks the availability of exif metadata tags.
      * @return Returns true if image contains exif metadata tags */
    virtual bool hasExif (unsigned int frame = 0) const = 0;
    /** Returns the directory of exif metadata tags.
      * @param root root number in the metadata tree
      * @return The directory of exif metadata tags */
    virtual rtexif::TagDirectory* getRootExifData (unsigned int root = 0) const = 0;
    /** Returns the directory of exif metadata tags.
      * @param frame frame number in the metadata tree
      * @return The directory of exif metadata tags */
    virtual rtexif::TagDirectory* getFrameExifData (unsigned int frame = 0) const = 0;
    /** Returns the directory of exif metadata tags containing at least the 'Make' tag for the requested frame.
      * If no usable metadata exist in the frame, send back the best TagDirectory describing the frame content.
      * @param imgSource rawimage that we want the metadata from
      * @param rawParams RawParams to select the frame number
      * @return The directory of exif metadata tags containing at least the 'Make' tag */
    virtual rtexif::TagDirectory* getBestExifData (ImageSource *imgSource, procparams::RAWParams *rawParams) const = 0;
    /** Checks the availability of IPTC tags.
      * @return Returns true if image contains IPTC tags */
    virtual bool hasIPTC (unsigned int frame = 0) const = 0;
    /** Returns the directory of IPTC tags.
      * @return The directory of IPTC tags */
    virtual procparams::IPTCPairs getIPTCData (unsigned int frame = 0) const = 0;
    /** @return a struct containing the date and time of the image */
    virtual tm getDateTime (unsigned int frame = 0) const = 0;
    /** @return a timestamp containing the date and time of the image */
    virtual time_t getDateTimeAsTS(unsigned int frame = 0) const = 0;
    /** @return the ISO of the image */
    virtual int getISOSpeed (unsigned int frame = 0) const = 0;
    /** @return the F number of the image */
    virtual double getFNumber  (unsigned int frame = 0) const = 0;
    /** @return the focal length used at the exposure */
    virtual double getFocalLen (unsigned int frame = 0) const = 0;
    /** @return the focal length in 35mm used at the exposure */
    virtual double getFocalLen35mm (unsigned int frame = 0) const = 0;
    /** @return the focus distance in meters, 0=unknown, 10000=infinity */
    virtual float getFocusDist (unsigned int frame = 0) const = 0;
    /** @return the shutter speed */
    virtual double getShutterSpeed (unsigned int frame = 0) const = 0;
    /** @return the exposure compensation */
    virtual double getExpComp (unsigned int frame = 0) const = 0;
    /** @return the maker of the camera */
    virtual std::string getMake     (unsigned int frame = 0) const = 0;
    /** @return the model of the camera */
    virtual std::string getModel    (unsigned int frame = 0) const = 0;

    std::string getCamera   (unsigned int frame = 0) const
    {
        return getMake(frame) + " " + getModel(frame);
    }

    /** @return the lens on the camera  */
    virtual std::string getLens     (unsigned int frame = 0) const = 0;
    /** @return the orientation of the image */
    virtual std::string getOrientation (unsigned int frame = 0) const = 0;

    /** @return true if the file is a PixelShift shot (Pentax and Sony bodies) */
    virtual bool getPixelShift (unsigned int frame = 0) const = 0;
    /** @return false: not an HDR file ; true: single or multi-frame HDR file (e.g. Pentax HDR raw file or 32 bit float DNG file or Log compressed) */
    virtual bool getHDR (unsigned int frame = 0) const = 0;
    /** @return the sample format based on MetaData */
    virtual IIOSampleFormat getSampleFormat (unsigned int frame = 0) const = 0;

    /** Functions to convert between floating point and string representation of shutter and aperture */
    static std::string apertureToString (double aperture);
    /** Functions to convert between floating point and string representation of shutter and aperture */
    static std::string shutterToString (double shutter);
    /** Functions to convert between floating point and string representation of shutter and aperture */
    static double apertureFromString (std::string shutter);
    /** Functions to convert between floating point and string representation of shutter and aperture */
    static double shutterFromString (std::string shutter);
    /** Functions to convert between floating point and string representation of exposure compensation */
    static std::string expcompToString (double expcomp, bool maskZeroexpcomp);

    virtual ~FramesMetaData () = default;

    /** Reads metadata from file.
      * @param fname is the name of the file
      * @param rml is a struct containing information about metadata location of the first frame.
      * Use it only for raw files. In caseof jpgs and tiffs pass a NULL pointer.
      * @param firstFrameOnly must be true to get the MetaData of the first frame only, e.g. for a PixelShift file.
      * @return The metadata */
    static FramesMetaData* fromFile (const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml, bool firstFrameOnly = false);
};

/** This listener interface is used to indicate the progress of time consuming operations */
class ProgressListener
{

public:
    virtual ~ProgressListener() {}
    /** This member function is called when the percentage of the progress has been changed.
      * @param p is a number between 0 and 1 */
    virtual void setProgress (double p) {}
    /** This member function is called when a textual information corresponding to the progress has been changed.
      * @param str is the textual information corresponding to the progress */
    virtual void setProgressStr (Glib::ustring str) {}
    /** This member function is called when the state of the processing has been changed.
      * @param inProcessing =true if the processing has been started, =false if it has been stopped */
    virtual void setProgressState (bool inProcessing) {}
    /** This member function is called when an error occurs during the operation.
      * @param descr is the error message */
    virtual void error (Glib::ustring descr) {}
};

class ImageSource;

/**
  * This class represents an image loaded into the memory. It is the basis of further processing.
  * The embedded icc profile and metadata information can be obtained through this class, too.
  */
class InitialImage
{

public:
    /** Returns the file name of the image.
      * @return The file name of the image */
    virtual Glib::ustring getFileName () = 0;
    /** Returns the embedded icc profile of the image.
      * @return The handle of the embedded profile */
    virtual cmsHPROFILE getEmbeddedProfile () = 0;
    /** Returns a class providing access to the exif and iptc metadata tags of all frames of the image.
      * @return An instance of the FramesMetaData class */
    virtual const FramesMetaData* getMetaData () = 0;
    /** This is a function used for internal purposes only. */
    virtual ImageSource* getImageSource () = 0;
    /** This class has manual reference counting. You have to call this function each time to make a new reference to an instance. */
    virtual void increaseRef () {}
    /** This class has manual reference counting. You have to call this function each time to remove a reference
      * (the last one deletes the instance automatically). */
    virtual void decreaseRef () {}

    virtual ~InitialImage () {}

    /** Loads an image into the memory.
      * @param fname the name of the file
      * @param isRaw shall be true if it is a raw file
      * @param errorCode is a pointer to a variable that is set to nonzero if an error happened (output)
      * @param pl is a pointer pointing to an object implementing a progress listener. It can be NULL, in this case progress is not reported.
      * @return an object representing the loaded and pre-processed image */
    static InitialImage* load (const Glib::ustring& fname, bool isRaw, int* errorCode, ProgressListener* pl = nullptr);
};

/** When the preview image is ready for display during staged processing (thus the changes have been updated),
  * the staged processor notifies the listener class implementing a PreviewImageListener.
  * It is important to note that the file passed to the listener can be used in a shared manner (copying it is not
  * needed) as long as the mutex corresponding to the image is used every time the image is accessed.
  * If the scale of the preview image is >1, no sharpening, no denoising and no cropping is applied, and
  * the transform operations (rotate, c/a, vignetting correction) are performed using a faster less quality algorithm.
  * The image you get with this listener is created to display on the monitor (monitor profile has been already applied). */
class PreviewImageListener
{
public:
    virtual ~PreviewImageListener() {}
    /** With this member function the staged processor notifies the listener that it allocated a new
     * image to store the end result of the processing. It can be used in a shared manner.
     * @param img is a pointer to the image
     * @param scale describes the current scaling applied compared to the 100% size (preview scale)
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
class DetailedCropListener
{
public:
    virtual ~DetailedCropListener() {}
    /** With this member function the staged processor notifies the listener that the detailed crop image has been updated.
      * @param img is a pointer to the detailed crop image */
    virtual void setDetailedCrop (IImage8* img, IImage8* imgtrue, procparams::ColorManagementParams cmp,
                                  procparams::CropParams cp, int cx, int cy, int cw, int ch, int skip) {}
    virtual bool getWindow       (int& cx, int& cy, int& cw, int& ch, int& skip)
    {
        return false;
    }
    //  virtual void    setPosition (int x, int y, bool update=true) {}

};

/** This listener is used when the full size of the final image has been changed (e.g. rotated by 90 deg.) */
class SizeListener
{
public:
    virtual ~SizeListener() {}
    /** This member function is called when the size of the final image has been changed
      * @param w is the width of the final image (without cropping)
      * @param h is the height of the final image (without cropping)
      * @param ow is the width of the final image (without resizing and cropping)
      * @param oh is the height of the final image (without resizing and cropping) */
    virtual void sizeChanged (int w, int h, int ow, int oh) {}
};

/** This listener is used when the histogram of the final image has changed. */
class HistogramListener
{
public:
    virtual ~HistogramListener() {}
    /** This member function is called when the histogram of the final image has changed.
      * @param histRed is the array of size 256 containing the histogram of the red channel
      * @param histGreen is the array of size 256 containing the histogram of the green channel
      * @param histBlue is the array of size 256 containing the histogram of the blue channel
      * @param histLuma is the array of size 256 containing the histogram of the luminance channel
      * other for curves backgrounds, histRAW is RAW without colors */
    virtual void histogramChanged (LUTu & histRed, LUTu & histGreen, LUTu & histBlue, LUTu & histLuma, LUTu & histToneCurve, LUTu & histLCurve, LUTu & histCCurve,/* LUTu & histCLurve,LUTu & histLLCurve, */LUTu & histLCAM, LUTu & histCCAM,
                                   LUTu & histRedRaw, LUTu & histGreenRaw, LUTu & histBlueRaw, LUTu & histChroma, LUTu & histLRETI) {}
};

/** This listener is used when the auto exposure has been recomputed (e.g. when the clipping ratio changed). */
class AutoExpListener
{
public:
    virtual ~AutoExpListener() {}
    /** This member function is called when the auto exposure has been recomputed.
      * @param brightness is the new brightness value (in logarithmic scale)
      * @param bright is the new ...
      * @param black is the new black level (measured in absolute pixel data)
      * @param contrast is the new contrast values
      * @param hlcompr is the new highlight recovery amount
      * @param hlcomprthresh is the new threshold for hlcompr
      * @param hlrecons set to true if HighLight Reconstruction is enabled */
    virtual void autoExpChanged (double brightness, int bright, int contrast, int black, int hlcompr, int hlcomprthresh, bool hlrecons) {}

    virtual void autoMatchedToneCurveChanged(procparams::ToneCurveParams::TcMode curveMode, const std::vector<double> &curve) {}
};

class AutoCamListener
{
public :
    virtual ~AutoCamListener() {}
    virtual void autoCamChanged (double ccam, double ccamout) {}
    virtual void adapCamChanged (double cadap) {}
    virtual void ybCamChanged (int yb) {}

};

class AutoChromaListener
{
public :
    virtual ~AutoChromaListener() {}
    virtual void chromaChanged (double autchroma, double autred, double autblue) {}
    virtual void noiseChanged (double nresid, double highresid) {}
    virtual void noiseTilePrev (int tileX, int tileY, int prevX, int prevY, int sizeT, int sizeP) {}

};

class RetinexListener
{
public :
    virtual ~RetinexListener() {}
    virtual void minmaxChanged (double cdma, double cdmin, double mini, double maxi, double Tmean, double Tsigma, double Tmin, double Tmax) {}

};

class AutoColorTonListener
{
public :
    virtual ~AutoColorTonListener() {}
    virtual void autoColorTonChanged (int bwct, int satthres, int satprot) {}
};

class AutoBWListener
{
public :
    virtual ~AutoBWListener() {}
    virtual void BWChanged (double redbw, double greenbw, double bluebw) {}

};

class localListener
{
public :
    virtual ~localListener() {}
    virtual void localChanged (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str,  std::string ex_str, int sp, int maxdat) {}
    virtual void localretChanged (int **datasp, std::string datastr, std::string ll_str, std::string lh_str, std::string cc_str, std::string hh_str, std::string sk_str, std::string ps_str,  std::string ex_str, int sp, int maxdat) {}

};


class AutoWBListener
{
public :
    virtual ~AutoWBListener() = default;
    virtual void WBChanged (double temp, double green) = 0;
};

class FrameCountListener
{
public :
    virtual ~FrameCountListener() = default;
    virtual void FrameCountChanged (int n, int frameNum) = 0;
};

class ImageTypeListener
{
public :
    virtual ~ImageTypeListener() = default;
    virtual void imageTypeChanged (bool isRaw, bool isBayer, bool isXtrans) = 0;
};

class WaveletListener
{
public :
    virtual ~WaveletListener() {}
    virtual void wavChanged (double nlevel) {}

};


/** This class represents a detailed part of the image (looking through a kind of window).
  * It can be created and destroyed with the appropriate members of StagedImageProcessor.
  * Several crops can be assigned to the same image.   */
class DetailedCrop
{
public:
    virtual ~DetailedCrop() {}
    /** Sets the window defining the crop. */
    virtual void setWindow   (int cx, int cy, int cw, int ch, int skip) {}

    /** First try to update (threadless update). If it returns false, make a full update */
    virtual bool tryUpdate  ()
    {
        return false;
    }
    /** Perform a full recalculation of the part of the image corresponding to the crop. */
    virtual void fullUpdate  () {}
    /** Sets the listener of the crop. */
    virtual void setListener (DetailedCropListener* il) {}
    /** Destroys the crop. */
    virtual void destroy () {}
};

/** This is a staged, cached image processing manager with partial image update support.  */
class StagedImageProcessor
{

public:
    /** Returns the inital image corresponding to the image processor.
      * @return the inital image corresponding to the image processor */
    virtual InitialImage* getInitialImage () = 0;
    /** Returns the current processing parameters.
      * @param dst is the location where the image processing parameters are copied (it is assumed that the memory is allocated by the caller) */
    virtual void        getParams (procparams::ProcParams* dst) = 0;
    /** An essential member function. Call this when a setting has been changed. This function returns a pointer to the
      * processing parameters, that you have to update to reflect the changed situation. When ready, call the paramsUpdateReady
      * function to start the image update.
      * @param change is the ID of the changed setting */
    virtual procparams::ProcParams* beginUpdateParams () = 0;
    /** An essential member function. This indicates that you are ready with the update of the processing parameters you got
      * with the beginUpdateParams call, so the image can be updated. This function returns immediately.
      * The image update starts immediately in the background. If it is ready, the result is passed to a PreviewImageListener
      * and to a DetailedCropListener (if enabled). */
    virtual void        endUpdateParams (ProcEvent change) = 0;
    void endUpdateParams(ProcEventCode change) { endUpdateParams(ProcEvent(change)); }
    virtual void        endUpdateParams (int changeFlags) = 0;
    // Starts a minimal update
    virtual void        startProcessing (int changeCode) = 0;
    /** Stops image processing. When it returns, the image processing is already stopped. */
    virtual void        stopProcessing () = 0;
    /** Sets the scale of the preview image. The larger the number is, the faster the image updates are (typical values are 4-5).
      * @param scale is the scale of the preview image */
    virtual void        setPreviewScale (int scale) = 0;
    /** Returns the scale of the preview image.
      * @return the current scale of the preview image */
    virtual int         getPreviewScale () = 0;
    /** Returns the full width of the resulting image (in 1:1 scale).
      * @return the width of the final image */
    virtual int         getFullWidth () = 0;
    /** Returns the full height of the resulting image (in 1:1 scale).
      * @return the height of the final image */
    virtual int         getFullHeight () = 0;
    /** Returns the width of the preview image.
      * @return the width of the preview image */
    virtual int         getPreviewWidth () = 0;
    /** Returns the height of the preview image.
      * @return the height of the preview image */
    virtual int         getPreviewHeight () = 0;

    virtual bool        updateTryLock() = 0;

    virtual void        updateUnLock() = 0;

    /** Creates and returns a Crop instance that acts as a window on the image
      * @param editDataProvider pointer to the EditDataProvider that communicates with the EditSubscriber
      * @return a pointer to the Crop object that handles the image data trough its own pipeline */
    virtual DetailedCrop* createCrop  (::EditDataProvider *editDataProvider, bool isDetailWindow) = 0;

    virtual bool        getAutoWB   (double& temp, double& green, double equal, double tempBias) = 0;
    virtual void        getCamWB    (double& temp, double& green) = 0;
    virtual void        getSpotWB  (int x, int y, int rectSize, double& temp, double& green) = 0;
    virtual void        getAutoCrop (double ratio, int &x, int &y, int &w, int &h) = 0;

    virtual void        saveInputICCReference (const Glib::ustring& fname, bool apply_wb) = 0;

    virtual void        setProgressListener     (ProgressListener* l) = 0;
    virtual void        setSizeListener         (SizeListener* l) = 0;
    virtual void        delSizeListener         (SizeListener* l) = 0;
    virtual void        setAutoExpListener      (AutoExpListener* l) = 0;
    virtual void        setHistogramListener    (HistogramListener *l) = 0;
    virtual void        setPreviewImageListener (PreviewImageListener* l) = 0;
    virtual void        setAutoCamListener      (AutoCamListener* l) = 0;
    virtual void        setFrameCountListener   (FrameCountListener* l) = 0;
    virtual void        setAutoBWListener       (AutoBWListener* l) = 0;
    virtual void        setlocalListener        (localListener* l) = 0;
    virtual void        setAutoWBListener       (AutoWBListener* l) = 0;
    virtual void        setAutoColorTonListener (AutoColorTonListener* l) = 0;
    virtual void        setAutoChromaListener   (AutoChromaListener* l) = 0;
    virtual void        setRetinexListener      (RetinexListener* l) = 0;
    virtual void        setWaveletListener      (WaveletListener* l) = 0;
    virtual void        setImageTypeListener    (ImageTypeListener* l) = 0;

    virtual void        setMonitorProfile       (const Glib::ustring& monitorProfile, RenderingIntent intent) = 0;
    virtual void        getMonitorProfile       (Glib::ustring& monitorProfile, RenderingIntent& intent) const = 0;
    virtual void        setSoftProofing         (bool softProof, bool gamutCheck) = 0;
    virtual void        getSoftProofing         (bool &softProof, bool &gamutCheck) = 0;

    virtual ~StagedImageProcessor () {}

    /** Returns a staged, cached image processing manager supporting partial updates
    * @param initialImage is a loaded and pre-processed initial image
    * @return the staged image processing manager */
    static StagedImageProcessor* create (InitialImage* initialImage);
    static void destroy (StagedImageProcessor* sip);
};


/**
  * @brief Initializes the RT engine
  * @param s is a struct of basic settings
  * @param baseDir base directory of RT's installation dir
  * @param userSettingsDir RT's base directory in the user's settings dir
  * @param loadAll if false, don't load the various dependencies (profiles, HALDClut files, ...), they'll be loaded from disk each time they'll be used (launching time improvement) */
int init (const Settings* s, Glib::ustring baseDir, Glib::ustring userSettingsDir, bool loadAll = true);

/** Cleanup the RT engine (static variables) */
void cleanup ();

/** This class  holds all the necessary informations to accomplish the full processing of the image */
class ProcessingJob
{

public:
    virtual ~ProcessingJob() {}

    /** Creates a processing job from a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately.
       * @param fname the name of the file
       * @param isRaw shall be true if it is a raw file
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */
    static ProcessingJob* create (const Glib::ustring& fname, bool isRaw, const procparams::ProcParams& pparams, bool fast = false);

    /** Creates a processing job from a file name. This function always succeeds. It only stores the data into the ProcessingJob class, it does not load
       * the image thus it returns immediately. This function increases the reference count of the initialImage. If you decide not the process the image you
       * have to cancel it by calling the member function void cancel(). If the image is processed the reference count of initialImage is decreased automatically, thus the ProcessingJob
       * instance gets invalid. You can not use a ProcessingJob instance to process an image twice.
       * @param initialImage is a loaded and pre-processed initial image
       * @param pparams is a struct containing the processing parameters
       * @return an object containing the data above. It can be passed to the functions that do the actual image processing. */
    static ProcessingJob* create (InitialImage* initialImage, const procparams::ProcParams& pparams, bool fast = false);

    /** Cancels and destroys a processing job. The reference count of the corresponding initialImage (if any) is decreased. After the call of this function the ProcessingJob instance
      * gets invalid, you must not use it any more. Dont call this function while the job is being processed.
      * @param job is the job to destroy */
    static void destroy (ProcessingJob* job);

    virtual bool fastPipeline() const = 0;
};

/** This function performs all the image processinf steps corresponding to the given ProcessingJob. It returns when it is ready, so it can be slow.
   * The ProcessingJob passed becomes invalid, you can not use it any more.
   * @param job the ProcessingJob to cancel.
   * @param errorCode is the error code if an error occured (e.g. the input image could not be loaded etc.)
   * @param pl is an optional ProgressListener if you want to keep track of the progress
   * @return the resulting image, with the output profile applied, exif and iptc data set. You have to save it or you can access the pixel data directly.  */
IImagefloat* processImage (ProcessingJob* job, int& errorCode, ProgressListener* pl = nullptr, bool flush = false);

/** This class is used to control the batch processing. The class implementing this interface will be called when the full processing of an
   * image is ready and the next job to process is needed. */
class BatchProcessingListener : public ProgressListener
{
public:
    /** This function is called when an image gets ready during the batch processing. It has to return with the next job, or with NULL if
                   * there is no jobs left.
                   * @param img is the result of the last ProcessingJob
                   * @return the next ProcessingJob to process */
    virtual ProcessingJob* imageReady (IImagefloat* img) = 0;
    virtual void error (Glib::ustring message) = 0;
};
/** This function performs all the image processinf steps corresponding to the given ProcessingJob. It runs in the background, thus it returns immediately,
   * When it finishes, it calls the BatchProcessingListener with the resulting image and asks for the next job. It the listener gives a new job, it goes on
   * with processing. If no new job is given, it finishes.
   * The ProcessingJob passed becomes invalid, you can not use it any more.
   * @param job the ProcessingJob to cancel.
   * @param bpl is the BatchProcessingListener that is called when the image is ready or the next job is needed. It also acts as a ProgressListener.
   **/
void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl);


extern MyMutex* lcmsMutex;
}

#endif


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
#ifndef __IMAGEDATA_H__
#define __IMAGEDATA_H__

#include <stdio.h>
#include <rawimage.h>
#include <string>
#include <glibmm.h>
#include <procparams.h>
#include <rtXmp.h>
#include <map>

namespace rtengine {

class SnapshotInfo{
	public:
		Glib::ustring   name; // User name of the snapshot
		int				id;   // identifier for the snapshot
		Glib::ustring   outputFilename; // If saved or in batch queue, this contains the relative path to the output file
		bool			queued; // processing pending
		bool			saved; // snapshot saved
		procparams::ProcParams params; // Processing parameters

		static const char *kCurrentSnapshotName;
};

typedef std::map< unsigned int , SnapshotInfo > snapshotsList_t;

typedef std::map< std::string, std::vector<Glib::ustring> > MetadataList;

class ExifPair
{
public:
	Glib::ustring group;
	Glib::ustring field;
	Glib::ustring value;
	Glib::ustring rawValue;
};

inline bool operator<(const ExifPair &v1, const ExifPair &v2 )
{
	int b = v1.group.compare( v2.group );
	if( b< 0 )
		return true;
	else if( b == 0 && v1.field.compare( v2.field )<0 )
		return true;
	return false;
}

class ImageMetaData {

  protected:
    Exiv2::XmpData  xmpData;  // This contains all metadata: EXIF(partial) and IPTC extracted and more XMP tags added later
    Exiv2::ExifData exifData; // This contains all EXIF data read from image (included makernotes): is valid only if exifExtracted is true
    Exiv2::IptcData iptcData; // Buffer for old IIM iptc data
    bool exifExtracted;       // Exif are fully read only first time, when xmp data is written with most important tiff/exif tags
    bool xmpEmbedded;         // Write xmp inside image file
    Glib::ustring fname;      // Filename of the original photo.
    Glib::ustring fnameMeta;  // Filename of persisted metadata
    Glib::ustring fnameMeta2; // Backup of persisted metadata  
  
    /** Init xmp RT data */
    int initRTXMP( );

    /** (Re)Init metadata from image */
    int readMetadataFromImage();

  public:
    // Equal to the definition in TIFF/EXIF tag
    typedef enum
    {
    	ePhotoOrientationNormal=1,
    	ePhotoOrientationMirrorHoriz=2,
    	ePhotoOrientationRotate180=3,
    	ePhotoOrientationMirrorVert=4,
    	ePhotoOrientationMirHorz270=5,
    	ePhotoOrientationRotate90=6,
    	ePhotoOrientationMirHorz90=7,
    	ePhotoOrientationRotate270=8
    }ePhotoOrientation;    
  

    ImageMetaData (const Glib::ustring &fname, const Glib::ustring &fnameMeta, const Glib::ustring &fnameMeta2="" );
    ImageMetaData ( const ImageMetaData &v);
    ~ImageMetaData ();

    /** Prepare metadata for writing */
    void            merge( bool syncExif=true, bool syncIPTC=true, bool removeProcessing=true );

    /** Set output image informations in xmp tags */
    void setOutputImageInfo(int w, int h, int bits );

    /** Write metadata to image */
    bool            writeToImage( const Glib::ustring &fname, bool writeExif=true, bool writeIPTC=true, bool writeXmp=true ) const;

    /** Load xmp data from sidecar file */
    int             loadXMP( );

    /** Save xmp data into sidecar file */
    int             saveXMP( ) const;

    /** Re-read metadata embedded in imagefile (delete previously modified data)*/
    int             resync( );

    /** Create a new snapshot with given name */
	int 			newSnapshot(const Glib::ustring &name, const rtengine::procparams::ProcParams& params, bool queued=false );

	/** Delete a snapshot */
	bool 			deleteSnapshot( int id );

	/** Delete all procparams saved inside XMP */
	bool            deleteAllSnapshots();

	/** Rename a snapshot*/
	bool 			renameSnapshot(int id, const Glib::ustring &newname );

	/** Return a l ist of all readable snapshots contained in xmp */
	snapshotsList_t getSnapshotsList();

	/** Get a snpashot by id*/
	SnapshotInfo    getSnapshot( int id );

	/** Get a snapshot identifier by its name*/
	int             getSnapshotId( const Glib::ustring &snapshotName=SnapshotInfo::kCurrentSnapshotName );

	/** Update porcessing parameters of a snapshot (usually current))*/
	bool            updateSnapshot( int id, const rtengine::procparams::ProcParams& params);

	/** Return number of snapshot queued*/
	int				getNumQueued();

	/** Set/reset the flag queued*/
	bool			setQueuedSnapshot( int id, bool inqueue=true );

	/** Set/reset the flag saved adn optionally the filename of the output file*/
	bool			setSavedSnapshot( int id, bool saved=true, const Glib::ustring &filename="" );

	/** Full EXIF is read only when thumbnail is loaded first time, so, to access full EXIF call updateExif*/
	int             updateExif();

	/** call updateExif() before getting ExifData */
	std::vector<ExifPair>     getExifData () const;

	/** return all IPTC metadata */
    const rtengine::MetadataList getIPTCData () const;

    /** return true if changes have been applied */
    bool  getIPTCDataChanged() const;

    /** change IPTC data */
    void setIPTCData( const rtengine::MetadataList &meta );

    /** if embedded, RT saves metadata inside source image and not in sidecar xmp file */
    bool isEmbedded() const { return xmpEmbedded; }

    /** @return a struct containing the date and time of the image */
    struct tm   getDateTime ();

    /** @return a timestamp containing the date and time of the image */
    time_t      getDateTimeAsTS();

    /** @return the ISO of the image */
    int         getISOSpeed ();

    /** @return the F number of the image */
    double      getFNumber  ();

    /** @return the focal length used at the exposure */
    double      getFocalLen () ;

    /** @return the shutter speed */
    double      getShutterSpeed ();

    /** @return the exposure compensation */
    double      getExpComp  ();

    /** @return if flash fired */
    bool        getFlashFired();

    /** @return the maker of the camera */
    std::string getMake     ();

    /** @return the model of the camera */
    std::string getModel    ();

    /** @return the complete name of the camera for displaying purpose: Maker and Model are simplified like dcraw */
    std::string getCamera   ();

    /** @return the lens on the camera  */
    std::string getLens     ();

    /** @return the serial number of the camera if available */
    std::string getSerialNumber ();

    /** @return the orientation of the image */
    ePhotoOrientation getOrientation ();

    /** Rank of the photo 0..1, -1= rejected*/
    void        setRank  (int rank);
    int         getRank  ();

    /** Label is used to stick a color flag */
    void        setLabel( const Glib::ustring &colorlabel );
    Glib::ustring getLabel();

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

    /** Reads metadata from file.
      * @param fname is the name of the file
      * @param rml is a struct containing information about metadata location. Use it only for raw files. In case
      * of jpgs and tiffs pass a NULL pointer.
      * @return The metadata */
      static ImageMetaData* fromFile (const Glib::ustring& fname, const Glib::ustring &fnameMeta, const Glib::ustring &fnameMeta2="", bool embed=false);
};

};
#endif

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
#include "cacheimagedata.h"
#include <vector>
#include <glib/gstdio.h>
#include "version.h"
#include <locale.h>

CacheImageData::CacheImageData ()
    : md5(""), supported(false), format(FT_Invalid), rankOld(-1), inTrashOld(false), recentlySaved(false),
      timeValid(false), year(0), month(0), day(0), hour(0), min(0), sec(0), exifValid(false),
      fnumber(0.0), shutter(0.0), focalLen(0.0), focalLen35mm(0.0), focusDist(0.f), iso(0), isHDR (false), isPixelShift (false),
      sampleFormat(rtengine::IIOSF_UNKNOWN), redAWBMul(-1.0), greenAWBMul(-1.0), blueAWBMul(-1.0), rotate(0), thumbImgType(0)
{
}

/*
 * Load the General, DateTime, ExifInfo, File info and ExtraRawInfo sections of the image data file
 */
int CacheImageData::load (const Glib::ustring& fname)
{
    setlocale(LC_NUMERIC, "C"); // to set decimal point to "."

    Glib::KeyFile keyFile;

    try {
        if (keyFile.load_from_file (fname)) {

            if (keyFile.has_group ("General")) {
                if (keyFile.has_key ("General", "MD5")) {
                    md5         = keyFile.get_string ("General", "MD5");
                }

                if (keyFile.has_key ("General", "Version")) {
                    version     = keyFile.get_string ("General", "Version");
                }

                if (keyFile.has_key ("General", "Supported")) {
                    supported   = keyFile.get_boolean ("General", "Supported");
                }

                if (keyFile.has_key ("General", "Format")) {
                    format      = (ThFileType)keyFile.get_integer ("General", "Format");
                }

                if (keyFile.has_key ("General", "Rank")) {
                    rankOld     = keyFile.get_integer ("General", "Rank");
                }

                if (keyFile.has_key ("General", "InTrash")) {
                    inTrashOld  = keyFile.get_boolean ("General", "InTrash");
                }

                if (keyFile.has_key ("General", "RecentlySaved")) {
                    recentlySaved = keyFile.get_boolean ("General", "RecentlySaved");
                }
            }

            timeValid = keyFile.has_group ("DateTime");

            if (timeValid) {
                if (keyFile.has_key ("DateTime", "Year")) {
                    year    = keyFile.get_integer ("DateTime", "Year");
                }

                if (keyFile.has_key ("DateTime", "Month")) {
                    month   = keyFile.get_integer ("DateTime", "Month");
                }

                if (keyFile.has_key ("DateTime", "Day")) {
                    day     = keyFile.get_integer ("DateTime", "Day");
                }

                if (keyFile.has_key ("DateTime", "Hour")) {
                    hour    = keyFile.get_integer ("DateTime", "Hour");
                }

                if (keyFile.has_key ("DateTime", "Min")) {
                    min     = keyFile.get_integer ("DateTime", "Min");
                }

                if (keyFile.has_key ("DateTime", "Sec")) {
                    sec     = keyFile.get_integer ("DateTime", "Sec");
                }
            }

            exifValid = false;

            if (keyFile.has_group ("ExifInfo")) {
                exifValid = true;

                if (keyFile.has_key ("ExifInfo", "Valid")) {
                    exifValid = keyFile.get_boolean ("ExifInfo", "Valid");
                }

                if (exifValid) {
                    if (keyFile.has_key ("ExifInfo", "FNumber")) {
                        fnumber     = keyFile.get_double ("ExifInfo", "FNumber");
                    }

                    if (keyFile.has_key ("ExifInfo", "Shutter")) {
                        shutter     = keyFile.get_double ("ExifInfo", "Shutter");
                    }

                    if (keyFile.has_key ("ExifInfo", "FocalLen")) {
                        focalLen    = keyFile.get_double ("ExifInfo", "FocalLen");
                    }

                    if (keyFile.has_key ("ExifInfo", "FocalLen35mm")) {
                        focalLen35mm = keyFile.get_double ("ExifInfo", "FocalLen35mm");
                    } else {
                        focalLen35mm = focalLen;    // prevent crashes on old files
                    }

                    if (keyFile.has_key ("ExifInfo", "FocusDist")) {
                        focusDist = keyFile.get_double ("ExifInfo", "FocusDist");
                    } else {
                        focusDist = 0;
                    }

                    if (keyFile.has_key ("ExifInfo", "ISO")) {
                        iso         = keyFile.get_integer ("ExifInfo", "ISO");
                    }

                    if (keyFile.has_key ("ExifInfo", "IsHDR")) {
                        isHDR = keyFile.get_integer ("ExifInfo", "IsHDR");
                    }

                    if (keyFile.has_key ("ExifInfo", "IsPixelShift")) {
                        isPixelShift = keyFile.get_boolean ("ExifInfo", "IsPixelShift");
                    }

                    if (keyFile.has_key ("ExifInfo", "ExpComp")) {
                        expcomp     = keyFile.get_string ("ExifInfo", "ExpComp");
                    }
                }

                if (keyFile.has_key ("ExifInfo", "Lens")) {
                    lens        = keyFile.get_string ("ExifInfo", "Lens");
                }

                if (keyFile.has_key ("ExifInfo", "CameraMake")) {
                    camMake     = keyFile.get_string ("ExifInfo", "CameraMake");
                }

                if (keyFile.has_key ("ExifInfo", "CameraModel")) {
                    camModel    = keyFile.get_string ("ExifInfo", "CameraModel");
                }
            }

            if (keyFile.has_group ("FileInfo")) {
                if (keyFile.has_key ("FileInfo", "Filetype")) {
                    filetype    = keyFile.get_string ("FileInfo", "Filetype");
                }
                if (keyFile.has_key ("FileInfo", "SampleFormat")) {
                    sampleFormat = (rtengine::IIO_Sample_Format)keyFile.get_integer ("FileInfo", "SampleFormat");
                }
            }

            if (format == FT_Raw && keyFile.has_group ("ExtraRawInfo")) {
                if (keyFile.has_key ("ExtraRawInfo", "ThumbImageType")) {
                    thumbImgType    = keyFile.get_integer ("ExtraRawInfo", "ThumbImageType");
                }
            } else {
                rotate = 0;
                thumbImgType = 0;
            }

            return 0;
        }
    } catch (Glib::Error &err) {
        if (options.rtSettings.verbose) {
            printf("CacheImageData::load / Error code %d while reading values from \"%s\":\n%s\n", err.code(), fname.c_str(), err.what().c_str());
        }
    } catch (...) {
        if (options.rtSettings.verbose) {
            printf("CacheImageData::load / Unknown exception while trying to load \"%s\"!\n", fname.c_str());
        }
    }

    return 1;
}

/*
 * Save the General, DateTime, ExifInfo, File info and ExtraRawInfo sections of the image data file
 */
int CacheImageData::save (const Glib::ustring& fname)
{

    Glib::ustring keyData;

    try {

    Glib::KeyFile keyFile;

    try {
        keyFile.load_from_file (fname);
    } catch (Glib::Error&) {}

    keyFile.set_string  ("General", "MD5", md5);
    keyFile.set_string  ("General", "Version", RTVERSION);
    keyFile.set_boolean ("General", "Supported", supported);
    keyFile.set_integer ("General", "Format", format);
    keyFile.set_boolean ("General", "RecentlySaved", recentlySaved);

    // remove the old implementation of Rank and InTrash from cache
    if (keyFile.has_key ("General", "Rank")) {
        keyFile.remove_key("General", "Rank");
    }

    if (keyFile.has_key ("General", "InTrash")) {
        keyFile.remove_key("General", "InTrash");
    }

    if (timeValid) {
        keyFile.set_integer ("DateTime", "Year", year);
        keyFile.set_integer ("DateTime", "Month", month);
        keyFile.set_integer ("DateTime", "Day", day);
        keyFile.set_integer ("DateTime", "Hour", hour);
        keyFile.set_integer ("DateTime", "Min", min);
        keyFile.set_integer ("DateTime", "Sec", sec);
    }

    keyFile.set_boolean  ("ExifInfo", "Valid", exifValid);

    if (exifValid) {
        keyFile.set_double  ("ExifInfo", "FNumber", fnumber);
        keyFile.set_double  ("ExifInfo", "Shutter", shutter);
        keyFile.set_double  ("ExifInfo", "FocalLen", focalLen);
        keyFile.set_double  ("ExifInfo", "FocalLen35mm", focalLen35mm);
        keyFile.set_double  ("ExifInfo", "FocusDist", focusDist);
        keyFile.set_integer ("ExifInfo", "ISO", iso);
        keyFile.set_boolean ("ExifInfo", "IsHDR", isHDR);
        keyFile.set_boolean ("ExifInfo", "IsPixelShift", isPixelShift);
        keyFile.set_string  ("ExifInfo", "ExpComp", expcomp);
    }

    keyFile.set_string  ("ExifInfo", "Lens", lens);
    keyFile.set_string  ("ExifInfo", "CameraMake", camMake);
    keyFile.set_string  ("ExifInfo", "CameraModel", camModel);
    keyFile.set_string  ("FileInfo", "Filetype", filetype);
    keyFile.set_integer ("FileInfo", "SampleFormat", sampleFormat);

    if (format == FT_Raw) {
        keyFile.set_integer ("ExtraRawInfo", "ThumbImageType", thumbImgType);
    }

    keyData = keyFile.to_data ();

    } catch (Glib::Error &err) {
        if (options.rtSettings.verbose) {
            printf("CacheImageData::save / Error code %d while reading values from \"%s\":\n%s\n", err.code(), fname.c_str(), err.what().c_str());
        }
    } catch (...) {
        if (options.rtSettings.verbose) {
            printf("CacheImageData::save / Unknown exception while trying to save \"%s\"!\n", fname.c_str());
        }
    }

    if (keyData.empty ()) {
        return 1;
    }

    FILE *f = g_fopen (fname.c_str (), "wt");

    if (!f) {
        if (options.rtSettings.verbose) {
            printf("CacheImageData::save / Error: unable to open file \"%s\" with write access!\n", fname.c_str());
        }

        return 1;
    } else {
        fprintf (f, "%s", keyData.c_str ());
        fclose (f);
        return 0;
    }
}


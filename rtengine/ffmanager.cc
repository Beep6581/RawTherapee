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
#include "ffmanager.h"
#include "../rtgui/options.h"
#include "rawimage.h"
#include "imagedata.h"
#include "median.h"
#include "utils.h"

namespace rtengine
{

extern const Settings* settings;

// *********************** class ffInfo **************************************

inline ffInfo& ffInfo::operator =(const ffInfo &o)
{
    if (this != &o) {
        pathname = o.pathname;
        maker = o.maker;
        model = o.model;
        lens = o.lens;
        focallength = o.focallength;
        timestamp = o.timestamp;
        aperture = o.aperture;

        if( ri ) {
            delete ri;
            ri = nullptr;
        }
    }

    return *this;
}

bool ffInfo::operator <(const ffInfo &e2) const
{
    if( this->maker.compare( e2.maker) >= 0 ) {
        return false;
    }

    if( this->model.compare( e2.model) >= 0 ) {
        return false;
    }

    if( this->lens.compare( e2.lens) >= 0 ) {
        return false;
    }

    if( this->focallength >= e2.focallength ) {
        return false;
    }

    if( this->timestamp >= e2.timestamp ) {
        return false;
    }

    return true;
}

std::string ffInfo::key(const std::string &mak, const std::string &mod, const std::string &len, double focal, double apert )
{
    std::ostringstream s;
    s << mak << " " << mod << " ";
    s.width(5);
    s << len << " ";
    s.precision( 2 );
    s.width(4);
    s << focal << "mm F" << apert;
    return s.str();
}

double ffInfo::distance(const std::string &mak, const std::string &mod, const std::string &len, double focallength, double aperture) const
{
    if( this->maker.compare( mak) != 0 ) {
        return INFINITY;
    }

    if( this->model.compare( mod) != 0 ) {
        return INFINITY;
    }

    if( this->lens.compare( len) != 0 ) {
        return INFINITY;
    }

    double dAperture = 2 * (log(this->aperture) - log(aperture)) / log(2); //more important for vignette
    double dfocallength = (log(this->focallength / 100.) - log(focallength / 100.)) / log(2); //more important for PRNU

    return sqrt( dfocallength * dfocallength + dAperture * dAperture);
}

RawImage* ffInfo::getRawImage()
{
    if(ri) {
        return ri;
    }

    updateRawImage();

    return ri;
}

/* updateRawImage() load into ri the actual pixel data from pathname if there is a single shot
 * otherwise load each file from the pathNames list and extract a template from the media;
 * the first file is used also for reading all information other than pixels
 */
void ffInfo::updateRawImage()
{
    typedef unsigned int acc_t;

    // averaging of flatfields if more than one is found matching the same key.
    // this may not be necessary, as flatfield is further blurred before being applied to the processed image.
    if( !pathNames.empty() ) {
        std::list<Glib::ustring>::iterator iName = pathNames.begin();
        ri = new RawImage(*iName, nullptr); // First file used also for extra pixels informations (width,height, shutter, filters etc.. )
        if( ri->loadRaw(true, false)) {
            delete ri;
            ri = nullptr;
        } else {
            int H = ri->get_height();
            int W = ri->get_width();
            ri->compress_image(0);
            int rSize = W * ((ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS) ? 1 : 3);
            acc_t **acc = new acc_t*[H];

            for( int row = 0; row < H; row++) {
                acc[row] = new acc_t[rSize ];
            }

            // copy first image into accumulators
            for (int row = 0; row < H; row++)
                for (int col = 0; col < rSize; col++) {
                    acc[row][col] = ri->data[row][col];
                }

            int nFiles = 1; // First file data already loaded

            for( ++iName; iName != pathNames.end(); ++iName) {
                RawImage* temp = new RawImage(*iName, nullptr);

                if( !temp->loadRaw(true, false)) {
                    temp->compress_image(0);     //\ TODO would be better working on original, because is temporary
                    nFiles++;

                    if( ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS ) {
                        for( int row = 0; row < H; row++) {
                            for( int col = 0; col < W; col++) {
                                acc[row][col] += temp->data[row][col];
                            }
                        }
                    } else {
                        for( int row = 0; row < H; row++) {
                            for( int col = 0; col < W; col++) {
                                acc[row][3 * col + 0] += temp->data[row][3 * col + 0];
                                acc[row][3 * col + 1] += temp->data[row][3 * col + 1];
                                acc[row][3 * col + 2] += temp->data[row][3 * col + 2];
                            }
                        }
                    }
                }

                delete temp;
            }

            for (int row = 0; row < H; row++) {
                for (int col = 0; col < rSize; col++) {
                    ri->data[row][col] = acc[row][col] / nFiles;
                }

                delete [] acc[row];
            }

            delete [] acc;
        }
    } else {
        ri = new RawImage(pathname, nullptr);
        if( ri->loadRaw(true, false)) {
            delete ri;
            ri = nullptr;
        } else {
            ri->compress_image(0);
        }
    }

    if(ri) {
        // apply median to avoid this step being executed each time a flat field gets applied
        int H = ri->get_height();
        int W = ri->get_width();
        float *cfatmp = (float (*)) malloc (H * W * sizeof * cfatmp);

#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic,16)
#endif

        for (int i = 0; i < H; i++) {
            int iprev = i < 2 ? i + 2 : i - 2;
            int inext = i > H - 3 ? i - 2 : i + 2;

            for (int j = 0; j < W; j++) {
                int jprev = j < 2 ? j + 2 : j - 2;
                int jnext = j > W - 3 ? j - 2 : j + 2;

                cfatmp[i * W + j] = median(ri->data[iprev][j], ri->data[i][jprev], ri->data[i][j], ri->data[i][jnext], ri->data[inext][j]);
            }
        }

        memcpy(ri->data[0], cfatmp, W * H * sizeof(float));

        free (cfatmp);

    }
}

// ************************* class FFManager *********************************

void FFManager::init( Glib::ustring pathname )
{
    std::vector<Glib::ustring> names;

    auto dir = Gio::File::create_for_path (pathname);

    if (!dir || !dir->query_exists()) {
        return;
    }

    try {

        auto enumerator = dir->enumerate_children ("standard::name");

        while (auto file = enumerator->next_file ()) {
            names.emplace_back (Glib::build_filename (pathname, file->get_name ()));
        }

    } catch (Glib::Exception&) {}

    ffList.clear();

    for (size_t i = 0; i < names.size(); i++) {
        try {
            addFileInfo(names[i]);
        } catch( std::exception& e ) {}
    }

    // Where multiple shots exist for same group, move filename to list
    for( ffList_t::iterator iter = ffList.begin(); iter != ffList.end(); ++iter ) {
        ffInfo &i = iter->second;

        if( !i.pathNames.empty() && !i.pathname.empty() ) {
            i.pathNames.push_back( i.pathname );
            i.pathname.clear();
        }

        if( settings->verbose ) {
            if( !i.pathname.empty() ) {
                printf( "%s:  %s\n", i.key().c_str(), i.pathname.c_str());
            } else {
                printf( "%s: MEAN of \n    ", i.key().c_str());

                for( std::list<Glib::ustring>::iterator iter = i.pathNames.begin(); iter != i.pathNames.end(); ++iter  ) {
                    printf( "%s, ", iter->c_str() );
                }

                printf("\n");
            }
        }
    }

    currentPath = pathname;
    return;
}

ffInfo* FFManager::addFileInfo (const Glib::ustring& filename, bool pool)
{
    auto ext = getFileExtension(filename);

    if (ext.empty() || !options.is_extention_enabled(ext)) {
        return nullptr;
    }

    auto file = Gio::File::create_for_path(filename);

    if (!file ) {
        return nullptr;
    }

    if (!file->query_exists()) {
        return nullptr;
    }

    try {

        auto info = file->query_info("standard::name,standard::type,standard::is-hidden");

        if (!info || info->get_file_type() == Gio::FILE_TYPE_DIRECTORY) {
            return nullptr;
        }

        if (!options.fbShowHidden && info->is_hidden()) {
            return nullptr;
        }

        RawImage ri(filename, nullptr);
        int res = ri.loadRaw(false, false); // Read informations about shot

        if (res != 0) {
            return nullptr;
        }

        ffList_t::iterator iter;

        if(!pool) {
            ffInfo n(filename, "", "", "", 0, 0, 0);
            iter = ffList.emplace("", n);
            return &(iter->second);
        }

        FramesData idata(filename, std::unique_ptr<RawMetaDataLocation>(new RawMetaDataLocation(ri.get_exifBase(), ri.get_ciffBase(), ri.get_ciffLen())), true);
        /* Files are added in the map, divided by same maker/model,lens and aperture*/
        std::string key(ffInfo::key(idata.getMake(), idata.getModel(), idata.getLens(), idata.getFocalLen(), idata.getFNumber()));
        iter = ffList.find(key);

        if(iter == ffList.end()) {
            ffInfo n(filename, idata.getMake(), idata.getModel(), idata.getLens(), idata.getFocalLen(), idata.getFNumber(), idata.getDateTimeAsTS());
            iter = ffList.emplace(key, n);
        } else {
            while(iter != ffList.end() && iter->second.key() == key && ABS(iter->second.timestamp - ri.get_timestamp()) > 60 * 60 * 6) { // 6 hour difference
                ++iter;
            }

            if(iter != ffList.end()) {
                iter->second.pathNames.push_back(filename);
            } else {
                ffInfo n(filename, idata.getMake(), idata.getModel(), idata.getLens(), idata.getFocalLen(), idata.getFNumber(), idata.getDateTimeAsTS());
                iter = ffList.emplace(key, n);
            }
        }

        return &(iter->second);

    } catch (Gio::Error&) {}

    return nullptr;
}

void FFManager::getStat( int &totFiles, int &totTemplates)
{
    totFiles = 0;
    totTemplates = 0;

    for( ffList_t::iterator iter = ffList.begin(); iter != ffList.end(); ++iter ) {
        ffInfo &i = iter->second;

        if( i.pathname.empty() ) {
            totTemplates++;
            totFiles += i.pathNames.size();
        } else {
            totFiles++;
        }
    }
}

/*  The search for the best match is twofold:
 *  if perfect matches for make and model are found, then the list is scanned for lesser distance in time
 *  otherwise if no match is found, the whole list is searched for lesser distance in lens and aperture
 */
ffInfo* FFManager::find( const std::string &mak, const std::string &mod, const std::string &len, double focal, double apert, time_t t )
{
    if( ffList.empty() ) {
        return nullptr;
    }

    std::string key( ffInfo::key(mak, mod, len, focal, apert) );
    ffList_t::iterator iter = ffList.find( key );

    if(  iter != ffList.end() ) {
        ffList_t::iterator bestMatch = iter;
        time_t bestDeltaTime = ABS(iter->second.timestamp - t);

        for(++iter; iter != ffList.end() && !key.compare( iter->second.key() ); ++iter ) {
            time_t d = ABS(iter->second.timestamp - t );

            if( d < bestDeltaTime ) {
                bestMatch = iter;
                bestDeltaTime = d;
            }
        }

        return &(bestMatch->second);
    } else {
        iter = ffList.begin();
        ffList_t::iterator bestMatch = iter;
        double bestD = iter->second.distance(  mak, mod, len, focal, apert );

        for( ++iter; iter != ffList.end(); ++iter ) {
            double d = iter->second.distance(  mak, mod, len, focal, apert );

            if( d < bestD ) {
                bestD = d;
                bestMatch = iter;
            }
        }

        return bestD != INFINITY ? &(bestMatch->second) : nullptr ;
    }
}

RawImage* FFManager::searchFlatField( const std::string &mak, const std::string &mod, const std::string &len, double focal, double apert, time_t t )
{
    ffInfo *ff = find( mak, mod, len, focal, apert, t );

    if( ff ) {
        return ff->getRawImage();
    } else {
        return nullptr;
    }
}

RawImage* FFManager::searchFlatField( const Glib::ustring filename )
{
    for ( ffList_t::iterator iter = ffList.begin(); iter != ffList.end(); ++iter ) {
        if( iter->second.pathname.compare( filename ) == 0  ) {
            return iter->second.getRawImage();
        }
    }

    ffInfo *ff = addFileInfo( filename , false);

    if(ff) {
        return ff->getRawImage();
    }

    return nullptr;
}


// Global variable
FFManager ffm;


}


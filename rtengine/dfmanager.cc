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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <cstdio>
#include <iostream>
#include <list>
#include <map>
#include <sstream>

#include <giomm.h>
#include <glibmm/ustring.h>

#include "dfmanager.h"

#include "imagedata.h"
#include "jaggedarray.h"
#include "noncopyable.h"
#include "pixelsmap.h"
#include "rawimage.h"
#include "utils.h"

#include "../rtgui/options.h"

namespace
{

std::string toUppercase(const std::string& string)
{
    return Glib::ustring(string).uppercase();
}

class dfInfo final
{
public:
    Glib::ustring pathname; // filename of dark frame
    std::list<Glib::ustring> pathNames; // other similar dark frames, used for average
    std::string maker; // manufacturer
    std::string model; // model
    int iso; // ISO (gain)
    double shutter; // shutter or exposure time in sec
    time_t timestamp; // seconds since 1 Jan 1970


    dfInfo(const Glib::ustring &name, const std::string &mak, const std::string &mod, int iso, double shut, time_t t)
        : pathname(name), maker(mak), model(mod), iso(iso), shutter(shut), timestamp(t), ri(nullptr) {}

    dfInfo(const dfInfo &o)
        : pathname(o.pathname), maker(o.maker), model(o.model), iso(o.iso), shutter(o.shutter), timestamp(o.timestamp), ri(nullptr) {}
    ~dfInfo();

    dfInfo &operator =(const dfInfo &o);

    // Calculate virtual distance between two shots; different model return infinite
    double distance(const std::string &mak, const std::string &mod, int iso, double shutter) const;

    static std::string key(const std::string &mak, const std::string &mod, int iso, double shut);
    std::string key() const
    {
        return key(maker, model, iso, shutter);
    }

    const rtengine::RawImage* getRawImage();
    const std::vector<rtengine::badPix>& getHotPixels();

private:
    rtengine::RawImage* ri; // Dark Frame raw data
    std::vector<rtengine::badPix> badPixels; // Extracted hot pixels

    void updateBadPixelList(const rtengine::RawImage* df);
    void updateRawImage();
};

dfInfo::~dfInfo()
{
    delete ri;
}

inline dfInfo& dfInfo::operator = (const dfInfo &o)
{
    if (this != &o) {
        pathname = o.pathname;
        maker = o.maker;
        model = o.model;
        iso = o.iso;
        shutter = o.shutter;
        timestamp = o.timestamp;

        if (ri) {
            delete ri;
            ri = nullptr;
        }
    }

    return *this;
}

std::string dfInfo::key(const std::string &mak, const std::string &mod, int iso, double shut)
{
    std::ostringstream s;
    s << mak << " " << mod << " ";
    s.width(5);
    s << iso << "ISO ";
    s.precision(2);
    s.width(4);
    s << shut << "s";
    return s.str();
}

double dfInfo::distance(const std::string &mak, const std::string &mod, int iso, double shutter) const
{
    if (this->maker.compare(mak) != 0) {
        return INFINITY;
    }

    if (this->model.compare(mod) != 0) {
        return INFINITY;
    }

    const double dISO = (log(this->iso / 100.) - log(iso / 100.)) / log(2);
    const double dShutter = (log(this->shutter) - log(shutter)) / log(2);
    return std::sqrt(dISO * dISO +  dShutter * dShutter);
}

const rtengine::RawImage* dfInfo::getRawImage()
{
    if (ri) {
        return ri;
    }

    updateRawImage();
    updateBadPixelList(ri);

    return ri;
}

const std::vector<rtengine::badPix>& dfInfo::getHotPixels()
{
    if (!ri) {
        updateRawImage();
        updateBadPixelList(ri);
    }

    return badPixels;
}

/* updateRawImage() load into ri the actual pixel data from pathname if there is a single shot
 * otherwise load each file from the pathNames list and extract a template from the media;
 * the first file is used also for reading all information other than pixels
 */
void dfInfo::updateRawImage()
{

    if (!pathNames.empty()) {
        std::list<Glib::ustring>::const_iterator iName = pathNames.begin();
        ri = new rtengine::RawImage(*iName); // First file used also for extra pixels information (width,height, shutter, filters etc.. )

        if (ri->loadRaw(true)) {
            delete ri;
            ri = nullptr;
        } else {
            const int H = ri->get_height();
            const int W = ri->get_width();
            ri->compress_image(0);
            const int rSize = W * ((ri->getSensorType() == rtengine::ST_BAYER || ri->getSensorType() == rtengine::ST_FUJI_XTRANS) ? 1 : 3);
            rtengine::JaggedArray<float> acc(W, H);

            // copy first image into accumulators
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < rSize; col++) {
                    acc[row][col] = ri->data[row][col];
                }
            }

            int nFiles = 1; // First file data already loaded

            for (++iName; iName != pathNames.end(); ++iName) {
                rtengine::RawImage temp(*iName);

                if (!temp.loadRaw(true)) {
                    temp.compress_image(0);     //\ TODO would be better working on original, because is temporary
                    nFiles++;

                    if (ri->getSensorType() == rtengine::ST_BAYER || ri->getSensorType() == rtengine::ST_FUJI_XTRANS) {
                        for (int row = 0; row < H; row++) {
                            for (int col = 0; col < W; col++) {
                                acc[row][col] += temp.data[row][col];
                            }
                        }
                    } else {
                        for (int row = 0; row < H; row++) {
                            for (int col = 0; col < W; col++) {
                                acc[row][3 * col + 0] += temp.data[row][3 * col + 0];
                                acc[row][3 * col + 1] += temp.data[row][3 * col + 1];
                                acc[row][3 * col + 2] += temp.data[row][3 * col + 2];
                            }
                        }
                    }
                }
            }
            const float factor = 1.f / nFiles;
            for (int row = 0; row < H; row++) {
                for (int col = 0; col < rSize; col++) {
                    ri->data[row][col] = acc[row][col] * factor;
                }
            }
        }
    } else {
        ri = new rtengine::RawImage(pathname);

        if (ri->loadRaw(true)) {
            delete ri;
            ri = nullptr;
        } else {
            ri->compress_image(0);
        }
    }
}

void dfInfo::updateBadPixelList(const rtengine::RawImage *df)
{
    if (!df) {
        return;
    }
    constexpr float threshold = 10.f / 8.f;

    if (df->getSensorType() == rtengine::ST_BAYER || df->getSensorType() == rtengine::ST_FUJI_XTRANS) {
        std::vector<rtengine::badPix> badPixelsTemp;

#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            std::vector<rtengine::badPix> badPixelsThread;
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int row = 2; row < df->get_height() - 2; ++row) {
                for (int col = 2; col < df->get_width() - 2; ++col) {
                    const float m = df->data[row - 2][col - 2] + df->data[row - 2][col] + df->data[row - 2][col + 2] +
                                    df->data[row][col - 2] + df->data[row][col + 2] +
                                    df->data[row + 2][col - 2] + df->data[row + 2][col] + df->data[row + 2][col + 2];

                    if (df->data[row][col] > m * threshold) {
                        badPixelsThread.emplace_back(col, row);
                    }
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            badPixelsTemp.insert(badPixelsTemp.end(), badPixelsThread.begin(), badPixelsThread.end());
        }
        badPixels.insert(badPixels.end(), badPixelsTemp.begin(), badPixelsTemp.end());
    } else {
        for (int row = 1; row < df->get_height() - 1; ++row) {
            for (int col = 1; col < df->get_width() - 1; ++col) {
                float m[3];

                for (int c = 0; c < 3; c++) {
                    m[c] = df->data[row - 1][3 * (col - 1) + c] + df->data[row - 1][3 * col + c] + df->data[row - 1][3 * (col + 1) + c] +
                           df->data[row]  [3 * (col - 1) + c] + df->data[row]  [3 * col + c] +
                           df->data[row + 1][3 * (col - 1) + c] + df->data[row + 1][3 * col + c] + df->data[row + 1][3 * (col + 1) + c];
                }

                if (df->data[row][3 * col] > m[0]*threshold || df->data[row][3 * col + 1] > m[1]*threshold || df->data[row][3 * col + 2] > m[2]*threshold) {
                    badPixels.emplace_back(col, row);
                }
            }
        }
    }

    if (rtengine::settings->verbose) {
        std::cout << "Extracted " << badPixels.size() << " pixels from darkframe:" << df->get_filename().c_str() << std::endl;
    }
}


}

class rtengine::DFManager::Implementation final :
    public NonCopyable
{
public:
    void init(const Glib::ustring& pathname);
    Glib::ustring getPathname() const
    {
        return currentPath;
    };
    void getStat(int& totFiles, int& totTemplates) const;
    const RawImage* searchDarkFrame(const std::string& mak, const std::string& mod, int iso, double shut, time_t t);
    const RawImage* searchDarkFrame(const Glib::ustring& filename);
    const std::vector<badPix>* getHotPixels(const std::string& mak, const std::string& mod, int iso, double shut, time_t t);
    const std::vector<badPix>* getHotPixels(const Glib::ustring& filename);
    const std::vector<badPix>* getBadPixels(const std::string& mak, const std::string& mod, const std::string& serial) const;

private:
    typedef std::multimap<std::string, dfInfo> dfList_t;
    typedef std::map<std::string, std::vector<badPix> > bpList_t;
    dfList_t dfList;
    bpList_t bpList;
    Glib::ustring currentPath;
    dfInfo* addFileInfo(const Glib::ustring &filename, bool pool = true);
    dfInfo* find(const std::string &mak, const std::string &mod, int isospeed, double shut, time_t t);
    int scanBadPixelsFile(const Glib::ustring &filename);
};


void rtengine::DFManager::Implementation::init(const Glib::ustring& pathname)
{
    if (pathname.empty()) {
        return;
    }
    std::vector<Glib::ustring> names;

    const auto dir = Gio::File::create_for_path(pathname);
    if (!dir || !dir->query_exists()) {
        return;
    }

    try {

        const auto enumerator = dir->enumerate_children("standard::name");

        while (const auto file = enumerator->next_file()) {
            names.emplace_back(Glib::build_filename(pathname, file->get_name()));
        }

    } catch (Glib::Exception&) {}

    dfList.clear();
    bpList.clear();

    for (const auto &name : names) {
        const auto lastdot = name.find_last_of('.');

        if (lastdot != Glib::ustring::npos && name.substr(lastdot) == ".badpixels") {
            const int n = scanBadPixelsFile(name);

            if (n > 0 && settings->verbose) {
                printf("Loaded %s: %d pixels\n", name.c_str(), n);
            }

            continue;
        }

        try {
            addFileInfo(name);
        } catch(std::exception& e) {}
    }

    // Where multiple shots exist for same group, move filename to list
    for (auto &df : dfList) {
        dfInfo &i = df.second;

        if (!i.pathNames.empty() && !i.pathname.empty()) {
            i.pathNames.push_back(i.pathname);
            i.pathname.clear();
        }

        if (settings->verbose) {
            if (!i.pathname.empty()) {
                printf("%s:  %s\n", i.key().c_str(), i.pathname.c_str());
            } else {
                printf("%s: MEAN of \n    ", i.key().c_str());

                for (std::list<Glib::ustring>::iterator path = i.pathNames.begin(); path != i.pathNames.end(); ++path) {
                    printf("%s, ", path->c_str());
                }

                printf("\n");
            }
        }
    }

    currentPath = pathname;
    return;
}

void rtengine::DFManager::Implementation::getStat(int& totFiles, int& totTemplates) const
{
    totFiles = 0;
    totTemplates = 0;

    for (const auto &df : dfList) {
        const dfInfo &i = df.second;

        if (i.pathname.empty()) {
            totTemplates++;
            totFiles += i.pathNames.size();
        } else {
            totFiles++;
        }
    }
}

/*  The search for the best match is twofold:
 *  if perfect matches for iso and shutter are found, then the list is scanned for lesser distance in time
 *  otherwise if no match is found, the whole list is searched for lesser distance in iso and shutter
 */
const rtengine::RawImage* rtengine::DFManager::Implementation::searchDarkFrame(const std::string& mak, const std::string& mod, int iso, double shut, time_t t)
{
    dfInfo* df = find(toUppercase(mak), toUppercase(mod), iso, shut, t);

    if (df) {
        return df->getRawImage();
    } else {
        return nullptr;
    }
}

const rtengine::RawImage* rtengine::DFManager::Implementation::searchDarkFrame(const Glib::ustring& filename)
{
    for (auto& df : dfList) {
        if (df.second.pathname.compare(filename) == 0) {
            return df.second.getRawImage();
        }
    }

    dfInfo *df = addFileInfo(filename, false);

    if (df) {
        return df->getRawImage();
    }

    return nullptr;
}

const std::vector<rtengine::badPix>* rtengine::DFManager::Implementation::getHotPixels(const Glib::ustring& filename)
{
    for (auto& df : dfList) {
        if (df.second.pathname.compare(filename) == 0) {
            return &df.second.getHotPixels();
        }
    }

    return nullptr;
}

const std::vector<rtengine::badPix>* rtengine::DFManager::Implementation::getHotPixels(const std::string& mak, const std::string& mod, int iso, double shut, time_t t)
{
    dfInfo* df = find(toUppercase(mak), toUppercase(mod), iso, shut, t);

    if (df) {
        if (settings->verbose) {
            if (!df->pathname.empty()) {
                printf("Searched hotpixels from %s\n", df->pathname.c_str());
            } else {
                if (!df->pathNames.empty()) {
                    printf("Searched hotpixels from template (first %s)\n", df->pathNames.begin()->c_str());
                }
            }
        }

        return &df->getHotPixels();
    } else {
        return nullptr;
    }
}

const std::vector<rtengine::badPix>* rtengine::DFManager::Implementation::getBadPixels(const std::string& mak, const std::string& mod, const std::string& serial) const
{
    bpList_t::const_iterator iter;
    bool found = false;

    if (!serial.empty()) {
        // search with serial number first
        std::ostringstream s;
        s << mak << " " << mod << " " << serial;
        iter = bpList.find(s.str());

        if (iter != bpList.end()) {
            found = true;
        }

        if (settings->verbose) {
            if (found) {
                printf("%s.badpixels found\n", s.str().c_str());
            } else {
                printf("%s.badpixels not found\n", s.str().c_str());
            }
        }
    }

    if (!found) {
        // search without serial number
        std::ostringstream s;
        s << mak << " " << mod;
        iter = bpList.find(s.str());

        if (iter != bpList.end()) {
            found = true;
        }

        if (settings->verbose) {
            if (found) {
                printf("%s.badpixels found\n", s.str().c_str());
            } else {
                printf("%s.badpixels not found\n", s.str().c_str());
            }
        }
    }

    if (!found) {
        return nullptr;
    } else {
        return &(iter->second);
    }
}

dfInfo* rtengine::DFManager::Implementation::addFileInfo(const Glib::ustring& filename, bool pool)
{
    const auto ext = getFileExtension(filename);

    if (ext.empty() || !options.is_extention_enabled(ext)) {
        return nullptr;
    }

    auto file = Gio::File::create_for_path(filename);

    if (!file) {
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

        RawImage ri(filename);

        if (ri.loadRaw(false) != 0) { // Read information about shot
            return nullptr;
        }

        if (!pool) {
            const dfInfo n(filename, "", "", 0, 0, 0);
            auto iter = dfList.emplace("", n);
            return &(iter->second);
        }

        FramesData idata(filename);
        /* Files are added in the map, divided by same maker/model,ISO and shutter*/
        std::string key(dfInfo::key(toUppercase(idata.getMake()), toUppercase(idata.getModel()), idata.getISOSpeed(), idata.getShutterSpeed()));
        auto iter = dfList.find(key);

        if (iter == dfList.end()) {
            dfInfo n(filename, toUppercase(idata.getMake()), toUppercase(idata.getModel()), idata.getISOSpeed(), idata.getShutterSpeed(), idata.getDateTimeAsTS());
            iter = dfList.emplace(key, n);
        } else {
            while(iter != dfList.end() && iter->second.key() == key && ABS(iter->second.timestamp - idata.getDateTimeAsTS()) > 60 * 60 * 6) { // 6 hour difference
                ++iter;
            }

            if (iter != dfList.end()) {
                iter->second.pathNames.push_back(filename);
            } else {
                dfInfo n(filename, toUppercase(idata.getMake()), toUppercase(idata.getModel()), idata.getISOSpeed(), idata.getShutterSpeed(), idata.getDateTimeAsTS());
                iter = dfList.emplace(key, n);
            }
        }

        return &(iter->second);

    } catch(Gio::Error&) {}

    return nullptr;
}

dfInfo* rtengine::DFManager::Implementation::find(const std::string& mak, const std::string& mod, int isospeed, double shut, time_t t)
{
    if (dfList.empty()) {
        return nullptr;
    }

    const std::string key(dfInfo::key(mak, mod, isospeed, shut));
    dfList_t::iterator iter = dfList.find(key);

    if (iter != dfList.end()) {
        dfList_t::iterator bestMatch = iter;
        time_t bestDeltaTime = ABS(iter->second.timestamp - t);

        for (++iter; iter != dfList.end() && !key.compare(iter->second.key()); ++iter) {
            const time_t d = ABS(iter->second.timestamp - t);

            if (d < bestDeltaTime) {
                bestMatch = iter;
                bestDeltaTime = d;
            }
        }

        return &(bestMatch->second);
    } else {
        iter = dfList.begin();
        dfList_t::iterator bestMatch = iter;
        double bestD = iter->second.distance(mak, mod, isospeed, shut);

        for (++iter; iter != dfList.end(); ++iter) {
            const double d = iter->second.distance(mak, mod, isospeed, shut);

            if (d < bestD) {
                bestD = d;
                bestMatch = iter;
            }
        }

        return bestD != RT_INFINITY ? &(bestMatch->second) : nullptr ;
    }
}

int rtengine::DFManager::Implementation::scanBadPixelsFile(const Glib::ustring& filename)
{
    FILE *file = ::fopen( filename.c_str(), "r" );

    if (!file) {
        return 0;
    }

    const auto lastdot = filename.find_last_of('.');
    auto dirpos1 = filename.find_last_of('/');
    auto dirpos2 = filename.find_last_of('\\');

    if (dirpos1 == Glib::ustring::npos && dirpos2 == Glib::ustring::npos) {
        dirpos1 = 0;
    } else if (dirpos1 != Glib::ustring::npos && dirpos2 != Glib::ustring::npos) {
        dirpos1 = (dirpos1 > dirpos2 ? dirpos1 : dirpos2);
    } else if (dirpos1 == Glib::ustring::npos) {
        dirpos1 = dirpos2;
    }

    const std::string makmodel(filename, dirpos1 + 1, lastdot - (dirpos1 + 1));
    std::vector<badPix> bp;
    char line[256];

    if (fgets(line, sizeof(line), file)) {
        int x, y;
        int offset = 0;
        int numparms = sscanf(line, "%d %d", &x, &y);

        if (numparms == 1) { // only one number in first line means, that this is the offset.
            offset = x;
        } else if (numparms == 2) {
            bp.emplace_back(x + offset, y + offset);
        }

        while(fgets(line, sizeof(line), file)) {
            if (sscanf(line, "%d %d", &x, &y) == 2) {
                bp.emplace_back(x + offset, y + offset);
            }
        }
    }

    const int numPixels = bp.size();

    if (numPixels > 0) {
        bpList[makmodel] = bp;
    }

    fclose(file);
    return numPixels;
}

rtengine::DFManager& rtengine::DFManager::getInstance()
{
    static DFManager instance;
    return instance;
}

void rtengine::DFManager::init(const Glib::ustring& pathname)
{
    implementation->init(pathname);
}

Glib::ustring rtengine::DFManager::getPathname() const
{
    return implementation->getPathname();
}

void rtengine::DFManager::getStat(int& totFiles, int& totTemplates) const
{
    implementation->getStat(totFiles, totTemplates);
}

const rtengine::RawImage* rtengine::DFManager::searchDarkFrame(const std::string& mak, const std::string& mod, int iso, double shut, time_t t)
{
    return implementation->searchDarkFrame(mak, mod, iso, shut, t);
}

const rtengine::RawImage* rtengine::DFManager::searchDarkFrame(const Glib::ustring& filename)
{
    return implementation->searchDarkFrame(filename);
}

const std::vector<rtengine::badPix>* rtengine::DFManager::getHotPixels(const std::string& mak, const std::string& mod, int iso, double shut, time_t t)
{
    return implementation->getHotPixels(mak, mod, iso, shut, t);
}

const std::vector<rtengine::badPix>* rtengine::DFManager::getHotPixels(const Glib::ustring& filename)
{
    return implementation->getHotPixels(filename);
}

const std::vector<rtengine::badPix>* rtengine::DFManager::getBadPixels(const std::string& mak, const std::string& mod, const std::string& serial) const
{
    return implementation->getBadPixels(mak, mod, serial);
}

rtengine::DFManager::DFManager() :
    implementation(new Implementation)
{
}

rtengine::DFManager::~DFManager() = default;

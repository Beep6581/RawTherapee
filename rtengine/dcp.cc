/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <www.oliverduis.de>
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

#include "dcp.h"
#include <cstring>
#include "safegtk.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "rawimagesource.h"
#include "improcfun.h"

using namespace rtengine;
using namespace rtexif;

#undef CLIP
#define MAXVAL 0xffff
#define CLIP(a) ((a)>0?((a)<MAXVAL?(a):MAXVAL):0)

DCPProfile::DCPProfile(Glib::ustring fname) {
    const int TagColorMatrix1=50721, TagColorMatrix2=50722, TagProfileHueSatMapDims=50937;
    const int TagProfileHueSatMapData1=50938, TagProfileHueSatMapData2=50939;
    const int TagCalibrationIlluminant1=50778, TagCalibrationIlluminant2=50779;
    const int TagProfileLookTableData=50982, TagProfileLookTableDims=50981;  // ProfileLookup is the low quality variant

    aDeltas=NULL; iHueDivisions=iSatDivisions=iValDivisions=iArrayCount=0;

    FILE *pFile = safe_g_fopen(fname, "rb");

    TagDirectory *tagDir=ExifManager::parseTIFF(pFile, false);

    // If there are two profiles, check what is the best target to take
    // We don't mix the profiles as adobe does, since with more and more non-tungsten light
    // it makes no sense. Take the daylight reference light
    Tag* tag = tagDir->getTag(TagCalibrationIlluminant2);
    bool use2nd = (tag!=NULL && tag->toInt(0,SHORT)>=20 && tag->toInt(0,SHORT)<=23);

    // Color Matrix
    tag = tagDir->getTag( use2nd ? TagColorMatrix2 : TagColorMatrix1);

    for (int row=0;row<3;row++) { 
        for (int col=0;col<3;col++) {
            mColorMatrix[col][row]=(float)tag->toDouble((col+row*3)*8);
        }
    }

    // LUT profile? Divisions counts
    bool useSimpleLookup=false;
    tag = tagDir->getTag(TagProfileHueSatMapDims);
    if (tag==NULL) {
        tag=tagDir->getTag(TagProfileLookTableDims);
        useSimpleLookup=true;
    }

    if (tag!=NULL) {
        iHueDivisions=tag->toInt(0); iSatDivisions=tag->toInt(4); iValDivisions=tag->toInt(8);

        // Saturation maps. Need to be unwinded.
        tag = tagDir->getTag(useSimpleLookup ? TagProfileLookTableData : ( use2nd ? TagProfileHueSatMapData2 : TagProfileHueSatMapData1));
        iArrayCount = tag->getCount()/3;

        aDeltas=new HSBModify[iArrayCount];

        const int TIFFFloatSize=4;
        for (int i=0;i<iArrayCount;i++) {
            aDeltas[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
            aDeltas[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
            aDeltas[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
        }
    }

    if (pFile!=NULL) fclose(pFile);
    delete tagDir;

    if (iArrayCount>0) {
        // Convert DNG color matrix to xyz_cam compatible matrix
        int i,j,k;

        double cam_xyz[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                for (k=0; k<3; k++)
                    cam_xyz[i][j] += mColorMatrix[j][k] * (i==k);


        // Multiply out XYZ colorspace 
        double cam_rgb[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (i=0; i < 3; i++)	
            for (j=0; j < 3; j++)
                for (k=0; k < 3; k++) 
                    cam_rgb[i][j] += cam_xyz[i][k] * xyz_sRGB[k][j];

        // Normalize cam_rgb so that:  cam_rgb * (1,1,1) is (1,1,1,1)
        double num;
        for (i=0; i<3; i++) {		
            for (num=j=0; j<3; j++) num += cam_rgb[i][j];
            for (j=0; j<3; j++) cam_rgb[i][j] /= num;
        }

        double rgb_cam[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        RawImageSource::inverse33 (cam_rgb, rgb_cam);

        memset(mXYZCAM,0,sizeof(mXYZCAM));
        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                for (k=0; k<3; k++)
                    mXYZCAM[i][j] += xyz_sRGB[i][k] * rgb_cam[k][j];
    }
}

DCPProfile::~DCPProfile() {
    delete[] aDeltas;
}

void DCPProfile::Apply(Imagefloat *pImg, Glib::ustring workingSpace) const {
    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    if (iArrayCount==0) {
        //===== No LUT- Calculate matrix for direct conversion raw>working space
        double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    mat[i][j] += mWork[i][k] * mXYZCAM[k][j];

        // Apply the matrix part
#pragma omp parallel for
        for (int y=0; y<pImg->height; y++) {
            float newr, newg, newb;
            for (int x=0; x<pImg->width; x++) {
                newr = mat[0][0]*pImg->r[y][x] + mat[0][1]*pImg->g[y][x] + mat[0][2]*pImg->b[y][x];
                newg = mat[1][0]*pImg->r[y][x] + mat[1][1]*pImg->g[y][x] + mat[1][2]*pImg->b[y][x];
                newb = mat[2][0]*pImg->r[y][x] + mat[2][1]*pImg->g[y][x] + mat[2][2]*pImg->b[y][x];

                pImg->r[y][x] = newr; pImg->g[y][x] = newg; pImg->b[y][x] = newb;
            }
        }
    }
    else {
        //===== LUT available- Calculate matrix for conversion raw>ProPhoto
        double m2ProPhoto[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    m2ProPhoto[i][j] += prophoto_xyz[i][k] * mXYZCAM[k][j];

        double m2Work[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i=0; i<3; i++) 
            for (int j=0; j<3; j++)
                for (int k=0; k<3; k++)
                    m2Work[i][j] += mWork[i][k] * xyz_prophoto[k][j];

        // Preperations for LUT
        float hScale = (iHueDivisions < 2) ? 0.0f : (iHueDivisions * (1.0f / 6.0f));
        float sScale = (float) (iSatDivisions - 1);
        float vScale = (float) (iValDivisions - 1);

        int maxHueIndex0 = iHueDivisions - 1;
        int maxSatIndex0 = iSatDivisions - 2;
        int maxValIndex0 = iValDivisions - 2;

        const HSBModify *tableBase = aDeltas;

        int hueStep = iSatDivisions;
        int valStep = iHueDivisions * hueStep;

        // Convert to prophoto and apply LUT
#pragma omp parallel for
        for (int y=0; y<pImg->height; y++) {
            float newr, newg, newb, h,s,v;
            for (int x=0; x<pImg->width; x++) {
                newr = m2ProPhoto[0][0]*pImg->r[y][x] + m2ProPhoto[0][1]*pImg->g[y][x] + m2ProPhoto[0][2]*pImg->b[y][x];
                newg = m2ProPhoto[1][0]*pImg->r[y][x] + m2ProPhoto[1][1]*pImg->g[y][x] + m2ProPhoto[1][2]*pImg->b[y][x];
                newb = m2ProPhoto[2][0]*pImg->r[y][x] + m2ProPhoto[2][1]*pImg->g[y][x] + m2ProPhoto[2][2]*pImg->b[y][x];

                ImProcFunctions::rgb2hsv(newr, newg, newb, h , s, v);
                h*=6.f;  // RT calculates in [0,1]

                // Apply the HueSatMap. Ported from Adobes reference implementation
                float hueShift, satScale, valScale;

                if (iValDivisions < 2)  // Optimize most common case of "2.5D" table.
                {
                    float hScaled = h * hScale;
                    float sScaled = s * sScale;

                    int hIndex0 = (int) hScaled;
                    int sIndex0 = (int) sScaled;

                    sIndex0 = MIN (sIndex0, maxSatIndex0);

                    int hIndex1 = hIndex0 + 1;

                    if (hIndex0 >= maxHueIndex0)
                    {
                        hIndex0 = maxHueIndex0;
                        hIndex1 = 0;
                    }

                    float hFract1 = hScaled - (float) hIndex0;
                    float sFract1 = sScaled - (float) sIndex0;

                    float hFract0 = 1.0f - hFract1;
                    float sFract0 = 1.0f - sFract1;

                    const HSBModify *entry00 = tableBase + hIndex0 * hueStep +
                        sIndex0;

                    const HSBModify *entry01 = entry00 + (hIndex1 - hIndex0) * hueStep;

                    float hueShift0 = hFract0 * entry00->fHueShift +
                        hFract1 * entry01->fHueShift;

                    float satScale0 = hFract0 * entry00->fSatScale +
                        hFract1 * entry01->fSatScale;

                    float valScale0 = hFract0 * entry00->fValScale +
                        hFract1 * entry01->fValScale;

                    entry00++;
                    entry01++;

                    float hueShift1 = hFract0 * entry00->fHueShift +
                        hFract1 * entry01->fHueShift;

                    float satScale1 = hFract0 * entry00->fSatScale +
                        hFract1 * entry01->fSatScale;

                    float valScale1 = hFract0 * entry00->fValScale +
                        hFract1 * entry01->fValScale;

                    hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
                    satScale = sFract0 * satScale0 + sFract1 * satScale1;
                    valScale = sFract0 * valScale0 + sFract1 * valScale1;

                } else {

                    float hScaled = h * hScale;
                    float sScaled = s * sScale;
                    float vScaled = v * vScale;

                    int hIndex0 = (int) hScaled;
                    int sIndex0 = (int) sScaled;
                    int vIndex0 = (int) vScaled;

                    sIndex0 = MIN (sIndex0, maxSatIndex0);
                    vIndex0 = MIN (vIndex0, maxValIndex0);

                    int hIndex1 = hIndex0 + 1;

                    if (hIndex0 >= maxHueIndex0)
                    {
                        hIndex0 = maxHueIndex0;
                        hIndex1 = 0;
                    }

                    float hFract1 = hScaled - (float) hIndex0;
                    float sFract1 = sScaled - (float) sIndex0;
                    float vFract1 = vScaled - (float) vIndex0;

                    float hFract0 = 1.0f - hFract1;
                    float sFract0 = 1.0f - sFract1;
                    float vFract0 = 1.0f - vFract1;

                    const HSBModify *entry00 = tableBase + vIndex0 * valStep + 
                        hIndex0 * hueStep +
                        sIndex0;

                    const HSBModify *entry01 = entry00 + (hIndex1 - hIndex0) * hueStep;

                    const HSBModify *entry10 = entry00 + valStep;
                    const HSBModify *entry11 = entry01 + valStep;

                    float hueShift0 = vFract0 * (hFract0 * entry00->fHueShift +
                        hFract1 * entry01->fHueShift) +
                        vFract1 * (hFract0 * entry10->fHueShift +
                        hFract1 * entry11->fHueShift);

                    float satScale0 = vFract0 * (hFract0 * entry00->fSatScale +
                        hFract1 * entry01->fSatScale) +
                        vFract1 * (hFract0 * entry10->fSatScale +
                        hFract1 * entry11->fSatScale);

                    float valScale0 = vFract0 * (hFract0 * entry00->fValScale +
                        hFract1 * entry01->fValScale) +
                        vFract1 * (hFract0 * entry10->fValScale +
                        hFract1 * entry11->fValScale);

                    entry00++;
                    entry01++;
                    entry10++;
                    entry11++;

                    float hueShift1 = vFract0 * (hFract0 * entry00->fHueShift +
                        hFract1 * entry01->fHueShift) +
                        vFract1 * (hFract0 * entry10->fHueShift +
                        hFract1 * entry11->fHueShift);

                    float satScale1 = vFract0 * (hFract0 * entry00->fSatScale +
                        hFract1 * entry01->fSatScale) +
                        vFract1 * (hFract0 * entry10->fSatScale +
                        hFract1 * entry11->fSatScale);

                    float valScale1 = vFract0 * (hFract0 * entry00->fValScale +
                        hFract1 * entry01->fValScale) +
                        vFract1 * (hFract0 * entry10->fValScale +
                        hFract1 * entry11->fValScale);

                    hueShift = sFract0 * hueShift0 + sFract1 * hueShift1;
                    satScale = sFract0 * satScale0 + sFract1 * satScale1;
                    valScale = sFract0 * valScale0 + sFract1 * valScale1;
                }

                hueShift *= (6.0f / 360.0f);	// Convert to internal hue range.

                h += hueShift; 
                s *= satScale;  // no clipping here, we are RT float :-)
                v *= valScale;

                // RT range correction
                if (h < 0.0f) h += 6.0f;
                if (h >= 6.0f) h -= 6.0f;
                h/=6.f;  
                ImProcFunctions::hsv2rgb( h, s, v, newr, newg, newb);

                pImg->r[y][x] = m2Work[0][0]*newr + m2Work[0][1]*newg + m2Work[0][2]*newb;
                pImg->g[y][x] = m2Work[1][0]*newr + m2Work[1][1]*newg + m2Work[1][2]*newb;
                pImg->b[y][x] = m2Work[2][0]*newr + m2Work[2][1]*newg + m2Work[2][2]*newb;
            }
        }
    }
}


// Integer variant is legacy, only used for thumbs. Simply take the matrix here
void DCPProfile::Apply(Image16 *pImg, Glib::ustring workingSpace) const {
    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    // Calculate matrix for direct conversion raw>working space
    double mat[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int i=0; i<3; i++) 
        for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
                mat[i][j] += mWork[i][k] * mXYZCAM[k][j];

    // Apply the matrix part
#pragma omp parallel for
    for (int y=0; y<pImg->height; y++) {
        float newr, newg, newb;
        for (int x=0; x<pImg->width; x++) {
            newr = mat[0][0]*pImg->r[y][x] + mat[0][1]*pImg->g[y][x] + mat[0][2]*pImg->b[y][x];
            newg = mat[1][0]*pImg->r[y][x] + mat[1][1]*pImg->g[y][x] + mat[1][2]*pImg->b[y][x];
            newb = mat[2][0]*pImg->r[y][x] + mat[2][1]*pImg->g[y][x] + mat[2][2]*pImg->b[y][x];

            pImg->r[y][x] = CLIP((int)newr); pImg->g[y][x] = CLIP((int)newg); pImg->b[y][x] = CLIP((int)newb);
        }
    }
}


// Generates as singleton
DCPStore* DCPStore::getInstance()
{
    static DCPStore* instance_ = 0;
    if ( instance_ == 0 )
    {
        static Glib::Mutex smutex_;
        Glib::Mutex::Lock lock(smutex_);
        if ( instance_ == 0 )
        {
            instance_ = new DCPStore();
        }
    }
    return instance_;
}

// Reads all profiles from the given profiles dir
void DCPStore::init (Glib::ustring rtProfileDir) {
    Glib::Mutex::Lock lock(mtx);

    fileStdProfiles.clear();

    Glib::ustring rootDirName=rtProfileDir;

    if (rootDirName!="") {
        std::deque<Glib::ustring> qDirs;

        qDirs.push_front(rootDirName);

        while (qDirs.size()) {
            // process directory
            Glib::ustring dirname = qDirs.back();
            qDirs.pop_back();

            Glib::Dir* dir = NULL;
            try {
                if (!safe_file_test (dirname, Glib::FILE_TEST_IS_DIR)) return;
                dir = new Glib::Dir (dirname);
            }
            catch (Glib::Exception& fe) {
                return;
            }
            dirname = dirname + "/";
            for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
                Glib::ustring fname = dirname + *i;
                Glib::ustring sname = *i;
                // ignore directories
                if (!safe_file_test (fname, Glib::FILE_TEST_IS_DIR)) {
                    int lastdot = sname.find_last_of ('.');
                    if (lastdot!=Glib::ustring::npos && lastdot<=(int)sname.size()-4 && (!sname.casefold().compare (lastdot, 4, ".dcp"))) {
                        Glib::ustring camShortName = sname.substr(0,lastdot).uppercase();
                        fileStdProfiles[camShortName]=fname;  // they will be loaded and cached on demand
                    }
                } else qDirs.push_front(fname);  // for later scanning
            }
            delete dir;
        }
    }
}

DCPProfile* DCPStore::getProfile (Glib::ustring filename) {
    Glib::Mutex::Lock lock(mtx);

    std::map<Glib::ustring, DCPProfile*>::iterator r = profileCache.find (filename);
    if (r!=profileCache.end()) return r->second;

    // Add profile
    profileCache[filename]=new DCPProfile(filename);

    return profileCache[filename];
}

DCPProfile* DCPStore::getStdProfile(Glib::ustring camShortName) {
    std::map<Glib::ustring, Glib::ustring>::iterator r = fileStdProfiles.find (camShortName.uppercase());
    if (r==fileStdProfiles.end()) return NULL;

    return getProfile(r->second);
}

bool DCPStore::isValidDCPFileName(Glib::ustring filename) const {
    if (!safe_file_test (filename, Glib::FILE_TEST_EXISTS) || safe_file_test (filename, Glib::FILE_TEST_IS_DIR)) return false;
    int pos=filename.find_last_of ('.');
    return pos>0 && !filename.casefold().compare (pos, 4, ".dcp");
}
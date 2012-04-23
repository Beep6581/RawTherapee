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
#include <cstring>

#include "dcp.h"
#include "safegtk.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "rt_math.h"

using namespace std;
using namespace rtengine;
using namespace rtexif;

DCPProfile::DCPProfile(Glib::ustring fname) {
    const int TIFFFloatSize=4;
    const int TagColorMatrix1=50721, TagColorMatrix2=50722, TagProfileHueSatMapDims=50937;
    const int TagProfileHueSatMapData1=50938, TagProfileHueSatMapData2=50939;
    const int TagCalibrationIlluminant1=50778, TagCalibrationIlluminant2=50779;
    const int TagProfileLookTableData=50982, TagProfileLookTableDims=50981;  // ProfileLookup is the low quality variant

    aDeltas1=aDeltas2=NULL; iHueDivisions=iSatDivisions=iValDivisions=iArrayCount=0;

    FILE *pFile = safe_g_fopen(fname, "rb");

    TagDirectory *tagDir=ExifManager::parseTIFF(pFile, false);

    Tag* tag = tagDir->getTag(TagCalibrationIlluminant1); iLightSource1 = (tag!=NULL ? tag->toInt(0,SHORT) : -1);  
    tag = tagDir->getTag(TagCalibrationIlluminant2); iLightSource2 = (tag!=NULL ? tag->toInt(0,SHORT) : -1);

    bool hasSecondHueSat = tagDir->getTag(TagProfileHueSatMapData2)!=NULL;  // some profiles have two matrices, but just one huesat

    // Color Matrix (1 is always there)
    tag = tagDir->getTag(TagColorMatrix1);

    for (int row=0;row<3;row++) { 
        for (int col=0;col<3;col++) {
            mColorMatrix1[col][row]=(float)tag->toDouble((col+row*3)*8);
        }
    }
    ConvertDNGMatrix2XYZCAM(mColorMatrix1,mXYZCAM1);

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
        tag = tagDir->getTag(useSimpleLookup ? TagProfileLookTableData : TagProfileHueSatMapData1);
        iArrayCount = tag->getCount()/3;

        aDeltas1=new HSBModify[iArrayCount];

        for (int i=0;i<iArrayCount;i++) {
            aDeltas1[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
            aDeltas1[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
            aDeltas1[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
        }
    }

    // For second profile, copy everything from first profile is no better data is available
    if (iLightSource2!=-1) {
        // Second matrix
        tag = tagDir->getTag(TagColorMatrix2);

        for (int row=0;row<3;row++) { 
            for (int col=0;col<3;col++) {
                mColorMatrix2[col][row]= (tag!=NULL ? (float)tag->toDouble((col+row*3)*8) : mColorMatrix1[col][row]);
            }
        }

        ConvertDNGMatrix2XYZCAM(mColorMatrix2,mXYZCAM2);

        // Second huesatmap, or copy of first
        if (hasSecondHueSat) {
            aDeltas2=new HSBModify[iArrayCount];

            // Saturation maps. Need to be unwinded.
            tag = tagDir->getTag(TagProfileHueSatMapData2);

            for (int i=0;i<iArrayCount;i++) {
                aDeltas2[i].fHueShift=tag->toDouble((i*3)*TIFFFloatSize);
                aDeltas2[i].fSatScale=tag->toDouble((i*3+1)*TIFFFloatSize);
                aDeltas2[i].fValScale=tag->toDouble((i*3+2)*TIFFFloatSize);
            }
        } else {
            if (aDeltas1!=NULL) {
                aDeltas2=new HSBModify[iArrayCount];
                for (int i=0;i<iArrayCount;i++) aDeltas2[i]=aDeltas1[i];
            }
        }
    }

    if (pFile!=NULL) fclose(pFile);
    delete tagDir;
}

DCPProfile::~DCPProfile() {
    delete[] aDeltas1; delete[] aDeltas2;
}

        // Convert DNG color matrix to xyz_cam compatible matrix
void DCPProfile::ConvertDNGMatrix2XYZCAM(const double (*mColorMatrix)[3], double (*mXYZCAM)[3]) {
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

    for (i=0; i<3; i++)
        for (j=0; j<3; j++) mXYZCAM[i][j]=0;

        for (i=0; i<3; i++)
            for (j=0; j<3; j++)
                for (k=0; k<3; k++)
                    mXYZCAM[i][j] += xyz_sRGB[i][k] * rgb_cam[k][j];
    }


const DCPProfile::HSBModify* DCPProfile::GetBestProfile(DCPLightType preferredProfile, double (*mXYZCAM)[3]) const {
    bool use2=false;

    if (iLightSource2!=-1) {
        DCPLightType t1=GetLightType(iLightSource1); DCPLightType t2=GetLightType(iLightSource2);

        // usually second is the daylight (default if nothing else found)
        if (t2==Daylight) use2=true;

        switch (preferredProfile) {
            case Tungsten:
                if (t1==Tungsten) use2=false; else if (t2==Tungsten) use2=true;
                break;

            case Fluorescent:
                if (t1==Fluorescent) use2=false; else if (t2==Fluorescent) use2=true;
                break;

            case Flash:
                if (t1==Flash) use2=false; else if (t2==Flash) use2=true;
                break;

            default: break;  // e.g. Daylight
        }
}

    // printf("DCP using LightSource %i: %i for requested %i\n", use2?2:1, use2?iLightSource2:iLightSource1, (int)preferredProfile);

    for (int row=0;row<3;row++) { 
        for (int col=0;col<3;col++) {
            mXYZCAM[col][row]= (use2 ? mXYZCAM2[col][row] : mXYZCAM1[col][row]);
        }
}

    return use2?aDeltas2:aDeltas1;
}

DCPLightType DCPProfile::GetLightType(short iLightSource) const {
    if (iLightSource==3 ||  iLightSource==17 || iLightSource==24) return Tungsten;
    if (iLightSource==2 || (iLightSource>=12 && iLightSource<=15)) return Fluorescent;
    if (iLightSource==4) return Flash;
    return Daylight;
}

void DCPProfile::Apply(Imagefloat *pImg, DCPLightType preferredProfile, Glib::ustring workingSpace) const {
    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    double mXYZCAM[3][3]; 
    const HSBModify* tableBase=GetBestProfile(preferredProfile,mXYZCAM);

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

                    int hIndex0 = max((int)hScaled, 0);
                    int sIndex0 = max(min((int)sScaled,maxSatIndex0),0);

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
                    int sIndex0 = max(min((int)sScaled,maxSatIndex0),0);
                    int vIndex0 = max(min((int)vScaled,maxValIndex0),0);

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
void DCPProfile::Apply(Image16 *pImg, DCPLightType preferredProfile, Glib::ustring workingSpace) const {
    TMatrix mWork = iccStore->workingSpaceInverseMatrix (workingSpace);

    double mXYZCAM[3][3]; 
    const HSBModify* tableBase=GetBestProfile(preferredProfile,mXYZCAM);

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

        while (!qDirs.empty()) {
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
                    size_t lastdot = sname.find_last_of ('.');
                    if (lastdot!=Glib::ustring::npos && lastdot<=sname.size()-4 && (!sname.casefold().compare (lastdot, 4, ".dcp"))) {
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
    Glib::ustring name2=camShortName.uppercase();

    // Warning: do NOT use map.find(), since it does not seem to work reliably here
    for (std::map<Glib::ustring, Glib::ustring>::iterator i=fileStdProfiles.begin();i!=fileStdProfiles.end();i++)
        if (name2==(*i).first) return getProfile((*i).second);

    return NULL;
}

bool DCPStore::isValidDCPFileName(Glib::ustring filename) const {
    if (!safe_file_test (filename, Glib::FILE_TEST_EXISTS) || safe_file_test (filename, Glib::FILE_TEST_IS_DIR)) return false;
    size_t pos=filename.find_last_of ('.');
    return pos>0 && !filename.casefold().compare (pos, 4, ".dcp");
}
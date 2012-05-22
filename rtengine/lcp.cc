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

#include "lcp.h"
#include "safegtk.h"
#include "iccmatrices.h"
#include "iccstore.h"
#include "rawimagesource.h"
#include "improcfun.h"
#include "rt_math.h"

#ifdef WIN32
#include <windows.h>
// for GCC32
#ifndef _WIN32_IE
#define _WIN32_IE 0x0600
#endif
#include <shlobj.h>
#endif


using namespace std;
using namespace rtengine;
using namespace rtexif;


LCPModelCommon::LCPModelCommon() {
    focLenX=focLenY=-1; imgXCenter=imgYCenter=0.5;
    for (int i=0;i<5;i++) param[i]=0;
    scaleFac=1;
}

bool LCPModelCommon::empty() const {
    return param[0]==0 && param[1]==0 && param[2]==0;
}

void LCPModelCommon::print() const {
    printf("focLen %g/%g; imgCenter %g/%g; scale %g\n",focLenX,focLenY,imgXCenter,imgYCenter,scaleFac);
    printf("param: %g/%g/%g/%g/%g\n",param[0],param[1],param[2],param[3],param[4]);
}

LCPPersModel::LCPPersModel() {
    focLen=focDist=aperture=0;
}

void LCPPersModel::print() const {
    printf("--- PersModel focLen %g; focDist %g; aperture %g\n", focLen, focDist, aperture);
    printf("Base:\n"); base.print();
    if (!chromRG.empty()) { printf("ChromRG:\n"); chromRG.print(); }
    if (!chromG.empty()) { printf("ChromG:\n"); chromG.print(); }
    if (!chromBG.empty()) { printf("ChromBG:\n"); chromBG.print(); }
    if (!vignette.empty()) { printf("Vignette:\n"); vignette.print(); }
    printf("\n");
}

// if !vignette then geometric
LCPMapper::LCPMapper(LCPProfile* pProf, float focalLength, float focalLength35mm, float aperture, bool vignette, int fullWidth, int fullHeight, const CoarseTransformParams& coarse, int rawRotationDeg)
{
    if (pProf==NULL) return;

    pProf-> calcParams(focalLength, aperture, vignette, mc);

    // determine in what the image with the RAW landscape in comparison (calibration target)
    int rot = (coarse.rotate+rawRotationDeg) % 360;

    bool swapXY  = (rot==90  || rot==270);
    bool mirrorX = (rot==90  || rot==180);
    bool mirrorY = (rot==180 || rot==270);

    // Mention that the Adobe technical paper has a bug here, the DMAX is handled differently for focLen and imgCenter
    int Dmax=fullWidth; if (fullHeight>fullWidth) Dmax=fullHeight;

    // correct focLens
    float focLenX=mc.focLenX; float focLenY=mc.focLenY;
    if (focLenX<0) {  // they may not be given
        // and 35mm may not be given either
        if (focalLength35mm<1) focalLength35mm = focalLength*pProf->sensorFormatFactor;

        focLenX=focLenY=focalLength / ( 35*focalLength/focalLength35mm);  // focLen must be calculated in pixels
    }

    if (swapXY) {
        x0 = (mirrorX ? 1-mc.imgYCenter : mc.imgYCenter) * fullWidth;
        y0 = (mirrorY ? 1-mc.imgXCenter : mc.imgXCenter) * fullHeight;
        fx = focLenY * Dmax;
        fy = focLenX * Dmax;
    } else {
        x0 = (mirrorX ? 1-mc.imgXCenter : mc.imgXCenter) * fullWidth; 
        y0 = (mirrorY ? 1-mc.imgYCenter : mc.imgYCenter) * fullHeight;
        fx = focLenX * Dmax; 
        fy = focLenY * Dmax;
    }

    //printf("\nMapping Dmax=%i f=%g/f35=%g vignette=%i using\n",Dmax, focalLength, focalLength35mm, vignette);
    //mc.print();
}

void LCPMapper::correctDistortion(double& x, double& y) const {
    double xd=(x-x0)/fx, yd=(y-y0)/fy;

    double rsqr      = xd*xd+yd*yd;
    double commonFac = (((mc.param[2]*rsqr + mc.param[1])*rsqr + mc.param[0])*rsqr + 1.)
        + 2. * (mc.param[3] * yd + mc.param[4] * xd);

    double xnew = xd * commonFac + mc.param[4] * rsqr;
    double ynew = yd * commonFac + mc.param[3] * rsqr;

    x = xnew * fx + x0;
    y = ynew * fy + y0;
}

float LCPMapper::correctVignette(int x, int y) const {
    double xd=((double)x-x0)/fx, yd=((double)y-y0)/fy;

    double rsqr      = xd*xd+yd*yd;
    double param0Sqr = mc.param[0]*mc.param[0];

    return 1. + rsqr * (-mc.param[0] + rsqr * ((param0Sqr - mc.param[1])
        - (param0Sqr*mc.param[0] - 2.*mc.param[0]*mc.param[1] + mc.param[2]) *rsqr
        + (param0Sqr*param0Sqr + mc.param[1]*mc.param[1]
    + 2.*mc.param[0]*mc.param[2] - 3.*param0Sqr*mc.param[1]) *rsqr*rsqr));
}

LCPProfile::LCPProfile(Glib::ustring fname) {
    const int BufferSize=8192;
    char buf[BufferSize];

    XML_Parser parser = XML_ParserCreate(NULL);
    if (!parser) throw "Couldn't allocate memory for XML parser";

    XML_SetElementHandler(parser, XmlStartHandler, XmlEndHandler);
    XML_SetCharacterDataHandler(parser, XmlTextHandler);
    XML_SetUserData(parser, (void *)this);


    isFisheye=inCamProfiles=firstLIDone=inPerspect=inAlternateLensID=false;
    sensorFormatFactor=1;
    for (int i=0;i<2000;i++) aPersModel[i]=NULL;
    persModelCount=0;

    FILE *pFile = safe_g_fopen(fname, "rb");

    bool done;
    do {
        int bytesRead = (int)fread(buf, 1, BufferSize, pFile);
        done=feof(pFile);
        if (XML_Parse(parser, buf, bytesRead, done) == XML_STATUS_ERROR) 
            throw "Invalid XML in LCP file";
    } while (!done);

    XML_ParserFree(parser);
}

void LCPProfile::calcParams(float focalLength, float aperture, bool vignette, LCPModelCommon& corr) const {
    // find the frames with the least distance, focal length wise
    LCPPersModel *pLow=NULL, *pHigh=NULL;

    float focalLengthLog=log(focalLength), apertureLog=aperture>0 ? log(aperture) : 0;

    // Pass 1: determining best focal length, if possible different apertures
    for (int pm=0;pm<persModelCount;pm++) {
        float f=aPersModel[pm]->focLen;

        if ((!vignette && !aPersModel[pm]->base.empty()) || (vignette && !aPersModel[pm]->vignette.empty())) {
            if (f <= focalLength && (pLow==NULL || f > pLow->focLen || (f==pLow->focLen && pLow->aperture>aPersModel[pm]->aperture))) {
                pLow=aPersModel[pm];
            }
            if (f >= focalLength && (pHigh==NULL || f < pHigh->focLen || (f==pHigh->focLen && pHigh->aperture<aPersModel[pm]->aperture))) {
                pHigh=aPersModel[pm];
            }
        }
    }

    if (!pLow) 
        pLow=pHigh;
    else if (!pHigh) 
        pHigh=pLow;
    else if (vignette) {
        // We have some, so take the best aperture for vignette (unimportant for distortion)
        float bestFocLenLow=pLow->focLen, bestFocLenHigh=pHigh->focLen;

        for (int pm=0;pm<persModelCount;pm++) {
            float aperLog=log(aPersModel[pm]->aperture);

            if ((!vignette && !aPersModel[pm]->base.empty()) || (vignette && !aPersModel[pm]->vignette.empty())) {
                if (fabs(aPersModel[pm]->focLen-bestFocLenLow)<0.01 && ((aperLog<=apertureLog && pLow->aperture>aperture)
                     || (aperLog<=apertureLog && fabs(apertureLog-aperLog)<fabs(apertureLog - log(pLow->aperture))))) {
                    pLow=aPersModel[pm];
                }
                if (fabs(aPersModel[pm]->focLen-bestFocLenHigh)<0.01 && ((aperLog>=apertureLog && pHigh->aperture<aperture)
                     || (aperLog>=apertureLog && fabs(apertureLog-aperLog)<fabs(apertureLog - log(pHigh->aperture))))) {
                    pHigh=aPersModel[pm];
                }
            }
        }
    }

    if (pLow!=NULL && pHigh!=NULL) {
        // average out the factors, linear interpolation in logarithmic scale
        float facLow=0, facHigh=0;

        // There is as foclen range, take that as basis
        if (pLow->focLen < pHigh->focLen) {
            float diff = log(pHigh->focLen) - log(pLow->focLen);
            facLow  = (log(pHigh->focLen)-focalLengthLog) / diff;
            facHigh = (focalLengthLog-log(pLow->focLen))  / diff;
        } else if (pLow->aperture < aperture && pHigh->aperture > aperture) {
            // FocLen is the same, take aperture (espc. used for vignetting)
            float diff = log(pHigh->aperture) - log(pLow->aperture);
            facLow  = (log(pHigh->aperture)-apertureLog) / diff;
            facHigh = (apertureLog-log(pLow->aperture))  / diff;
        } else {
            facLow=facHigh=0.5;
        }
        
        LCPModelCommon& mLow =(vignette?pLow->vignette :pLow->base);
        LCPModelCommon& mHigh=(vignette?pHigh->vignette:pHigh->base);

        corr.focLenX    = facLow * mLow.focLenX    + facHigh * mHigh.focLenX;
        corr.focLenY    = facLow * mLow.focLenY    + facHigh * mHigh.focLenY;
        corr.imgXCenter = facLow * mLow.imgXCenter + facHigh * mHigh.imgXCenter;
        corr.imgYCenter = facLow * mLow.imgYCenter + facHigh * mHigh.imgYCenter;
        corr.scaleFac   = facLow * mLow.scaleFac   + facHigh * mHigh.scaleFac;

        for (int i=0;i<5;i++) corr.param[i]= facLow * mLow.param[i] + facHigh * mHigh.param[i];

        //printf("LCP ( %i frames) vignette=%i for Fno %g - %g; FocLen %g - %g with fac %g - %g:\n", persModelCount, vignette, pLow->aperture, pHigh->aperture, pLow->focLen, pHigh->focLen, facLow, facHigh); 
    } else printf("Error: LCP file contained no parameters\n");
}

void LCPProfile::print() const {
    printf("=== Profile %s\n", profileName.c_str());
    printf("Frames: %i, RAW: %i; Fisheye: %i; Sensorformat: %f\n",persModelCount,isRaw,isFisheye,sensorFormatFactor);
    for (int pm=0;pm<persModelCount;pm++) aPersModel[pm]->print();
}

void XMLCALL LCPProfile::XmlStartHandler(void *pLCPProfile, const char *el, const char **attr) {
    LCPProfile *pProf=static_cast<LCPProfile*>(pLCPProfile);
    bool parseAttr=false;

    // clean up tagname
    const char* src=strrchr(el,':');
    if (src==NULL) src=const_cast<char*>(el); else src++;

    strcpy(pProf->lastTag,src);

    if (!strcmp("CameraProfiles",src)) pProf->inCamProfiles=true;
    if (!strcmp("AlternateLensIDs",src)) pProf->inAlternateLensID=true;
    if (!pProf->inCamProfiles || pProf->inAlternateLensID) return;

    if (!strcmp("li",src)) {
        pProf->pCurPersModel=new LCPPersModel();
        pProf->pCurCommon=&pProf->pCurPersModel->base;  // iterated to next tags within persModel
        return;
    }

    if (!strcmp("PerspectiveModel",src)) {
        pProf->firstLIDone=true; pProf->inPerspect=true;
        return;
    } else if (!strcmp("FisheyeModel",src)) {
        pProf->firstLIDone=true; pProf->inPerspect=true;
        pProf->isFisheye=true;  // just misses third param, and different path, rest is the same
        return;
    } else if (!strcmp("Description",src)) parseAttr=true;

    // Move pointer to general section
    if (pProf->inPerspect) {
        if (!strcmp("ChromaticRedGreenModel",src)) {
            pProf->pCurCommon=&pProf->pCurPersModel->chromRG;
            parseAttr=true;
        } else if (!strcmp("ChromaticGreenModel",src)) {
            pProf->pCurCommon=&pProf->pCurPersModel->chromG;
            parseAttr=true;
        } else if (!strcmp("ChromaticBlueGreenModel",src)) {
            pProf->pCurCommon=&pProf->pCurPersModel->chromBG;
            parseAttr=true;
        } else if (!strcmp("VignetteModel",src)) {
            pProf->pCurCommon=&pProf->pCurPersModel->vignette;
            parseAttr=true;
        }
    }

    // some profiles (espc. Pentax) have a different structure that is attributes based
    // simulate tags by feeding them in
    if (parseAttr && attr!=NULL) {
        for (int i = 0; attr[i]; i += 2) {
            const char* nameStart=strrchr(attr[i],':');
            if (nameStart==NULL) nameStart=const_cast<char*>(attr[i]); else nameStart++;

            strcpy(pProf->lastTag, nameStart);
            XmlTextHandler(pLCPProfile, attr[i+1], strlen(attr[i+1]));
        }
    }
}

void XMLCALL LCPProfile::XmlTextHandler(void *pLCPProfile, const XML_Char *s, int len) {
    LCPProfile *pProf=static_cast<LCPProfile*>(pLCPProfile);

    if (!pProf->inCamProfiles || pProf->inAlternateLensID) return;

    // Check if it contains non-whitespaces (there are several calls to this for one tag unfortunately)
    bool onlyWhiteSpace=true;
    int i=0;
    while (i<len && onlyWhiteSpace) { onlyWhiteSpace=isspace(s[i]); i++; }
    if (onlyWhiteSpace) return;

    // convert to null terminated
    char raw[len+1]; memcpy(raw,s,len); raw[len]=0;
    char* tag=pProf->lastTag;

    //printf("%s : %s\n",tag,raw);

    // Common data section
    if (!pProf->firstLIDone) {
        // Generic tags are the same for all
        if (!strcmp("ProfileName",tag)) 
            pProf->profileName=Glib::ustring(raw);
        else if (!strcmp("Model",tag)) 
            pProf->camera=Glib::ustring(raw);
        else if (!strcmp("Lens",tag)) 
            pProf->lens=Glib::ustring(raw);
        else if (!strcmp("CameraPrettyName",tag)) 
            pProf->cameraPrettyName=Glib::ustring(raw);
        else if (!strcmp("LensPrettyName",tag)) 
            pProf->lensPrettyName=Glib::ustring(raw);
        else if (!strcmp("CameraRawProfile",tag)) 
            pProf->isRaw=!strcmp("True",raw);
    }

    // --- Now all floating points. Must replace local dot characters
    // WARNING: called by different threads, that may run on different local settings,
    // so don't use system params
    if (atof("1,2345")==1.2345) {
        char* p=raw;
        while (*p) { 
            if (*p=='.') *p=',';
            p++;
        }
    }

    if (!pProf->firstLIDone) {
        if (!strcmp("SensorFormatFactor",tag)) 
            pProf->sensorFormatFactor=atof(raw);   
    }

    // Perspective model base data
    if (!strcmp("FocalLength",tag)) 
        pProf->pCurPersModel->focLen=atof(raw);
    else if (!strcmp("FocusDistance",tag)) 
        pProf->pCurPersModel->focDist=atof(raw);
    else if (!strcmp("ApertureValue",tag)) 
        pProf->pCurPersModel->aperture=atof(raw);

    // Section depended
    if (!strcmp("FocalLengthX",tag)) 
        pProf->pCurCommon->focLenX=atof(raw);
    else if (!strcmp("FocalLengthY",tag)) 
        pProf->pCurCommon->focLenY=atof(raw);
    else if (!strcmp("ImageXCenter",tag)) 
        pProf->pCurCommon->imgXCenter=atof(raw);
    else if (!strcmp("ImageYCenter",tag)) 
        pProf->pCurCommon->imgYCenter=atof(raw);
    else if (!strcmp("ScaleFactor",tag)) 
        pProf->pCurCommon->scaleFac=atof(raw);
    else if (!strcmp("RadialDistortParam1",tag) || !strcmp("VignetteModelParam1",tag)) 
        pProf->pCurCommon->param[0]=atof(raw);
    else if (!strcmp("RadialDistortParam2",tag) || !strcmp("VignetteModelParam2",tag)) 
        pProf->pCurCommon->param[1]=atof(raw);
    else if (!strcmp("RadialDistortParam3",tag) || !strcmp("VignetteModelParam3",tag)) 
        pProf->pCurCommon->param[2]=atof(raw);
    else if (!strcmp("TangentialDistortParam1",tag)) 
        pProf->pCurCommon->param[3]=atof(raw);
    else if (!strcmp("TangentialDistortParam2",tag)) 
        pProf->pCurCommon->param[4]=atof(raw);    
}

void XMLCALL LCPProfile::XmlEndHandler(void *pLCPProfile, const char *el) {
    LCPProfile *pProf=static_cast<LCPProfile*>(pLCPProfile);

    if (strstr(el,":CameraProfiles")) pProf->inCamProfiles=false;
    if (strstr(el,":AlternateLensIDs")) pProf->inAlternateLensID=false;

    if (!pProf->inCamProfiles || pProf->inAlternateLensID) return;

    if (strstr(el,":PerspectiveModel") || strstr(el,":FisheyeModel"))
        pProf->inPerspect=false;
    else if (strstr(el, ":li")) {
        pProf->aPersModel[pProf->persModelCount]=pProf->pCurPersModel;
        pProf->pCurPersModel=NULL;
        pProf->persModelCount++;
    }
}

// Generates as singleton
LCPStore* LCPStore::getInstance()
{
    static LCPStore* instance_ = 0;
    if ( instance_ == 0 )
    {
        static Glib::Mutex smutex_;
        Glib::Mutex::Lock lock(smutex_);
        if ( instance_ == 0 )
        {
            instance_ = new LCPStore();
        }
    }
    return instance_;
}

LCPProfile* LCPStore::getProfile (Glib::ustring filename) {
    if (filename.length()==0 || !isValidLCPFileName(filename)) return NULL;

    Glib::Mutex::Lock lock(mtx);

    std::map<Glib::ustring, LCPProfile*>::iterator r = profileCache.find (filename);
    if (r!=profileCache.end()) return r->second;

    // Add profile (if exists)
    profileCache[filename]=new LCPProfile(filename);
    //profileCache[filename]->print();
    return profileCache[filename];
}

bool LCPStore::isValidLCPFileName(Glib::ustring filename) const {
    if (!safe_file_test (filename, Glib::FILE_TEST_EXISTS) || safe_file_test (filename, Glib::FILE_TEST_IS_DIR)) return false;
    size_t pos=filename.find_last_of ('.');
    return pos>0 && !filename.casefold().compare (pos, 4, ".lcp");
}

Glib::ustring LCPStore::getDefaultCommonDirectory() const {
    Glib::ustring dir;

#ifdef WIN32
    WCHAR pathW[MAX_PATH]={0}; char pathA[MAX_PATH];

    if (SHGetSpecialFolderPathW(NULL,pathW,CSIDL_COMMON_APPDATA,false)) {
        char pathA[MAX_PATH];
        WideCharToMultiByte(CP_UTF8,0,pathW,-1,pathA,MAX_PATH,0,0);
        Glib::ustring fullDir=Glib::ustring(pathA)+Glib::ustring("\\Adobe\\CameraRaw\\LensProfiles\\1.0");
        if (safe_file_test(fullDir, Glib::FILE_TEST_IS_DIR)) dir=fullDir;
    }
#endif

    // TODO: Add Mac paths here

    return dir;
}
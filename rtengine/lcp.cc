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
    return focLenX<0 && focLenY<0;
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
LCPMapper::LCPMapper(LCPProfile* pProf, float focalLength, bool vignette, int fullWidth, int fullHeight, const CoarseTransformParams& coarse, int rawRotationDeg)
{
    if (pProf==NULL) return;

    pProf-> calcBasePerspectiveParams(focalLength, vignette, mc);
    
    // determine in what the image with the RAW landscape in comparison (calibration target)
    int rot = (coarse.rotate+rawRotationDeg) % 360;

    swapXY  = (rot==90  || rot==270);
    mirrorX = (rot==90  || rot==180);
    mirrorY = (rot==180 || rot==270);

    // Mention that the Adobe technical paper has a bug here, the DMAX is handled differently for focLen and imgCenter
    int Dmax=fullWidth; if (fullHeight>fullWidth) Dmax=fullHeight;
    
    if (swapXY) {
        x0 = (mirrorX ? 1-mc.imgYCenter : mc.imgYCenter) * fullWidth;
        y0 = (mirrorY ? 1-mc.imgXCenter : mc.imgXCenter) * fullHeight;
        fx = mc.focLenY * Dmax;
        fy = mc.focLenX * Dmax;
    } else {
        x0 = (mirrorX ? 1-mc.imgXCenter : mc.imgXCenter) * fullWidth; 
        y0 = (mirrorY ? 1-mc.imgYCenter : mc.imgYCenter) * fullHeight;
        fx = mc.focLenX * Dmax; 
        fy = mc.focLenY * Dmax;
    }
}


void LCPMapper::correctDistortion(double& x, double& y) const {
    double xd=(x-x0)/fx, yd=(y-y0)/fy;

    double rsqr      = xd*xd+yd*yd;
    double rsqrrsqr  = rsqr*rsqr; // speed
    double commonFac = (mc.param[2]*rsqrrsqr*rsqr + mc.param[1]*rsqrrsqr + mc.param[0]*rsqr + 1.)
                       + 2 * (mc.param[3] * yd + mc.param[4] * xd);

    double xnew = xd * commonFac + mc.param[4] * rsqr;
    double ynew = yd * commonFac + mc.param[3] * rsqr;

    x = xnew * fx + x0;
    y = ynew * fy + y0;
}

float LCPMapper::correctVignette(double x, double y) const {
    double xd=(x-x0)/fx, yd=(y-y0)/fy;

    double rsqr      = xd*xd+yd*yd;
    double rsqrrsqr  = rsqr*rsqr; // speed
    double param0Sqr = mc.param[0]*mc.param[0];

    return 1. - mc.param[0]*rsqr + (param0Sqr - mc.param[1]) *rsqrrsqr
        - (param0Sqr*mc.param[0] - 2*mc.param[0]*mc.param[1] + mc.param[2]) *rsqrrsqr*rsqr
        + (param0Sqr*param0Sqr + mc.param[1]*mc.param[1]
           + 2*mc.param[0]*mc.param[2] - 3*mc.param[0]*mc.param[0]*mc.param[1]) *rsqrrsqr*rsqrrsqr;
}


LCPProfile::LCPProfile(Glib::ustring fname) {
    const int BufferSize=8192;
    char buf[BufferSize];

    XML_Parser parser = XML_ParserCreate(NULL);
    if (!parser) throw "Couldn't allocate memory for XML parser";

    XML_SetElementHandler(parser, XmlStartHandler, XmlEndHandler);
    XML_SetCharacterDataHandler(parser, XmlTextHandler);
    XML_SetUserData(parser, (void *)this);


    isFisheye=inCamProfiles=firstLIDone=inPerspect=false;

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


void LCPProfile::calcBasePerspectiveParams(float focalLength, bool vignette, LCPModelCommon& corr) {
    // find the frames with the least distance, focal length wise
    LCPPersModel *pLow=NULL, *pHigh=NULL;

    std::list<LCPPersModel*>::const_iterator it;
    for (it=persModels.begin(); it!=persModels.end(); ++it) {
        float f=(*it)->focLen;

        if ((!vignette && !(*it)->base.empty()) || (vignette && !(*it)->vignette.empty())) {
            if (f <= focalLength && (pLow ==NULL || f > pLow->focLen))  pLow= (*it);
            if (f >= focalLength && (pHigh==NULL || f < pHigh->focLen)) pHigh=(*it);
        }
    }

    if (!pLow) pLow=pHigh; if (!pHigh) pHigh=pLow;

    // average out the factors, linear interpolation
    float facLow = 0, facHigh=0;
    if (pLow && pHigh && pLow->focLen < pHigh->focLen) {
        facLow  = (pHigh->focLen-focalLength) / (pHigh->focLen-pLow->focLen);
        facHigh = (focalLength-pLow->focLen) / (pHigh->focLen-pLow->focLen);
    } else {
        facLow=pHigh?0:1; facHigh=pHigh?1:0;
    }

    LCPModelCommon& mLow =(vignette?pLow->vignette :pLow->base);
    LCPModelCommon& mHigh=(vignette?pHigh->vignette:pHigh->base);

    corr.focLenX    = facLow * mLow.focLenX    + facHigh * mHigh.focLenX;
    corr.focLenY    = facLow * mLow.focLenY    + facHigh * mHigh.focLenY;
    corr.imgXCenter = facLow * mLow.imgXCenter + facHigh * mHigh.imgXCenter;
    corr.imgYCenter = facLow * mLow.imgYCenter + facHigh * mHigh.imgYCenter;
    corr.scaleFac   = facLow * mLow.scaleFac   + facHigh * mHigh.scaleFac;

    for (int i=0;i<5;i++) corr.param[i]= facLow * mLow.param[i] + facHigh * mHigh.param[i];
}

void LCPProfile::print() const {
    printf("=== Profile %s\n", profileName.c_str());
    printf("RAW: %i; Fisheye: %i\n",isRaw,isFisheye);
    std::list<LCPPersModel*>::const_iterator it;
    for (it=persModels.begin(); it!=persModels.end(); ++it) (*it)->print();
}

void XMLCALL LCPProfile::XmlStartHandler(void *pLCPProfile, const char *el, const char **attr) {
    LCPProfile *pProf=static_cast<LCPProfile*>(pLCPProfile);

    // clean up tagname
    char* src=strrchr(el,':');
    if (src==NULL) src=const_cast<char*>(el); else src++;

    strcpy(pProf->lastTag,src);

    if (!strcmp("CameraProfiles",src)) pProf->inCamProfiles=true;
    if (!pProf->inCamProfiles) return;

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
    }

    // Move pointer to general section
    if (pProf->inPerspect) {
        if (!strcmp("ChromaticRedGreenModel",src)) 
            pProf->pCurCommon=&pProf->pCurPersModel->chromRG;
        else if (!strcmp("ChromaticGreenModel",src)) 
            pProf->pCurCommon=&pProf->pCurPersModel->chromG;
        else if (!strcmp("ChromaticBlueGreenModel",src)) 
            pProf->pCurCommon=&pProf->pCurPersModel->chromBG;
        else if (!strcmp("VignetteModel",src)) 
            pProf->pCurCommon=&pProf->pCurPersModel->vignette;
    }
}

void XMLCALL LCPProfile::XmlTextHandler(void *pLCPProfile, const XML_Char *s, int len) {
    LCPProfile *pProf=static_cast<LCPProfile*>(pLCPProfile);

    if (!pProf->inCamProfiles) return;

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
    struct lconv * lc=localeconv();

    char* p=raw;
    while (*p) { 
        if (*p=='.') *p=lc->decimal_point[0];
        p++;
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

    if (!pProf->inCamProfiles) return;

    if (strstr(el,":PerspectiveModel") || strstr(el,":FisheyeModel"))
        pProf->inPerspect=false;
    else if (strstr(el, ":li")) {
        pProf->persModels.push_back(pProf->pCurPersModel);
        pProf->pCurPersModel=NULL;
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
    if (filename.length()==0) return NULL;

    Glib::Mutex::Lock lock(mtx);

    std::map<Glib::ustring, LCPProfile*>::iterator r = profileCache.find (filename);
    if (r!=profileCache.end()) return r->second;

    // Add profile
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
#else
    printf("Sorry, default LCP directory are currently only configured on Windows\n");
#endif

    return dir;
}
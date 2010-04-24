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
#include <iccstore.h>
#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif
#include <iccmatrices.h>
#include <glib/gstdio.h>

namespace rtengine {

ICCStore iccStore;

const double (*wprofiles[])[3]  = {sRGB_d50, adobe_d50, prophoto_d50, widegamut_d50, bruce_d50, beta_d50, best_d50};
const double (*iwprofiles[])[3] = {d50_sRGB, d50_adobe, d50_prophoto, d50_widegamut, d50_bruce, d50_beta, d50_best};
const char* wpnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB"};         

std::vector<std::string> getWorkingProfiles () {

    std::vector<std::string> res;
    for (int i=0; i<sizeof(wpnames)/sizeof(wpnames[0]); i++) 
        res.push_back (wpnames[i]);
    return res;
}

std::vector<std::string> getOutputProfiles () {

    return iccStore.getOutputProfiles ();
}

std::vector<std::string> ICCStore::getOutputProfiles () {

    std::vector<std::string> res;
    for (std::map<std::string, cmsHPROFILE>::iterator i=fileProfiles.begin(); i!=fileProfiles.end(); i++)
        res.push_back (i->first);
    return res;
}


ICCStore::ICCStore () {

    cmsErrorAction (LCMS_ERROR_SHOW);

    int N = sizeof(wpnames)/sizeof(wpnames[0]);
    for (int i=0; i<N; i++) {
        wProfiles[wpnames[i]] = iccStore.createFromMatrix (wprofiles[i]);
        wProfilesGamma[wpnames[i]] = iccStore.createFromMatrix (wprofiles[i], true);
        wMatrices[wpnames[i]] = wprofiles[i];
        iwMatrices[wpnames[i]] = iwprofiles[i];
    }
    
    double mat[3][3]={1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0};
    xyz  = createFromMatrix (mat, false, "XYZ");
    srgb = cmsCreate_sRGBProfile ();
}

int ICCStore::numOfWProfiles () {

    return sizeof(wpnames)/sizeof(wpnames[0]);
}

TMatrix ICCStore::workingSpaceMatrix (Glib::ustring name) {

    std::map<std::string, TMatrix>::iterator r = wMatrices.find (name);
    if (r!=wMatrices.end()) 
        return r->second;
    else
        return wMatrices["sRGB"];
}

TMatrix ICCStore::workingSpaceInverseMatrix (Glib::ustring name) {

    std::map<std::string, TMatrix>::iterator r = iwMatrices.find (name);
    if (r!=iwMatrices.end()) 
        return r->second;
    else 
        return iwMatrices["sRGB"];
}

cmsHPROFILE ICCStore::workingSpace (Glib::ustring name) {

    std::map<std::string, cmsHPROFILE>::iterator r = wProfiles.find (name);
    if (r!=wProfiles.end()) 
        return r->second;
    else
        return wProfiles["sRGB"];
}

cmsHPROFILE ICCStore::workingSpaceGamma (Glib::ustring name) {

    std::map<std::string, cmsHPROFILE>::iterator r = wProfilesGamma.find (name);
    if (r!=wProfilesGamma.end()) 
        return r->second;
    else
        return wProfilesGamma["sRGB"];
}

cmsHPROFILE ICCStore::getProfile (Glib::ustring name) {


    std::map<std::string, cmsHPROFILE>::iterator r = fileProfiles.find (name);
    if (r!=fileProfiles.end()) 
        return r->second;
    else {
        if (!name.compare (0, 5, "file:") && Glib::file_test (name.substr(5), Glib::FILE_TEST_EXISTS) && !Glib::file_test (name.substr(5), Glib::FILE_TEST_IS_DIR)) {
            ProfileContent pc (name.substr(5));
            if (pc.data) {
                cmsHPROFILE profile = pc.toProfile ();
                if (profile) {
                    fileProfiles[name] = profile;
                    fileProfileContents[name] = pc;
                    return profile;
                }
            }
        }
    }
    return NULL;
}

ProfileContent ICCStore::getContent (Glib::ustring name) {

    return fileProfileContents[name];
}

std::vector<std::string> ICCStore::parseDir (Glib::ustring pdir) {

    fileProfiles.clear ();
    fileProfileContents.clear ();
    std::vector<std::string> result;
    if (pdir!="") {
        // process directory
        Glib::ustring dirname = pdir;
        Glib::Dir* dir = NULL;
        try {
	    if (!Glib::file_test (dirname, Glib::FILE_TEST_IS_DIR))
                return result;
            dir = new Glib::Dir (dirname);
        }
        catch (Glib::Exception& fe) {
            return result;
        }
        dirname = dirname + "/";
        for (Glib::DirIterator i = dir->begin(); i!=dir->end(); ++i) {
            Glib::ustring fname = dirname + *i;
            Glib::ustring sname = *i;
            // ignore directories
            if (!Glib::file_test (fname, Glib::FILE_TEST_IS_DIR)) {
                int lastdot = sname.find_last_of ('.');
                if (lastdot!=Glib::ustring::npos && lastdot<=sname.size()-4 && (!sname.casefold().compare (lastdot, 4, ".icm") || !sname.casefold().compare (lastdot, 4, ".icc"))) {
//                    printf ("processing file %s...\n", fname.c_str());
                    Glib::ustring name = sname.substr(0,lastdot);
                    ProfileContent pc (fname);
                    if (pc.data) {
                        cmsHPROFILE profile = pc.toProfile ();
                        if (profile) {
                            fileProfiles[name] = profile;
                            fileProfileContents[name] = pc;
                            result.push_back (name);
                        }
                    }
                }
            }
        }
        delete dir;
    }
    return result;
}

ProfileContent::ProfileContent (Glib::ustring fileName) {

    data = NULL;
    FILE* f = g_fopen (fileName.c_str(), "rb");
    if (!f)
        return;
    fseek (f, 0, SEEK_END);
    length = ftell (f);
    fseek (f, 0, SEEK_SET);
    data = new char[length+1];
    fread (data, length, 1, f);
	data[length] = 0;
    fclose (f);
}

ProfileContent::ProfileContent (const ProfileContent& other) {

    length = other.length;
    if (other.data) {
        data = new char[length+1];
        memcpy (data, other.data, length+1);
    }
    else
        data = NULL;
}

ProfileContent& ProfileContent::operator= (const ProfileContent other) {

    length = other.length;
    if (data)
        delete [] data;
    if (other.data) {
        data = new char[length+1];
        memcpy (data, other.data, length+1);
    }
    else
        data = NULL;
    return *this;
}

ProfileContent::~ProfileContent () {

    if (data)
        delete [] data;
}
        
cmsHPROFILE ProfileContent::toProfile () {

    if (data)
        return cmsOpenProfileFromMem (data, length);
    else
        return NULL;
}

cmsHPROFILE ICCStore::createFromMatrix (const double matrix[3][3], bool gamma, Glib::ustring name) {

    static const unsigned phead[] =
        { 1024, 0, 0x2100000, 0x6d6e7472, 0x52474220, 0x58595a20, 0, 0, 0,
          0x61637370, 0, 0, 0, 0, 0, 0, 0, 0xf6d6, 0x10000, 0xd32d };
    unsigned pbody[] =
      { 10, 0x63707274, 0, 36,	/* cprt */
	    0x64657363, 0, 40,	/* desc */
    	0x77747074, 0, 20,	/* wtpt */
	    0x626b7074, 0, 20,	/* bkpt */
	    0x72545243, 0, 14,	/* rTRC */
	    0x67545243, 0, 14,	/* gTRC */
	    0x62545243, 0, 14,	/* bTRC */
	    0x7258595a, 0, 20,	/* rXYZ */
	    0x6758595a, 0, 20,	/* gXYZ */
	    0x6258595a, 0, 20 };	/* bXYZ */
    static const unsigned pwhite[] = { 0xf351, 0x10000, 0x116cc };
    // 0x63757276 : curveType, 0 : reserved, 1 : entries (1=gamma, 0=identity), 0x1000000=1.0 
    unsigned pcurve[] = { 0x63757276, 0, 0, 0x1000000 };
//    unsigned pcurve[] = { 0x63757276, 0, 1, 0x1000000 };

	if (gamma) {
        pcurve[2] = 1;
        pcurve[3] = 0x1f00000;
    }

    // constructing profile header
    unsigned* oprof = new unsigned [phead[0]/sizeof(unsigned)];
    memset (oprof, 0, phead[0]);
    memcpy (oprof, phead, sizeof(phead));

    oprof[0] = 132 + 12*pbody[0];
    
    // constructing tag directory (pointers inside the file), and types
    // 0x74657874 : text
    // 0x64657363 : description tag
    for (int i=0; i < pbody[0]; i++) {
      oprof[oprof[0]/4] = i ? (i > 1 ? 0x58595a20 : 0x64657363) : 0x74657874;
      pbody[i*3+2] = oprof[0];
      oprof[0] += (pbody[i*3+3] + 3) & -4;
    }   
    memcpy (oprof+32, pbody, sizeof(pbody));

    // wtpt
    memcpy ((char *)oprof+pbody[8]+8, pwhite, sizeof(pwhite));

    // r/g/b TRC
    for (int i=4; i < 7; i++)
      memcpy ((char *)oprof+pbody[i*3+2], pcurve, sizeof(pcurve));
      
    // r/g/b XYZ
//    pseudoinverse ((double (*)[3]) out_rgb[output_color-1], inverse, 3);
    for (int i=0; i < 3; i++)
      for (int j=0; j < 3; j++) {
        oprof[pbody[j*3+23]/4+i+2] = matrix[j][i] * 0x10000 + 0.5;
//    	for (num = k=0; k < 3; k++)
//    	  num += xyzd50_srgb[i][k] * inverse[j][k];
      }

    // convert to network byte order
    for (int i=0; i < phead[0]/4; i++)
      oprof[i] = htonl(oprof[i]);
      
    // cprt
    strcpy ((char *)oprof+pbody[2]+8, "--rawtherapee profile--");

    // desc
    oprof[pbody[5]/4+2] = name.size() + 1;   
    strcpy ((char *)oprof+pbody[5]+12, name.c_str());


    return cmsOpenProfileFromMem (oprof, ntohl(oprof[0]));
}
}

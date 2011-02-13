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
#ifdef QTBUILD
#include <QDir>
#include <QFile>
#include <QFileInfo>
#else
#include <glib/gstdio.h>
#endif
#include <string>
#include "string.h"

namespace rtengine {

ICCStore* iccStore;

// Members belonging to ProfileContent
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ProfileContent::ProfileContent (const String& fileName) {

#ifdef QTBUILD
	QFile f (fileName);
	if (f.open (QFile::ReadOnly)) {
		QDataStream in (&f);
		length = QFileInfo(f).size();
		data = new unsigned char[length+1];
		in.readRawData ((char*)data, length);
		data[length] = 0;
		f.close ();
	}
#else
    data = NULL;
    FILE* f = g_fopen (fileName.c_str(), "rb");
    if (!f)
        return;
    fseek (f, 0, SEEK_END);
    length = ftell (f);
    fseek (f, 0, SEEK_SET);
    data = new unsigned char[length+1];
    fread (data, length, 1, f);
	data[length] = 0;
    fclose (f);
#endif
}

ProfileContent::ProfileContent (const ProfileContent& other) {

    length = other.length;
    if (other.data) {
        data = new unsigned char[length+1];
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
        data = new unsigned char[length+1];
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

// Members belonging to ICCStore
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const Matrix33 sRGB_d50 (0.435859, 0.385336, 0.143023, 0.222385, 0.717021, 0.0605936, 0.0139162, 0.0971389, 0.713817);
const Matrix33 adobe_d50 (0.6097395054954, 0.2052518325737, 0.1492308013399, 0.3111142944042, 0.6256618480686, 0.0632241329247, 0.0194773131652, 0.0608872306106, 0.744846530711);
const Matrix33 prophoto_d50 (0.797675, 0.135192, 0.0313534, 0.288040, 0.711874, 0.000086, 0.000000, 0.000000, 0.825210);
const Matrix33 widegamut_d50 (0.716105, 0.100930, 0.147186, 0.258187, 0.724938, 0.0168748, 0.000000, 0.0517813, 0.773429);
const Matrix33 bruce_d50 (0.4941607255908, 0.3204990468435, 0.1495612990809, 0.2521412970174, 0.684494580042, 0.0633643619597, 0.0157852934504, 0.062927176507, 0.746498914581);
const Matrix33 beta_d50 (0.671254, 0.174583, 0.118383, 0.303273, 0.663786, 0.0329413, 0.000000, 0.040701, 0.784509);
const Matrix33 best_d50 (0.632670, 0.204556, 0.126995, 0.228457, 0.737352, 0.0341908, 0.000000, 0.00951424, 0.815696);


Matrix33 wprofiles[] = {sRGB_d50, adobe_d50, prophoto_d50, widegamut_d50, bruce_d50, beta_d50, best_d50};
const char* wpnames[] = {"sRGB", "Adobe RGB", "ProPhoto", "WideGamut", "BruceRGB", "Beta RGB", "BestRGB"};

ICCStore::ICCStore () {

    int N = sizeof(wpnames)/sizeof(wpnames[0]);

    for (int i=0; i<N; i++) {
        wMatrices[wpnames[i]] = wprofiles[i];
        wProfiles[wpnames[i]] = createFromMatrix (wprofiles[i]);
        iwMatrices[wpnames[i]] = wprofiles[i].inverse ();
        wProfilesGamma[wpnames[i]] = createFromMatrix (wprofiles[i].inverse (), true);
    }
    float mat[3][3]={1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0};
    xyz  = createFromMatrix (mat, false, "XYZ");
    srgb = cmsCreate_sRGBProfile ();
}

StringList ICCStore::getWorkingProfiles () {

    StringList res;
    for (int i=0; i<sizeof(wpnames)/sizeof(wpnames[0]); i++)
        res.push_back (wpnames[i]);
    return res;
}

StringList ICCStore::getOutputProfiles () {

	StringList res;
    for (std::map<String, cmsHPROFILE>::iterator i=fileProfiles.begin(); i!=fileProfiles.end(); i++)
        res.push_back (i->first);
    return res;
}

Matrix33 ICCStore::workingSpaceMatrix (const String& name) {

    std::map<String, Matrix33>::iterator r = wMatrices.find (name);
    if (r!=wMatrices.end()) 
        return r->second;
    else
        return wMatrices["sRGB"];
}

Matrix33 ICCStore::workingSpaceInverseMatrix (const String& name) {

    std::map<String, Matrix33>::iterator r = iwMatrices.find (name);
    if (r!=iwMatrices.end()) 
        return r->second;
    else 
        return iwMatrices["sRGB"];
}

cmsHPROFILE ICCStore::workingSpace (const String& name) {

    std::map<String, cmsHPROFILE>::iterator r = wProfiles.find (name);
    if (r!=wProfiles.end()) 
        return r->second;
    else
        return wProfiles["sRGB"];
}

cmsHPROFILE ICCStore::workingSpaceGamma (const String& name) {

    std::map<String, cmsHPROFILE>::iterator r = wProfilesGamma.find (name);
    if (r!=wProfilesGamma.end()) 
        return r->second;
    else
        return wProfilesGamma["sRGB"];
}

cmsHPROFILE ICCStore::getProfile (const String& name) {

    std::map<String, cmsHPROFILE>::iterator r = fileProfiles.find (name);
    if (r!=fileProfiles.end()) 
        return r->second;
    else {
#ifdef QTBUILD
    	QFileInfo finfo;
    	if (name.left(5) == "file:" && (finfo=QFileInfo (name.mid(5))).exists() && !finfo.isDir() && finfo.isReadable()) {
        	ProfileContent pc (name.mid(5));
#else
        if (!name.compare (0, 5, "file:") && Glib::file_test (name.substr(5), Glib::FILE_TEST_EXISTS) && !Glib::file_test (name.substr(5), Glib::FILE_TEST_IS_DIR)) {
        	ProfileContent pc (name.substr(5));
#endif
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

ProfileContent ICCStore::getContent (const String& name) {

    return fileProfileContents[name];
}

StringList ICCStore::parseDir (const String& pdir) {

    fileProfiles.clear ();
    fileProfileContents.clear ();
    StringList result;
    if (pdir!="") {
#ifdef QTBUILD
    	QDir dir (pdir);
    	if (!dir.exists())
    		return result;

        QStringList filters;
        filters << "*.icm" << "*.icc";

        QFileInfoList entries = dir.entryInfoList (filters, QDir::Files | QDir::Readable | QDir::NoDotAndDotDot);
    	foreach (const QFileInfo& fi, entries) {
    		ProfileContent pc (dir.absoluteFilePath (fi.absoluteFilePath ()));
			if (pc.data) {
				cmsHPROFILE profile = pc.toProfile ();
				if (profile) {
					QString name = fi.completeBaseName ();
					fileProfiles[name] = profile;
					fileProfileContents[name] = pc;
					result.push_back (name);
				}
			}
    	}
#else
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
                    String name = sname.substr(0,lastdot);
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
#endif
    }
    return result;
}


cmsHPROFILE ICCStore::createFromMatrix (const Matrix33& matrix, bool gamma, const String& name) {

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
        oprof[pbody[j*3+23]/4+i+2] = matrix.data[i][j] * 0x10000 + 0.5;
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
    strcpy ((char *)oprof+pbody[5]+12, String2PChar(name));


    return cmsOpenProfileFromMem (oprof, ntohl(oprof[0]));
}

StringList getWorkingProfiles () {

    return getWorkingProfiles ();
}

StringList getOutputProfiles () {

    return getOutputProfiles ();
}

}

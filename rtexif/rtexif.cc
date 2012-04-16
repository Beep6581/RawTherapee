/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Some parts of the source code (e.g. ciff support) are taken from dcraw
 *  that is copyrighted by Dave Coffin
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
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <stdint.h>

#include "rtexif.h"

using namespace std;

namespace rtexif {

StdInterpreter stdInterpreter;

//--------------- class TagDirectory ------------------------------------------
// this class is a collection (an array) of tags
//-----------------------------------------------------------------------------

#define TAG_SUBFILETYPE 0x00fe

TagDirectory::TagDirectory () 
  : attribs(ifdAttribs), order(INTEL), parent(NULL) {}

TagDirectory::TagDirectory (TagDirectory* p, const TagAttrib* ta, ByteOrder border) 
    : attribs(ta), order(border), parent(p) {}

TagDirectory::TagDirectory (TagDirectory* p, FILE* f, int base, const TagAttrib* ta, ByteOrder border, bool skipIgnored) {

  attribs = ta;
  order = border;
  parent = p;
  
  int numOfTags = get2 (f, order);
  if (numOfTags<=0 || numOfTags>200)
    return;
    
  bool thumbdescr = false;
  for (int i=0; i<numOfTags; i++) {

    Tag* newTag = new Tag (this, f, base);

    // filter out tags with unknown type
    if ((int)newTag->getType()==0) {
        delete newTag;
        continue;
    }
        
    if (skipIgnored) {
    int id = newTag->getID();

    // detect and possibly ignore tags of directories belonging to the embedded thumbnail image
    if (attribs==ifdAttribs && id==TAG_SUBFILETYPE && newTag->toInt()!=0)
        thumbdescr = true;

    const TagAttrib* attrib = getAttrib (id);

    if (!attrib || attrib->ignore==1 || (thumbdescr && attrib->ignore==2)) 
      delete newTag;
    else 
        addTag (newTag);
    } else addTag (newTag);
  }  
}
 
TagDirectory::~TagDirectory () {

  for (int i=0; i<tags.size(); i++) 
    delete tags[i];
}

class CompareTags {
public:
    int operator() (Tag* const& a, Tag* const& b) const {
        return a->getID() < b->getID();
    }
};

void TagDirectory::sort () {

    std::sort (tags.begin(), tags.end(), CompareTags());
    for (int i=0; i<tags.size(); i++) 
        if (tags[i]->isDirectory())
            for (int j=0; tags[i]->getDirectory(j); j++)
                tags[i]->getDirectory(j)->sort ();
}
TagDirectory*  TagDirectory::getRoot()
{
	if(parent) return parent->getRoot();
	else return this;
}

const TagAttrib* TagDirectory::getAttrib (int id) {

    if (attribs)
        for (int i=0; attribs[i].ignore!=-1; i++)
            if (attribs[i].ID==id)
                return &attribs[i];

    return NULL;
}

const TagAttrib* TagDirectory::getAttrib (const char* name) {

    if (attribs)
        for (int i=0; attribs[i].ignore!=-1; i++)
            if (!strcmp (attribs[i].name, name))
                return &attribs[i];

    return NULL;
}

void TagDirectory::printAll () const {
  
  for (int i=0; i<tags.size(); i++) {
    std::string name = tags[i]->nameToString ();
    if (tags[i]->isDirectory())
        for (int j=0; tags[i]->getDirectory(j); j++) {
            printf ("==== DIRECTORY %s[%d]: ====\n", name.c_str(), j);
            tags[i]->getDirectory(j)->printAll ();
            printf ("==== END OF DIRECTORY %s[%d] ====\n", name.c_str(), j);
        }
    else { 
        std::string value = tags[i]->valueToString ();
        printf ("%s: %s\n", name.c_str(), value.c_str());
     }
  }
}

void TagDirectory::addTag (Tag* tag) {

  // look up if it already exists:
  if (getTag (tag->getID()))
    delete tag;
  else
    tags.push_back (tag);
}

void TagDirectory::addTagFront (Tag* tag) {

  // look up if it already exists:
  if (getTag (tag->getID()))
    delete tag;
  else
    tags.insert (tags.begin(), tag);
}

void TagDirectory::replaceTag (Tag* tag) {

  // look up if it already exists:
  for (int i=0; i<tags.size(); i++) 
    if (tags[i]->getID()==tag->getID()) {
      delete tags[i];
      tags[i] = tag;
      return;
    }
  tags.push_back (tag);
}

Tag* TagDirectory::getTag (int ID) const {

  for (int i=0; i<tags.size(); i++) 
    if (tags[i]->getID()==ID)
      return tags[i];
  return NULL;
}

Tag* TagDirectory::getTag (const char* name) const {

  if (attribs) {
    for (int i=0; attribs[i].ignore!=-1; i++)
        if (!strcmp (attribs[i].name, name))
            return getTag (attribs[i].ID);
  }
  return NULL;
}

Tag* TagDirectory::findTag (const char* name) const {
	  if (attribs) {
		for (int i=0; attribs[i].ignore!=-1; i++)
			if (!strcmp (attribs[i].name, name)){
				Tag* t= getTag (attribs[i].ID);
				if(t) return t;
				else break;
			}
	  }
	  for (int i=0; i<tags.size(); i++)
	     if(tags[i]->isDirectory()){
	    	 TagDirectory *dir = tags[i]->getDirectory();
	    	 Tag* t=dir->findTag(name);
	    	 if(t) return t;
	     }
	  return NULL;
}

void TagDirectory::keepTag (int ID) {
    for (int i=0; i<tags.size(); i++) 
        if (tags[i]->getID()==ID) tags[i]->setKeep(true);
}

int TagDirectory::calculateSize () {

  int size = 2; // space to store the number of tags
  for (int i=0; i<tags.size(); i++) 
    if (tags[i]->getKeep()) 
        size += 12 + tags[i]->calculateSize ();
    
  size += 4; // next ifd pointer
  return size;
}
 
TagDirectory* TagDirectory::clone (TagDirectory* parent) {
    
    TagDirectory* td = new TagDirectory (parent, attribs, order);
    for (int i=0; i<tags.size(); i++)
        td->tags.push_back (tags[i]->clone (td));
    return td;
}

int TagDirectory::write (int start, unsigned char* buffer) {

    int size = calculateSize ();
    int tagnum = 0;
    int nondirspace = 0;
    for (int i=0; i<tags.size(); i++)
        if (tags[i]->getKeep()) {
            tagnum++;
            if (!tags[i]->isDirectory())
                nondirspace += tags[i]->calculateSize();
        }
    int nextValOffs = start + 2 + tagnum * 12 + 4;
    int nextDirOffs = nextValOffs + nondirspace;
    int pos = start;
    sset2 (tagnum, buffer+start, order);
    pos += 2;
    int maxPos = start + size;
    for (int i=0; i<tags.size(); i++) {
        if (tags[i]->getKeep()) {
            if (!tags[i]->isDirectory())
                nextValOffs = tags[i]->write (pos, nextValOffs, buffer);  // pos: where to put the tag, dataoffset: the place where the value can be put. return: next data offset
            else
                nextDirOffs = tags[i]->write (pos, nextDirOffs, buffer);  // pos: where to put the tag, dataoffset: the place where the value can be put. return: next data offset
            
            pos += 12;
        }
    }
    sset4 (0, buffer+pos, order);
    return maxPos;
}

void TagDirectory::applyChange (std::string name, std::string value) {

    std::string::size_type dp = name.find_first_of ('.');
    std::string fseg = name.substr (0,dp);
    // this is a final segment: apply change
    if (dp==std::string::npos) {

        Tag* t = NULL;
        for (int i=0; i<tags.size(); i++)
            if (tags[i]->nameToString()==fseg) {
                t = tags[i];
                break;
            }
        
        if (value=="#keep" && t)
            t->setKeep (true);
        else if (value=="#delete" && t) 
            t->setKeep (false);
        else if (t && !t->isDirectory()) 
            t->valueFromString (value);
        else {
            const TagAttrib* attrib = NULL;
            for (int i=0; attribs[i].ignore!=-1; i++)
                if (!strcmp (attribs[i].name, fseg.c_str())) {
                    attrib = &attribs[i];
                    break;
                }            
            if (attrib) {
                Tag* nt = new Tag (this, attrib);
                nt->initString (value.c_str());
                addTag (nt);
            }
        }
    }
    // this is a subdirectory
    else {
        // try to find it
        std::string::size_type dp1 = fseg.find_first_of ('[');
        std::string basename = fseg.substr (0,dp1);
        Tag* t = NULL;
        int dirnum = -1;
        for (int i=0; i<tags.size(); i++)
            if (tags[i]->isDirectory()) {
                for (int j=0; tags[i]->getDirectory(j); j++) {
                    if (tags[i]->nameToString(j) == fseg) {
                        t = tags[i];
                        dirnum = j;
                        break;
                    }
                }
                if (!t && tags[i]->nameToString() == basename) { // found it, but that directory index does not exist
                    t = tags[i];
                    dirnum = -1;
                }
            }
        if (!t && value!="#keep" && value!="#delete") {
            const TagAttrib* attrib = NULL;
            for (int i=0; attribs[i].ignore!=-1; i++)
                if (!strcmp (attribs[i].name, fseg.c_str())) {
                    attrib = &attribs[i];
                    break;
                }            
            if (attrib && attrib->subdirAttribs) {
                t = new Tag (this, attrib);
                t->initSubDir ();
                addTag (t);
            }
            dirnum = 0;
        }
        if (t && dirnum>=0)
            t->getDirectory(dirnum)->applyChange (name.substr (dp+1, std::string::npos), value);
    }
}

TagDirectoryTable::TagDirectoryTable ()
:zeroOffset(0),valuesSize(0)
{
}

TagDirectoryTable::TagDirectoryTable (TagDirectory* p, unsigned char *v,int memsize,int offs, TagType type, const TagAttrib* ta, ByteOrder border)
:TagDirectory(p,ta,border),zeroOffset(offs),valuesSize(memsize),defaultType( type )
{
	  values = new unsigned char[valuesSize];
	  memcpy(values,v,valuesSize);
	  for( const TagAttrib* tattr = ta; tattr->ignore != -1; tattr++){
		  Tag* newTag = new Tag (this, tattr, (values + zeroOffset+ tattr->ID*getTypeSize(type)),type);
		  tags.push_back(newTag); // Here we can insert more tag in the same offset because of bitfield meaning
	  }
}

TagDirectoryTable::TagDirectoryTable (TagDirectory* p, FILE* f, int memsize,int offs, TagType type, const TagAttrib* ta, ByteOrder border)
:TagDirectory(p,ta,border),zeroOffset(offs),valuesSize(memsize),defaultType( type )
{
	  values = new unsigned char[valuesSize];
	  fread (values, 1, valuesSize, f);

	  for( const TagAttrib* tattr = ta; tattr->ignore != -1; tattr++){
		  Tag* newTag = new Tag (this, tattr, (values + zeroOffset+ tattr->ID*getTypeSize(type)),type);
		  tags.push_back(newTag); // Here we can insert more tag in the same offset because of bitfield meaning
	  }
}
TagDirectory* TagDirectoryTable::clone (TagDirectory* parent) {

    TagDirectory* td = new TagDirectoryTable (parent,values,valuesSize,zeroOffset,defaultType, attribs, order);
    return td;
}

TagDirectoryTable::~TagDirectoryTable()
{
	if(values)
		delete [] values;
}
int TagDirectoryTable::calculateSize ()
{
	return valuesSize;
}

int TagDirectoryTable::write (int start, unsigned char* buffer) {
	if( values && valuesSize){
        memcpy(buffer+start,values,valuesSize);
        return start+valuesSize;
	}else
		return start;
}

//--------------- class Tag ---------------------------------------------------
// this class represents a tag stored in the directory
//-----------------------------------------------------------------------------

Tag::Tag (TagDirectory* p, FILE* f, int base) 
  : type(INVALID), count(0), value(NULL), allocOwnMemory(true), attrib(NULL), parent(p), directory(NULL) {

  tag   = get2 (f, getOrder());
  type  = (TagType)get2 (f, getOrder());
  count = get4 (f, getOrder());

  makerNoteKind = NOMK;
  keep = false;

  // filter out invalid tags
  if ((int)type<1 || (int)type>14 || count>900000) {
    type = INVALID;
    return;
  }

  // save file position
  int save = ftell(f) + 4;

  // load value field (possibly seek before)
  valuesize = count * getTypeSize(type);

  if (valuesize > 4) 
    fseek (f, get4(f, getOrder()) + base, SEEK_SET);

  attrib = parent->getAttrib (tag); 

  if (attrib && (attrib->action==1 || attrib->action==3))
    keep = true;

  // if this tag is the makernote, it needs special treatment (brand specific parsing)
  if (tag==0x927C && attrib && !strcmp (attrib->name, "MakerNote")) {
    value = NULL;
    // select format of makernote
    char make[128], model[128];
    Tag* tmake = parent->getParent()->getTag ("Make");
    if (tmake) 
        tmake->toString (make);
    else
        make[0] = 0;
    Tag* tmodel = parent->getParent()->getTag ("Model");
    if (tmodel)
        tmodel->toString (model);
    else
        model[0] = 0;
    if (!strncmp(make, "NIKON", 5)) {
        if (!strncmp(model, "NIKON E700",10)||!strncmp(model, "NIKON E800",10)||!strncmp(model, "NIKON E900",10)||!strncmp(model, "NIKON E900S",11)||!strncmp(model, "NIKON E910", 10)||!strncmp(model, "NIKON E950", 10)) {
            makerNoteKind = HEADERIFD;
            valuesize = 8;
            value = new unsigned char[8];
            fread (value, 1, 8, f);
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, base, nikon2Attribs, getOrder());
            directory[1] = NULL;
        }
        else if (!strncmp(model, "NIKON E990",10)||(!strncmp(model, "NIKON D1",8) && model[8]!='0')) {
            makerNoteKind = IFD;
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, base, nikon3Attribs, getOrder());
            directory[1] = NULL;
        }
        else {
            // needs refinement! (embedded tiff header parsing)
            makerNoteKind = NIKON3;
            valuesize = 18;
            value = new unsigned char[18];
            int basepos = ftell (f);
            fread (value, 1, 18, f);
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, basepos+10, nikon3Attribs, getOrder());
            directory[1] = NULL;
        }
    }
    else if (!strncmp(make, "Canon", 5)) {
        makerNoteKind = IFD;
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, canonAttribs, getOrder());
        directory[1] = NULL;
    }
    else if (!strncmp(make, "PENTAX", 6)) {
        makerNoteKind = HEADERIFD;
        valuesize = 6;
        value = new unsigned char[6];
        fread (value, 1, 6, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, pentaxAttribs, getOrder());
        directory[1] = NULL;
    }    
    else if (!strncmp(make, "FUJIFILM", 8)) {
        makerNoteKind = FUJI;
        valuesize = 12;
        value = new unsigned char[12];
        fread (value, 1, 12, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, ftell(f)-12, fujiAttribs, INTEL);
        directory[1] = NULL;
    }    
    else if (!strncmp(make, "KONICA MINOLTA", 14) || !strncmp(make, "Minolta", 7)) {
        makerNoteKind = IFD;
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, minoltaAttribs, getOrder());
        directory[1] = NULL;
    }    
    else if (!strncmp(make, "SONY", 4)) {
        valuesize = 12;
        value = new unsigned char[12];
        fread (value, 1, 12, f);
        if (!strncmp((char*)value, "SONY DSC", 8)) 
            makerNoteKind = HEADERIFD;
        else {
            makerNoteKind = IFD;
            fseek (f, -12, SEEK_CUR);
        }
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, sonyAttribs, getOrder());
        directory[1] = NULL;
    }    
    else if (!strncmp(make, "OLYMPUS", 7)) {
        makerNoteKind = HEADERIFD;
        valuesize = 8;
        value = new unsigned char[12];
        fread (value, 1, 8, f);
        directory = new TagDirectory*[2];
        directory[1] = NULL;
        if (!strncmp((char*)value, "OLYMPUS", 7)) {
            makerNoteKind = OLYMPUS2;
            fread (value+8, 1, 4, f);
            valuesize = 12;
            directory[0] = new TagDirectory (parent, f, ftell(f)-12, olympusAttribs, value[8]=='I' ? INTEL : MOTOROLA);
        }
        else 
            directory[0] = new TagDirectory (parent, f, base, olympusAttribs, getOrder());
    }    
    else {
        type = INVALID;
        fseek (f, save, SEEK_SET);
        return;
    }
  }
  else if (attrib && attrib->subdirAttribs) {
	    // Some subdirs are specific of maker and model
	    char make[128], model[128];
	    Tag* tmake = parent->getRoot()->getTag ("Make");
	    if (tmake) tmake->toString (make);
	    else       make[0] = 0;
	    Tag* tmodel = parent->getRoot()->getTag ("Model");
	    if (tmodel) tmodel->toString (model);
	    else        model[0] = 0;


        if (!strncmp(make, "SONY", 4)) {
			switch( tag ){
			case 0x0114:
				{
			        directory = new TagDirectory*[2];
			        directory[1] = NULL;
					if( strstr(model, "A330")  || strstr(model, "A380") )
					   directory[0] = new TagDirectoryTable (parent, f, valuesize*2,0,SHORT , sonyCameraSettingsAttribs2, MOTOROLA);
					else
					   directory[0] = new TagDirectoryTable (parent, f, valuesize*2,0,SHORT , sonyCameraSettingsAttribs, MOTOROLA);
					makerNoteKind = TABLESUBDIR;
				}
				break;
			default:
				goto defsubdirs;
			}
        }else if (!strncmp(make, "PENTAX", 6)) {
			switch( tag ){
			case 0x005c:
			case 0x0205:
			case 0x0206:
			case 0x0208:
			case 0x0216:
		        directory = new TagDirectory*[2];
		        directory[1] = NULL;
		        directory[0] = new TagDirectoryTable (parent, f, valuesize,0,BYTE , attrib->subdirAttribs, getOrder());
		        makerNoteKind = TABLESUBDIR;
				break;
			case 0x0215:
		        directory = new TagDirectory*[2];
		        directory[1] = NULL;
		        directory[0] = new TagDirectoryTable (parent, f, valuesize,0,LONG , attrib->subdirAttribs, getOrder());
		        makerNoteKind = TABLESUBDIR;
				break;
			case 0x0207:
				{   // There are 2 format pentaxLensDataAttribs
					int offsetFirst = 4;
					if( strstr(model, "*ist") || strstr(model, "GX-1") || strstr(model, "K100D") || strstr(model, "K110D") )
						offsetFirst = 3;
					if( strstr(model, "K-5") || strstr(model, "K-r") )
						offsetFirst = 12;
			        directory = new TagDirectory*[2];
			        directory[1] = NULL;
					directory[0] = new TagDirectoryTable (parent, f, valuesize,offsetFirst,BYTE , attrib->subdirAttribs, getOrder());
					makerNoteKind = TABLESUBDIR;
				}
				break;
			default:
				goto defsubdirs;
			}
        }else if (!strncmp(make, "Canon", 5)) {
			switch( tag ){
			case 0x0001:
			case 0x0002:
			case 0x0004:
			case 0x0005:
			case 0x0093:
			case 0x0098:
			case 0x00a0:
		        directory = new TagDirectory*[2];
		        directory[1] = NULL;
		        directory[0] = new TagDirectoryTable (parent, f, valuesize,0,SSHORT , attrib->subdirAttribs, getOrder());
		        makerNoteKind = TABLESUBDIR;
				break;
			case 0x009a:
			case 0x4013:
		        directory = new TagDirectory*[2];
		        directory[1] = NULL;
		        directory[0] = new TagDirectoryTable (parent, f, valuesize,0,LONG , attrib->subdirAttribs, getOrder());
		        makerNoteKind = TABLESUBDIR;
				break;
			default:
				goto defsubdirs;
			}
        }else if(type==UNDEFINED){
            count = 1;
            type = LONG;
		directory = new TagDirectory*[2];
        	directory[0] = new TagDirectory (parent, f, base, attrib->subdirAttribs, getOrder());
        	directory[1] = NULL;
        }else
        	goto defsubdirs;
  }
  else {
    // read value
    value = new unsigned char [valuesize];
    fread (value, 1, valuesize, f);
  }
  // seek back to the saved position
  fseek (f, save, SEEK_SET);
  return;
defsubdirs:
		// read value
		value = new unsigned char [valuesize];
		fread (value, 1, valuesize, f);
          int pos = ftell (f);
          // count the number of valid subdirs
          int sdcount = count;
          if (sdcount>0) {
              if (parent->getAttribTable()==olympusAttribs)
                  sdcount = 1;
              // allocate space
              directory = new TagDirectory*[sdcount+1];
              // load directories
              for (int j=0,i=0; j<count; j++,i++) {
                  int newpos = base + toInt(j*4, LONG);
                  fseek (f, newpos, SEEK_SET);
                  directory[i] = new TagDirectory (parent, f, base, attrib->subdirAttribs, getOrder());
                  fseek (f, pos, SEEK_SET);
              }
              // set the terminating NULL
              directory[sdcount] = NULL;
           }
           else
              type = INVALID;

// seek back to the saved position
fseek (f, save, SEEK_SET);
return;

}

Tag* Tag::clone (TagDirectory* parent) {

    Tag* t = new Tag (parent, attrib);
    
    t->tag = tag;
    t->type = type;
    t->count = count;
    t->keep = keep;
    t->valuesize = valuesize;
    if (value) {
        t->value = new unsigned char [valuesize];
        memcpy (t->value, value, valuesize);
    }
    else
        value = NULL;
    t->makerNoteKind = makerNoteKind;
    if (directory) {
        int ds = 0;
        for (; directory[ds]; ds++);
        t->directory = new TagDirectory*[ds+1];
        for (int i=0; i<ds; i++)
            t->directory[i] = directory[i]->clone (parent);
        t->directory[ds] = NULL;
    }
    else
        t->directory = NULL;
    return t;
}

Tag::~Tag () {

  // delete value
  if (value && allocOwnMemory)
    delete [] value;
  
  // if there are directories behind the tag, delete them
  int i = 0;
  if (directory) {
    while (directory[i])
        delete directory[i++];
    delete [] directory;
  }  
}

void Tag::setInt (int v, int ofs, TagType astype) {

    if (astype==SHORT)
        sset2 (v, value+ofs, getOrder());
    else if (astype==RATIONAL) {
        sset4 (v, value+ofs, getOrder());
        sset4 (1, value+ofs+4, getOrder());
    }
    else
        sset4 (v, value+ofs, getOrder());
}

void Tag::fromInt (int v) {

    if (type==SHORT)
        sset2 (v, value, getOrder());
    else
        sset4 (v, value, getOrder());
}

void Tag::fromString (const char* v, int size) {

	if( value && allocOwnMemory)
        delete [] value;
    if (size<0)
        valuesize = strlen (v) + 1;
    else
        valuesize = size;
    count = valuesize;
    if( allocOwnMemory )
       value = new unsigned char [valuesize];
    memcpy ((char*)value, v, valuesize);
}

int Tag::toInt (int ofs, TagType astype) {

  int a;
  if (astype == INVALID)
    astype = type;
  switch (astype) {
    case BYTE:  return value[ofs];
    case ASCII: return 0;
    case SSHORT:return (int)int2_to_signed(sget2 (value+ofs, getOrder())); 
    case SHORT: return (int)sget2 (value+ofs, getOrder());
    case SLONG:  
    case LONG:  return (int)sget4 (value+ofs, getOrder());
    case SRATIONAL: 
    case RATIONAL: a = (int)sget4 (value+ofs+4, getOrder()); return a==0 ? 0 : (int)sget4 (value+ofs, getOrder()) / a;
    case FLOAT:     return (int)toDouble(ofs);
    case UNDEFINED: return 0;
    default: return 0; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
  }
  return 0;
}

double Tag::toDouble (int ofs) {
  union IntFloat { uint32_t i; float f; } conv;

  double ud, dd;
  switch (type) {
    case BYTE:  return (double)((int)value[ofs]);
    case ASCII: return 0.0;
    case SSHORT:return (double)int2_to_signed(sget2 (value+ofs, getOrder()));  
    case SHORT: return (double)((int)sget2 (value+ofs, getOrder()));
    case SLONG:  
    case LONG:  return (double)((int)sget4 (value+ofs, getOrder()));
    case SRATIONAL: 
    case RATIONAL: ud = (int)sget4 (value+ofs, getOrder()); dd = (int)sget4 (value+ofs+4, getOrder()); return dd==0. ? 0. : (double)ud / (double)dd;
    case FLOAT:     
        conv.i=sget4 (value+ofs, getOrder());
        return conv.f;  // IEEE FLOATs are already C format, they just need a recast

    case UNDEFINED: return 0.;
    default: return 0.; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
  }
  return 0.;
}

void Tag::toRational (int& num, int& denom, int ofs) {

  switch (type) {
    case BYTE:   num = (int)value[ofs]; denom = 1; break;
    case ASCII:  num = 0; denom = 0; break;
    case SSHORT: 
    case SHORT:  num = (int)sget2 (value+ofs, getOrder()); denom = 1; break;
    case SLONG:  
    case LONG:   num = (int)sget4 (value+ofs, getOrder()); denom = 1; break;
    case SRATIONAL: 
    case RATIONAL: num = (int)sget4 (value+ofs, getOrder()); denom = (int)sget4 (value+ofs+4, getOrder()); break;
    case FLOAT:   num = 0; denom = 0; break;
    case UNDEFINED: num = 0; denom = 0; break;
    default: num = 0; denom = 0; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
  }
}

void Tag::toString (char* buffer, int ofs) {

  if (type==UNDEFINED && !directory) {
      bool isstring = true;
      int i=0;
      for (i=0; i+ofs<count && i<64 && value[i+ofs]; i++)
        if (value[i+ofs]<32 || value[i+ofs]>126)
            isstring = false;
      if (isstring) {
        int j = 0;
        for (i=0; i+ofs<count && i<64 && value[i+ofs]; i++) {
            if (value[i+ofs]=='<' || value[i+ofs]=='>')
                buffer[j++] = '\\';
            buffer[j++] = value[i+ofs];
        }
        buffer[j++] = 0;     
        return;
      }
  }
  else if (type==ASCII) {
    sprintf (buffer, "%.64s", value+ofs);
    return; 
  }

  int maxcount = 4;
  if (count<4)
    maxcount = count;
  
  strcpy (buffer, "");
  for (int i=0; i<maxcount; i++) {
    if (i>0)
        strcat (buffer, ", ");
    char* b = buffer + strlen(buffer);
        
    switch (type) {
        case UNDEFINED: 
        case BYTE:  sprintf (b, "%d", value[i+ofs]); break;
        case SSHORT: sprintf (b, "%d", toInt(2*i+ofs)); break;
        case SHORT: sprintf (b, "%u", toInt(2*i+ofs)); break;
        case SLONG: sprintf (b, "%ld", (long)toInt(4*i+ofs)); break;
        case LONG:  sprintf (b, "%lu", (unsigned long)toInt(4*i+ofs)); break;
        case SRATIONAL: 
        case RATIONAL: sprintf (b, "%d/%d", (int)sget4 (value+8*i+ofs, getOrder()), (int)sget4 (value+8*i+ofs+4, getOrder())); break; 
        case FLOAT:    sprintf (b, "%g", toDouble(8*i+ofs)); break;
	default: break;
    }
  }
  if (count > maxcount)
    strcat (buffer, "...");
}

std::string Tag::nameToString (int i) {

    char buffer[1024];
    if (attrib) 
        strcpy (buffer, attrib->name);
    else 
        sprintf (buffer, "0x%x", tag);
    if (i>0)
        sprintf (buffer+strlen(buffer)-1, "[%d]", i);
    return buffer;        
}

std::string Tag::valueToString () {

    char buffer[1024];
    if (attrib && attrib->interpreter) 
        return attrib->interpreter->toString (this);
    else {
        toString (buffer);
        return buffer;
    }
}

void Tag::valueFromString (const std::string& value) {

    if (attrib && attrib->interpreter) 
        attrib->interpreter->fromString (this, value);
}

int Tag::calculateSize () {
    int size = 0;
    
    if (directory) {
        int j;
        for (j=0; directory[j]; j++)
            size += directory[j]->calculateSize ();
        if (j>1)
            size += 4*j;
    }
    else if (valuesize > 4) 
      size += valuesize + (valuesize%2); // we align tags to even byte positions 

   if (makerNoteKind!=NOMK) 
        count = directory[0]->calculateSize ();
        
   if (makerNoteKind==NIKON3 || makerNoteKind==OLYMPUS2 || makerNoteKind==FUJI) 
        size += valuesize;
   else if (makerNoteKind==HEADERIFD)
        size += valuesize;
    return size;
}    

int Tag::write (int offs, int dataOffs, unsigned char* buffer) {

  if ((int)type==0 || offs>65500)
    return dataOffs;

  sset2 (tag, buffer+offs, parent->getOrder());
  offs += 2;
  unsigned short typ = type;
  sset2 (typ, buffer+offs, parent->getOrder());
  offs += 2;
  sset4 (count, buffer+offs, parent->getOrder());
  offs += 4;
  if (!directory) {
    if (valuesize>4) {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        memcpy (buffer+dataOffs, value, valuesize);
        if (valuesize%2)
            buffer[dataOffs+valuesize] = 0; // zero padding required by the exif standard
        return dataOffs + valuesize + (valuesize%2);
    }
    else {
        memcpy (buffer+offs, value, valuesize);
        return dataOffs;
    }
  }
  else {
    if (makerNoteKind==NIKON3) {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        memcpy (buffer+dataOffs, value, 18);
        dataOffs += 10;
        dataOffs += directory[0]->write (8, buffer+dataOffs);
        return dataOffs;
    }
    else if (makerNoteKind==OLYMPUS2 || makerNoteKind==FUJI) {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        memcpy (buffer+dataOffs, value, valuesize);
        dataOffs += valuesize + directory[0]->write (valuesize, buffer+dataOffs);
        return dataOffs;
    }
    else if (makerNoteKind==HEADERIFD) {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        memcpy (buffer+dataOffs, value, valuesize);
        dataOffs += valuesize;
        dataOffs += directory[0]->write (dataOffs, buffer);
        return dataOffs;
    }
    else if( makerNoteKind==TABLESUBDIR){
    	sset4 (dataOffs, buffer+offs, parent->getOrder());
        dataOffs = directory[0]->write (dataOffs, buffer);
        return dataOffs;
    }
    else if (!directory[1]) {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        return directory[0]->write (dataOffs, buffer);
    }
    else {
        sset4 (dataOffs, buffer+offs, parent->getOrder());
        int linkOffs = dataOffs;
        for (int i=0; directory[i]; i++)
            dataOffs += 4;              
        for (int i=0; directory[i]; i++) {
            sset4 (dataOffs, buffer+linkOffs, parent->getOrder());
            linkOffs += 4;
            dataOffs = directory[i]->write (dataOffs, buffer);
        }
        return dataOffs;
    }
  }
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr)
 : tag(attr ? attr->ID : -1), type(INVALID), count(0), value(NULL), valuesize(0), keep(true), allocOwnMemory(true), attrib(attr), parent(p), directory(NULL), makerNoteKind (NOMK) {
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, int data, TagType t) 
 : tag(attr ? attr->ID : -1), type(t), count(1), value(NULL), valuesize(0), keep(true), allocOwnMemory(true), attrib(attr), parent(p), directory(NULL), makerNoteKind (NOMK) {

    initInt (data, t);
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, unsigned char *data, TagType t)
 : tag(attr ? attr->ID : -1), type(t), count(1), value(NULL), valuesize(0), keep(true), allocOwnMemory(false), attrib(attr), parent(p), directory(NULL), makerNoteKind (NOMK) {

    initType (data, t);
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, const char* text) 
 : tag(attr ? attr->ID : -1), type(ASCII), count(1), value(NULL), valuesize(0), keep(true), allocOwnMemory(true), attrib(attr), parent(p), directory(NULL), makerNoteKind (NOMK) {

    initString (text);
}

void Tag::initType (unsigned char *data, TagType type)
{
    valuesize = getTypeSize(type);
    if( allocOwnMemory ){
       value = new unsigned char[valuesize];
       memcpy ((char*)value, data, valuesize);
    }else
    	value = data;
}

void Tag::initInt (int data, TagType t, int cnt) {
  
    type = t;
    if (t==LONG) 
        valuesize = 4;
    else if (t==SHORT) 
        valuesize = 2;
    else if (t==RATIONAL) 
        valuesize = 8;

    count = cnt;
    valuesize *= count;
    value = new unsigned char[valuesize];
    setInt (data, 0, t);
}
  
void Tag::initString (const char* text) {

    type = ASCII;
    count = strlen(text)+1;
    valuesize = count;
    value = new unsigned char[valuesize];
    strcpy ((char*)value, text);
}

void Tag::initSubDir () {
    type = LONG;
    valuesize = 4;
    count = 1;
    value = new unsigned char[4];
    setInt (0);
    directory = new TagDirectory*[2];
    directory[0] = new TagDirectory (parent, attrib ? attrib->subdirAttribs : NULL, parent->getOrder());
    directory[1] = NULL;
}

void Tag::initMakerNote (MNKind mnk, const TagAttrib* ta) {
    type = UNDEFINED;
    valuesize = 4;
    count = 1;
    value = new unsigned char[4];
    setInt (0);
    directory = new TagDirectory*[2];
    directory[0] = new TagDirectory (parent, ta, parent->getOrder());
    directory[1] = NULL;
    makerNoteKind = mnk;
}

void Tag::initUndefArray (const char* data, int len) {
    type = UNDEFINED;
    count = valuesize = len;
    value = new unsigned char[valuesize];
    memcpy (value, data, len);
}

void Tag::initLongArray (const char* data, int len) {
    type = LONG;
    count = (len+3)/4;
    valuesize = count * 4;
    value = new unsigned char[valuesize];
    memcpy (value, data, len);
}

void Tag::initRational (int num, int den) {
    count = 1;
    valuesize = 8;
    value = new unsigned char[8];
    type = RATIONAL;
    setInt (num, 0);
    setInt (den, 4);
}

//--------------- class IFDParser ---------------------------------------------
// static functions to read tag directories from different kinds of files
//-----------------------------------------------------------------------------


const TagAttrib* lookupAttrib (const TagAttrib* dir, const char* field) {

    for (int i=0; dir[i].ignore!=-1; i++)
        if (!strcmp (dir[i].name, field))
            return &dir[i];
    return 0;
}


TagDirectory* ExifManager::parseCIFF (FILE* f, int base, int length) {

    TagDirectory* root = new TagDirectory (NULL, ifdAttribs, INTEL);
    Tag* exif = new Tag (root, lookupAttrib(ifdAttribs,"Exif"));
    exif->initSubDir ();
    Tag* mn = new Tag (exif->getDirectory(), lookupAttrib(exifAttribs,"MakerNote"));
    mn->initMakerNote (IFD, canonAttribs);
    root->addTag (exif);
    exif->getDirectory()->addTag (mn);
    parseCIFF (f, base, length, root);
    root->sort ();
    return root;
}

Tag* ExifManager::saveCIFFMNTag (FILE* f, TagDirectory* root, int len, const char* name) {
    int s = ftell (f);
    char* data = new char [len];
    fread (data, len, 1, f);
    TagDirectory* mn = root->getTag ("Exif")->getDirectory()->getTag("MakerNote")->getDirectory();
    Tag* cs = new Tag (mn, lookupAttrib(canonAttribs, name));
    cs->initUndefArray (data, len);
    mn->addTag (cs);
    fseek (f, s, SEEK_SET);
    return cs;
}

void ExifManager::parseCIFF (FILE* f, int base, int length, TagDirectory* root) {

  char buffer[1024];
  Tag* t;

  fseek (f, base+length-4, SEEK_SET);

  int dirStart = get4 (f, INTEL) + base;
  fseek (f, dirStart, SEEK_SET);

  int numOfTags = get2 (f, INTEL);

  if (numOfTags > 100) return;
  
  float exptime, shutter, aperture, fnumber, ev;
  exptime = fnumber = shutter = aperture = ev = -1000;
  int focal_len, iso;
  focal_len = iso = -1;
  
  TagDirectory* exif = root->getTag("Exif")->getDirectory();
  
  time_t timestamp = time (NULL);
  
  for (int i=0; i<numOfTags; i++) {

    int type = get2 (f, INTEL);
    int len  = get4 (f, INTEL);
    int nextPos = ftell (f) + 4;
    
    // seek to the location of the value
    fseek (f, base + get4 (f, INTEL), SEEK_SET);
  
    if ((((type >> 8) + 8) | 8) == 0x38)
      parseCIFF (f, ftell(f), len, root);	// Parse a sub-table 
     
    if (type == 0x0810) {
      fread (buffer, 64, 1, f);
      t = new Tag (root, lookupAttrib(ifdAttribs,"Artist"));
      t->initString (buffer);
      root->addTag (t);
    }
    if (type == 0x080a) {
      fread (buffer, 64, 1, f);
      t = new Tag (root, lookupAttrib(ifdAttribs,"Make"));
      t->initString (buffer);
      root->addTag (t);
      fseek (f, strlen(buffer) - 63, SEEK_CUR);
      fread (buffer, 64, 1, f);
      t = new Tag (root, lookupAttrib(ifdAttribs,"Model"));
      t->initString (buffer);
      root->addTag (t);
    }
    if (type == 0x1818) {
      ev = int_to_float(get4(f, INTEL));
      shutter = int_to_float(get4(f, INTEL));
      exptime = pow (2, -shutter);
      aperture = int_to_float(get4(f, INTEL));
      fnumber = pow (2, aperture/2);
      
    }
    if (type == 0x102d) {
        Tag* t = saveCIFFMNTag (f, root, len, "CanonCameraSettings");
        int mm = t->toInt (34, SHORT);
        Tag* nt = new Tag (exif, lookupAttrib(exifAttribs,"MeteringMode"));
        switch (mm) {
            case 0: nt->initInt (5, SHORT); break;
            case 1: nt->initInt (3, SHORT); break;
            case 2: nt->initInt (1, SHORT); break;
            case 3: nt->initInt (5, SHORT); break;
            case 4: nt->initInt (6, SHORT); break;
            case 5: nt->initInt (2, SHORT); break;
        }
        exif->addTag (nt);
        nt = new Tag (exif, lookupAttrib(exifAttribs,"MaxApertureValue"));
        nt->initRational (t->toInt(52,SHORT), 32);
        exif->addTag (nt);
        int em = t->toInt(40,SHORT);
        nt = new Tag (exif, lookupAttrib(exifAttribs,"ExposureProgram"));
        switch (em) {
            case 0: nt->initInt (2, SHORT); break;
            case 1: nt->initInt (2, SHORT); break;
            case 2: nt->initInt (4, SHORT); break;
            case 3: nt->initInt (3, SHORT); break;
            case 4: nt->initInt (1, SHORT); break;
            default: nt->initInt (0, SHORT); break;
        }
        exif->addTag (nt);
        nt = new Tag (exif, lookupAttrib(exifAttribs,"Flash"));
        if (t->toInt(8,SHORT)==0)
            nt->initInt (0, SHORT);
        else
            nt->initInt (1, SHORT);
        exif->addTag (nt);        
        nt = new Tag (exif, lookupAttrib(exifAttribs,"MaxApertureValue"));
        nt->initRational (t->toInt(52,SHORT), 32);
        exif->addTag (nt);        
    }
    if (type == 0x1029) 
        saveCIFFMNTag (f, root, len, "CanonFocalLength");
    if (type == 0x1031) 
        saveCIFFMNTag (f, root, len, "SensorInfo");
    if (type == 0x1033) 
        saveCIFFMNTag (f, root, len, "CustomFunctions");
    if (type == 0x1038) 
        saveCIFFMNTag (f, root, len, "CanonAFInfo");
    if (type == 0x1093) 
        saveCIFFMNTag (f, root, len, "CanonFileInfo");
    if (type == 0x10a9) 
        saveCIFFMNTag (f, root, len, "ColorBalance");
    if (type == 0x102a) {
        saveCIFFMNTag (f, root, len, "CanonShotInfo");

        iso = pow (2, (get4(f, INTEL),get2(f, INTEL))/32.0 - 4) * 50;
        aperture  = ((get2(f, INTEL),(short)get2(f, INTEL))/32.0);
        fnumber = pow (2, aperture/2);
        shutter = ((short)get2(f, INTEL))/32.0;
        ev = ((short)get2(f, INTEL))/32.0;
        fseek (f, 34, SEEK_CUR);
        if (shutter > 1e6) shutter = get2 (f, INTEL) / 10.0; 
        exptime   = pow (2,-shutter);
    }
    if (type == 0x5029) {
      focal_len = len >> 16;
      if ((len & 0xffff) == 2) focal_len /= 32;
    }
//    if (type == 0x5813) flash_used = int_to_float(len);
    if (type == 0x580e) timestamp  = len;
    if (type == 0x180e) timestamp  = get4 (f, INTEL); 
    if ((type | 0x4000) == 0x580e) 
      timestamp = mktime (gmtime (&timestamp));
    fseek (f, nextPos, SEEK_SET);
  }
  if (shutter>-999) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"ShutterSpeedValue"));
    t->initRational ((int)(shutter*10000), 10000);
    exif->addTag (t);
  }
  if (exptime>-999) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"ExposureTime"));
    t->initRational ((int)(exptime*10000), 10000);
    exif->addTag (t);
  }
  if (aperture>-999) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"ApertureValue"));
    t->initRational ((int)(aperture*10), 10);
    exif->addTag (t);
  }
  if (fnumber>-999) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"FNumber"));
    t->initRational ((int)(fnumber*10), 10);
    exif->addTag (t);
  }
  if (ev>-999) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"ExposureBiasValue"));
    t->initRational ((int)(ev*1000), 1000);
    exif->addTag (t);
  }
  if (iso>0) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"ISOSpeedRatings"));
    t->initInt (iso, LONG);
    exif->addTag (t);
  }
  if (focal_len>0) {
    t = new Tag (exif, lookupAttrib(exifAttribs,"FocalLength"));
    t->initRational (focal_len*32, 32);
    exif->addTag (t);
  }

  if (timestamp!=time(NULL)) {
    struct tm* tim = localtime (&timestamp);
    strftime (buffer, 20, "%Y:%m:%d %H:%M:%S", tim);
    t = new Tag (exif, lookupAttrib(exifAttribs,"DateTimeOriginal"));
    t->initString (buffer);
    exif->addTag (t);
    t = new Tag (exif, lookupAttrib(exifAttribs,"DateTimeDigitized"));
    t->initString (buffer);
    exif->addTag (t);
    t = new Tag (root, lookupAttrib(ifdAttribs,"DateTime"));
    t->initString (buffer);
    root->addTag (t);
  }
}

TagDirectory* ExifManager::parse (FILE* f, int base, bool skipIgnored) {
  setlocale(LC_NUMERIC, "C"); // to set decimal point in sscanf
  // read tiff header
  fseek (f, base, SEEK_SET);
  unsigned short bo; 
  fread (&bo, 1, 2, f);
  ByteOrder order = (ByteOrder)((int)bo);
  get2 (f, order);
  int firstifd = get4 (f, order);

  // seek to IFD0
  fseek (f, base+firstifd, SEEK_SET);

  // first read the IFD directory
  TagDirectory* root =  new TagDirectory (NULL, f, base, ifdAttribs, order, skipIgnored);

  // fix ISO issue with nikon and panasonic cameras
  Tag* exif = root->getTag ("Exif");
  if (exif && !exif->getDirectory()->getTag("ISOSpeedRatings")) {
    Tag* make = root->getTag ("Make");
    if (make && !strncmp((char*)make->getValue(), "NIKON", 5)) {
        Tag* mn   = exif->getDirectory()->getTag("MakerNote");
        if (mn) {
            Tag* iso = mn->getDirectory()->getTag("ISOSpeed");
            if (iso) {
                std::string isov = iso->valueToString ();
                Tag* niso = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
                niso->initInt (atoi(isov.c_str()), SHORT);
                exif->getDirectory()->addTagFront (niso);
            }
        }
    }
    else if (make && (!strncmp((char*)make->getValue(), "Panasonic", 9) || !strncmp((char*)make->getValue(), "LEICA", 5))) {
        Tag* iso = root->getTag("PanaISO");
        if (iso) {
            std::string isov = iso->valueToString ();
            Tag* niso = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
            niso->initInt (atoi(isov.c_str()), SHORT);
            exif->getDirectory()->addTagFront (niso);
        }
    }
  }

// root->printAll ();

  return root;
}

TagDirectory* ExifManager::parseJPEG (FILE* f) {

  fseek (f, 0, SEEK_SET);
  unsigned char markerl = 0xff;
  unsigned char c;
  fread (&c, 1, 1, f);
  const char exifid[] = "Exif\0\0";
  char idbuff[8];
  int tiffbase = -1;
  while (fread (&c, 1, 1, f)) {
    if (c!=markerl) continue;
    if (fread (&c, 1, 1, f) && c==0xe1) {  // APP1 marker found
      if (fread (idbuff, 1, 8, f)<8)
        return NULL;
      if (!memcmp(idbuff+2, exifid, 6)) {     // Exif info found
        tiffbase = ftell (f);        
        return parse (f, tiffbase);
      }
    }
  }
  return NULL;
}

TagDirectory* ExifManager::parseTIFF (FILE* f, bool skipIgnored) {

  return parse (f, 0, skipIgnored);
}

std::vector<Tag*> ExifManager::defTags; 

// forthis: the byte order will be taken from directory "forthis"
const std::vector<Tag*>& ExifManager::getDefaultTIFFTags (TagDirectory* forthis) {

    for (int i=0; i<defTags.size(); i++)
        delete defTags[i];
    defTags.clear ();
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"ImageWidth"), 0, LONG));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"ImageHeight"), 0, LONG));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"XResolution"), 300, RATIONAL));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"YResolution"), 300, RATIONAL));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"ResolutionUnit"), 2, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"Software"), "RawTherapee"));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"Orientation"), 1, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"SamplesPerPixel"), 3, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"BitsPerSample"), 8, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"PlanarConfiguration"), 1, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"PhotometricInterpretation"), 2, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib(ifdAttribs,"Compression"), 1, SHORT));
    
    return defTags;
}



int ExifManager::createJPEGMarker (const TagDirectory* root, const rtengine::procparams::ExifPairs& changeList, int W, int H, unsigned char* buffer) {

  // write tiff header
  int offs = 6;
  memcpy (buffer, "Exif\0\0", 6);
  ByteOrder order = INTEL;
  if (root)
    order = root->getOrder ();
  sset2 ((unsigned short)order, buffer+offs, order); offs += 2;
  sset2 (42, buffer+offs, order); offs += 2;
  sset4 (8, buffer+offs, order); offs += 4;

  TagDirectory* cl;
  if (root)
    //FIXME: static_cast needed here
    cl = ((TagDirectory*)root)->clone (NULL);
  else
    cl = new TagDirectory (NULL, ifdAttribs, INTEL);

  for (rtengine::procparams::ExifPairs::const_iterator i=changeList.begin(); i!=changeList.end(); i++)
     cl->applyChange (i->first, i->second);
  
   getDefaultTIFFTags (cl);
  
  defTags[0]->setInt (W, 0, LONG);
  defTags[1]->setInt (H, 0, LONG);
  defTags[8]->setInt (8, 0, SHORT);
  
  for (int i=defTags.size()-1; i>=0; i--) 
    cl->replaceTag (defTags[i]->clone (cl));
  cl->sort ();
  int size = cl->write (8, buffer+6);

  delete cl;
  
  return size + 6;
}

int ExifManager::createTIFFHeader (const TagDirectory* root, const rtengine::procparams::ExifPairs& changeList, int W, int H, int bps, const char* profiledata, int profilelen, const char* iptcdata, int iptclen, unsigned char* buffer) {

// write tiff header
    int offs = 0;
    ByteOrder order = INTEL;
    if (root)
        order = root->getOrder ();
    sset2 ((unsigned short)order, buffer+offs, order); offs += 2;
    sset2 (42, buffer+offs, order); offs += 2;
    sset4 (8, buffer+offs, order); offs += 4;

    TagDirectory* cl;
    if (root)
	//FIXME: static_cast needed here
        cl = ((TagDirectory*)root)->clone (NULL);
    else
        cl = new TagDirectory (NULL, ifdAttribs, INTEL);

// add tiff strip data
    int rps = 8;
    int strips = ceil((double)H/rps);
    cl->replaceTag (new Tag (cl, lookupAttrib(ifdAttribs,"RowsPerStrip"), rps, LONG));
    Tag* stripBC   = new Tag (cl, lookupAttrib(ifdAttribs,"StripByteCounts"));
    stripBC->initInt (0, LONG, strips);
    cl->replaceTag (stripBC);
    Tag* stripOffs = new Tag (cl, lookupAttrib(ifdAttribs,"StripOffsets"));
    stripOffs->initInt (0, LONG, strips);
    cl->replaceTag (stripOffs);
    for (int i=0; i<strips-1; i++)
        stripBC->setInt (rps*W*3*bps/8, i*4);
    int remaining = (H-rps*floor((double)H/rps))*W*3*bps/8;
    if (remaining)
        stripBC->setInt (remaining, (strips-1)*4);   
    else
        stripBC->setInt (rps*W*3*bps/8, (strips-1)*4);
    if (profiledata) {
        Tag* icc = new Tag (cl, lookupAttrib(ifdAttribs,"ICCProfile"));
        icc->initUndefArray (profiledata, profilelen);
        cl->replaceTag (icc);
    }
    if (iptcdata) {
        Tag* iptc = new Tag (cl, lookupAttrib(ifdAttribs,"IPTCData"));
        iptc->initLongArray (iptcdata, iptclen);
        cl->replaceTag (iptc);
    }
    
// apply list of changes
    for (rtengine::procparams::ExifPairs::const_iterator i=changeList.begin(); i!=changeList.end(); i++)
        cl->applyChange (i->first, i->second);
  
  // append default properties   
        getDefaultTIFFTags (cl);

    defTags[0]->setInt (W, 0, LONG);
    defTags[1]->setInt (H, 0, LONG);
    defTags[8]->initInt(0, SHORT, 3);
    for (int i=0;i<3;i++) defTags[8]->setInt(bps, i*2, SHORT);
 
    for (int i=defTags.size()-1; i>=0; i--) 
        cl->replaceTag (defTags[i]->clone (cl));
   
// calculate strip offsets  
    int size = cl->calculateSize ();
    int byps = bps / 8;
    for (int i=0; i<strips; i++)
        stripOffs->setInt (size + 8 + i*rps*W*3*byps, i*4);

    cl->sort ();
    int endOffs = cl->write (8, buffer);
  
//  cl->printAll();
    delete cl;

    return endOffs;
}

//-----------------------------------------------------------------------------
// global functions to read byteorder dependent data
//-----------------------------------------------------------------------------
unsigned short sget2 (unsigned char *s, rtexif::ByteOrder order) {

  if (order == rtexif::INTEL)   return s[0] | s[1] << 8;
  else			return s[0] << 8 | s[1];
}

int sget4 (unsigned char *s, rtexif::ByteOrder order) {

  if (order == rtexif::INTEL)   return s[0] | s[1] << 8 | s[2] << 16 | s[3] << 24;
  else                  return s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
}

unsigned short get2 (FILE* f, rtexif::ByteOrder order) {

  unsigned char str[2] = { 0xff,0xff };
  fread (str, 1, 2, f);
  return rtexif::sget2 (str, order);
}

int get4 (FILE* f, rtexif::ByteOrder order) {

  unsigned char str[4] = { 0xff,0xff,0xff,0xff };
  fread (str, 1, 4, f);
  return rtexif::sget4 (str, order);
}

void sset2 (unsigned short v, unsigned char *s, rtexif::ByteOrder order) {

  if (order == rtexif::INTEL) {
    s[0] = v & 0xff; v >>= 8;
    s[1] = v;
  }
  else {
    s[1] = v & 0xff; v >>= 8;
    s[0] = v;
  }
}

void sset4 (int v, unsigned char *s, rtexif::ByteOrder order) {

  if (order == rtexif::INTEL) {
    s[0] = v & 0xff; v >>= 8;
    s[1] = v & 0xff; v >>= 8;
    s[2] = v & 0xff; v >>= 8;
    s[3] = v;
  }
  else {
    s[3] = v & 0xff; v >>= 8;
    s[2] = v & 0xff; v >>= 8;
    s[1] = v & 0xff; v >>= 8;
    s[0] = v;
  }
}

float int_to_float (int i) {
  union { int i; float f; } u;
  u.i = i;
  return u.f;
}

short int int2_to_signed (short unsigned int i) {
  union { short unsigned int i; short int s; } u;
  u.i = i;
  return u.s;
}

int getTypeSize( TagType type ){
	return ("11124811248484"[type<14?type:0]-'0');
}

/* Function to parse and extract focal length and aperture information from description
 * @fullname must conform to the following formats
 * <focal>mm f/<aperture>
 * <focal>-<focal>mm f/<aperture>
 * <focal>-<focal>mm f/<aperture>-<aperture>
 * NB: no space between separator '-'; no space between focal length and 'mm'
 */
bool extractLensInfo(std::string &fullname,double &minFocal, double &maxFocal, double &maxApertureAtMinFocal, double &maxApertureAtMaxFocal)
{
	 minFocal=0.0;
	 maxFocal=0.0;
	 maxApertureAtMinFocal=0.0;
	 maxApertureAtMaxFocal=0.0;
	 char buffer[1024];
	 strcpy(buffer,fullname.c_str());
	 char *pF = strstr(buffer,"f/" );
	 if( pF ){
		 sscanf(pF+2,"%lf-%lf",&maxApertureAtMinFocal,&maxApertureAtMaxFocal);
		 if(maxApertureAtMinFocal >0. && maxApertureAtMaxFocal==0.)
			 maxApertureAtMaxFocal= maxApertureAtMinFocal;

		 char *pMM = pF-3;
		 while( pMM[0]!= 'm' && pMM[1]!= 'm' && pMM>buffer) pMM--;
		 if( pMM[0]== 'm' && pMM[1]== 'm' ){
		    char *sp=pMM;
		    while( *sp != ' ' && sp > buffer )sp--;
		    sscanf(sp+1,"%lf-%lf",&minFocal,&maxFocal);
		    return true;
		 }
	 }
     return false;
}

}

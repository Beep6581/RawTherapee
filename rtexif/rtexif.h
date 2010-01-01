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
#ifndef _MEXIF3_
#define _MEXIF3_

#include <stdio.h>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

namespace rtexif {

enum TagType {INVALID=0, BYTE=1, ASCII=2, SHORT=3, LONG=4, RATIONAL=5, UNDEFINED=7, SSHORT=8, SLONG=9, SRATIONAL=10, FLOAT=11, DOUBLE=12, OLYUNDEF=13, SUBDIR=99};
enum ActionCode {DONTWRITE=0, WRITE=1, SYSTEM=2};
enum ByteOrder {INTEL=0x4949, MOTOROLA=0x4D4D};
enum MNKind {NOMK, IFD, HEADERIFD, NIKON3, OLYMPUS2, FUJI};

struct TIFFHeader {

  unsigned short byteOrder;
  unsigned short fixed;
  unsigned int   ifdOffset;
};

class Tag;
class Interpreter;

// structure of informations describing an exif tag
struct TagAttrib {
    int                 ignore;   // =0: never ignore, =1: always ignore, =2: ignore if the subdir type is reduced image, =-1: end of table
    int                 action;    //=0: dont write it to the output, =1: write it to the output, =2: dont write, dont show, =3: write, dont show
    int                 editable;
    const  TagAttrib*   subdirAttribs;     // =0 ->not subdir
    unsigned short      ID;
    const char*         name;
    Interpreter*        interpreter;
};

// a directory of tags
class TagDirectory {

  protected:
    std::vector<Tag*> tags;     // tags in the directory
    const TagAttrib*  attribs;  // descriptor table to decode the tags
    ByteOrder         order;    // byte order
    TagDirectory*     parent;   // parent directory (NULL if root)
    
  public:
    TagDirectory ();
    TagDirectory (TagDirectory* p, FILE* f, int base, const TagAttrib* ta, ByteOrder border);
    TagDirectory (TagDirectory* p, const TagAttrib* ta, ByteOrder border);
   ~TagDirectory ();
   
    inline ByteOrder getOrder      () const { return order; }   
    TagDirectory*    getParent     () { return parent; }
    inline int       getCount      () const { return tags.size (); }
    const TagAttrib* getAttrib     (int id);
    const TagAttrib* getAttrib     (const char* name);
    const TagAttrib* getAttribTable() { return attribs; }
    Tag*             getTag        (const char* name);
    Tag*             getTag        (int ID);
    void             addTag        (Tag* a);
    void             addTagFront   (Tag* a);
    void             replaceTag    (Tag* a);
    inline Tag*      getTagByIndex (int ix) { return tags[ix]; }
    inline void      setOrder      (ByteOrder bo) { order = bo; }   

    int           calculateSize ();
    int           write         (int start, unsigned char* buffer);
    TagDirectory* clone         (TagDirectory* parent);
    void          applyChange   (std::string field, std::string value);

    void printAll () const;
    void sort     ();
};

// a class representing a single tag
class Tag {

  protected:
    unsigned short tag;
    TagType        type;
    unsigned int   count;
    unsigned char* value;
    int            valuesize;
    bool           keep;

    const TagAttrib* attrib;
    TagDirectory*    parent;
    TagDirectory**   directory;
    MNKind           makerNoteKind;
    
  public:
    Tag (TagDirectory* parent, FILE* f, int base);                          // parse next tag from the file
    Tag (TagDirectory* parent, const TagAttrib* attr); 
    Tag (TagDirectory* parent, const TagAttrib* attr, int data, TagType t);  // create a new tag from array (used
    Tag (TagDirectory* parent, const TagAttrib* attr, const char* data);  // create a new tag from array (used
   ~Tag ();
   void initInt        (int data, TagType t, int count=1);
   void initString     (const char* text);
   void initSubDir     ();
   void initMakerNote  (MNKind mnk, const TagAttrib* ta);
   void initUndefArray (const char* data, int len);
   void initLongArray  (const char* data, int len);
   void initRational   (int num, int den);

    // get basic tag properties
    int                  getID     () const { return tag; }
    int                  getCount  () const { return count; }
    TagType              getType   () const { return type; }
    unsigned char*       getValue  () const { return value; }
    const TagAttrib*     getAttrib () const { return attrib; }
    inline ByteOrder     getOrder  () const { return parent ? parent->getOrder() : INTEL; }
    inline TagDirectory* getParent () const { return parent; }
    int                  getValueSize () const { return valuesize; }

    // read/write value
    int    toInt        (int ofs=0, TagType astype=INVALID);
    void   fromInt      (int v);
    double toDouble     (int ofs=0);
    void   toRational   (int& num, int& denom, int ofs=0);
    void   toString     (char* buffer, int ofs=0);
    void   fromString   (const char* v, int size=-1);
    void   setInt       (int v, int ofs=0, TagType astype=LONG);

    // additional getter/setter for more confortable use
    std::string valueToString   ();
    std::string nameToString    (int i=0);
    void        valueFromString (const std::string& value);
    
    // functions for writing
    int  calculateSize ();
    int  write         (int offs, int dataOffs, unsigned char* buffer);
    Tag* clone         (TagDirectory* parent);

    // to control if the tag shall be written 
    bool getKeep ()       { return keep; }
    void setKeep (bool k) { keep = k; }   
   
    // get subdirectory (there can be several, the last is NULL)
    bool           isDirectory  ()        { return directory!=NULL; }
    TagDirectory*  getDirectory (int i=0) { return directory[i]; }

    MNKind getMakerNoteFormat () { return makerNoteKind; }
 };

class ExifManager {

    static std::vector<Tag*> defTags;
  
    static Tag* saveCIFFMNTag (FILE* f, TagDirectory* root, int len, const char* name);
  public:
    static TagDirectory* parse (FILE*f, int base);
    static TagDirectory* parseJPEG (FILE*f);
    static TagDirectory* parseTIFF (FILE*f);
    static TagDirectory* parseCIFF (FILE* f, int base, int length);
    static void          parseCIFF (FILE* f, int base, int length, TagDirectory* root);
    
    static const std::vector<Tag*>& getDefaultTIFFTags (TagDirectory* forthis);
    static int    createJPEGMarker (const TagDirectory* root, const std::vector< std::pair<std::string,std::string> >& changeList, int W, int H, unsigned char* buffer);
    static int    createTIFFHeader (const TagDirectory* root, const std::vector< std::pair<std::string,std::string> >& changeList, int W, int H, int bps, const char* profiledata, int profilelen, const char* iptcdata, int iptclen, unsigned char* buffer);
};

class Interpreter {
    protected:
        char buffer[1024];
    public:
        Interpreter () {}
        virtual std::string toString (Tag* t) { return ""; }
        virtual void fromString (Tag* t, const std::string& value) {}
};

class StdInterpreter : public Interpreter {
    public:
        StdInterpreter () {}
        virtual std::string toString (Tag* t) {
            t->toString (buffer);
            return std::string (buffer);
        }
        virtual void fromString (Tag* t, const std::string& value) {
            if (t->getType()==SHORT || t->getType()==LONG)
                t->fromInt (atoi(value.c_str()));
            else
                t->fromString (value.c_str());
        }
};
extern StdInterpreter stdInterpreter; 
class ChoiceInterpreter : public Interpreter {
    protected:
        std::map<int,std::string> choices;
    public:
        ChoiceInterpreter () {};
        virtual std::string toString (Tag* t) {
            std::map<int,std::string>::iterator r = choices.find (t->toInt());
            if (r!=choices.end())
                return r->second;
            else {
                t->toString (buffer);
                return std::string (buffer);
            }
        }
};

inline unsigned short sget2 (unsigned char *s, ByteOrder order);
inline int sget4 (unsigned char *s, ByteOrder order);
inline unsigned short get2 (FILE* f, ByteOrder order);
inline int get4 (FILE* f, ByteOrder order);
inline void sset2 (unsigned short v, unsigned char *s, ByteOrder order);
inline void sset4 (int v, unsigned char *s, ByteOrder order);
inline float int_to_float (int i);
inline short int int2_to_signed (short unsigned int i);


extern const TagAttrib exifAttribs[];
extern const TagAttrib gpsAttribs[];
extern const TagAttrib iopAttribs[];
extern const TagAttrib ifdAttribs[];
extern const TagAttrib nikon2Attribs[];
extern const TagAttrib nikon3Attribs[];
extern const TagAttrib canonAttribs[];
extern const TagAttrib pentaxAttribs[];
extern const TagAttrib fujiAttribs[];
extern const TagAttrib minoltaAttribs[];
extern const TagAttrib sonyAttribs[];
extern const TagAttrib olympusAttribs[];
};
#endif

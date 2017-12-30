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
#include <iostream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <sstream>
#include <stdint.h>
#include <tiff.h>

#include <glib/gstdio.h>
#include <glib/gunicode.h>

#include "rtexif.h"

#include "../rtgui/cacheimagedata.h"
#include "../rtgui/version.h"
#include "../rtgui/ppversion.h"

// see end of ExifManager::parse(bool, bool)
#define PRINT_METADATA_TREE 0

using namespace std;

namespace rtexif
{

Interpreter stdInterpreter;

//--------------- class TagDirectory ------------------------------------------
// this class is a collection (an array) of tags
//-----------------------------------------------------------------------------

TagDirectory::TagDirectory ()
    : attribs (ifdAttribs), order (HOSTORDER), parent (nullptr) {}

TagDirectory::TagDirectory (TagDirectory* p, const TagAttrib* ta, ByteOrder border)
    : attribs (ta), order (border), parent (p) {}

TagDirectory::TagDirectory (TagDirectory* p, FILE* f, int base, const TagAttrib* ta, ByteOrder border, bool skipIgnored)
    : attribs (ta), order (border), parent (p)
{

    int numOfTags = get2 (f, order);

    if (numOfTags <= 0 || numOfTags > 1000) { // KodakIfd has lots of tags, thus 1000 as the limit
        return;
    }

    bool thumbdescr = false;

    for (int i = 0; i < numOfTags; i++) {

        Tag* newTag = new Tag (this, f, base);

        // filter out tags with unknown type
        if ((int)newTag->getType() == 0) {
            delete newTag;
            continue;
        }

        if (skipIgnored) {
            int id = newTag->getID();

            // detect and possibly ignore tags of directories belonging to the embedded thumbnail image
            if (attribs == ifdAttribs && id == TIFFTAG_SUBFILETYPE && newTag->toInt() != 0) {
                thumbdescr = true;
            }

            const TagAttrib* attrib = getAttrib (id);

            if (!attrib || attrib->ignore == 1 || (thumbdescr && attrib->ignore == 2)) {
                delete newTag;
            } else {
                addTag (newTag);
            }
        } else {
            addTag (newTag);
        }
    }
}

TagDirectory::~TagDirectory ()
{

    for (size_t i = 0; i < tags.size(); i++) {
        delete tags[i];
    }
}

class CompareTags
{
public:
    int operator() (Tag* const& a, Tag* const& b) const
    {
        return a->getID() < b->getID();
    }
};

void TagDirectory::sort ()
{

    std::sort (tags.begin(), tags.end(), CompareTags());

    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->isDirectory())
            for (int j = 0; tags[i]->getDirectory (j); j++) {
                tags[i]->getDirectory (j)->sort ();
            }
}
TagDirectory*  TagDirectory::getRoot()
{
    if (parent) {
        return parent->getRoot();
    } else {
        return this;
    }
}

const TagAttrib* TagDirectory::getAttrib (int id)
{

    if (attribs)
        for (int i = 0; attribs[i].ignore != -1; i++)
            if (attribs[i].ID == id) {
                return &attribs[i];
            }

    return nullptr;
}

const TagAttrib* TagDirectory::getAttrib (const char* name)
{

    if (attribs)
        for (int i = 0; attribs[i].ignore != -1; i++)
            if (!strcmp (attribs[i].name, name)) {
                return &attribs[i];
            }

    return nullptr;
}

const TagAttrib* TagDirectory::getAttribP (const char* name)
{

    if (attribs)
        for (int i = 0; attribs[i].ignore != -1; i++) {
            // Yeah, self made comparison!
            const char *n = name;
            const char *a = attribs[i].name;

            while (*n && *a && *n == *a) {
                n++;
                a++;
            };

            if (!*a && (!*n || *n == '/')) {
                // we reached the end of the subpart of name and the end of attribs->name, so they match
                if (*n == '/') {
                    Tag* tag = getTag (attribs[i].ID);
                    TagDirectory *tagDir;

                    if (attribs[i].subdirAttribs && tag && (tagDir = tag->getDirectory())) {
                        return tagDir->getAttribP (n + 1);
                    } else {
                        return nullptr;
                    }
                } else {
                    return &attribs[i];
                }
            }
        }

    return nullptr;
}

void TagDirectory::printAll (unsigned int level) const
{

    // set the spacer prefix string
    char prefixStr[level * 4 + 1];
    unsigned int i;

    for (i = 0; i < level * 4; i++) {
        prefixStr[i] = ' ';
    }

    prefixStr[i] = '\0';

    // recursively iterate over the tag list
    for (size_t i = 0; i < tags.size(); i++) {
        std::string name = tags[i]->nameToString ();

        TagDirectory* currTagDir;
        if (tags[i]->isDirectory()) {
            for (int j = 0; (currTagDir = tags[i]->getDirectory (j)) != nullptr; j++) {
                printf ("%s+-- DIRECTORY %s[%d]:\n", prefixStr, name.c_str(), j);
                currTagDir->printAll (level + 1);
            }
        } else {
            printf ("%s- %s\n", prefixStr, name.c_str());
        }
    }
}

/** @brief Dump the TagDirectory and its sub-directories to the file 'fname'
 *
 * This method has been created to dump the metadata for the Custom Profile Builders.
 * It contains an [RT General] section to communicate some parameters, then the TagDirectory follows.
 *
 * The key is composed as follow: "010F_Make", i.e. "tag number or ID _ tag name"
 * Entries like:
 *
 * 927C_MakerNotesSony=$subdir
 *
 * indicates that this tag refer to a sub-directory. RT's Keywords begins with $, where & is the first char of the value.
 * $subdir is the only keyword so far.
 *
 * You'll have then to check for the [EXIF/927C_MakerNotesSony] section, given that the root section
 * is named [EXIF].
 *
 * WARNING: Some string will be sanitized, i.e. the new line char will be replaced by "\n". You'll
 * have to check for this escape string if you want a correct display of the value, but your KeyFile module
 * will most likely handle that automatically for you.
 *
 * @param commFNname Absolute path of the temporary communication file's name
 * @param commFNname Absolute path of the image's file name
 * @param commFNname Absolute path of the output profiles's file name
 * @param defaultPParams absolute or relative path (to the application's folder) of the default ProcParams to use
 * @param cfs pointer to a CacheImageData object that will contain common values
 * @param flagMode will tell whether the Custom Profile Builder is called for on flagging event or for real development
 * @param keyfile The KeyFile object to dump to. Has to be NULL (default value) on first call!
 * @param tagDirName Name of the current TagDirectory (full path, i.e. "EXIF/MakerNotes/LensInfo"). Can be empty on first call, "EXIF" will then be used
 *
 * @return True if everything went fine, false otherwise
 */
bool TagDirectory::CPBDump (const Glib::ustring &commFName, const Glib::ustring &imageFName, const Glib::ustring &profileFName, const Glib::ustring &defaultPParams,
                            const CacheImageData* cfs, const bool flagMode, Glib::KeyFile *keyFile, Glib::ustring tagDirName) const
{
    const auto kf = keyFile ? keyFile : new Glib::KeyFile;

    if (!kf) {
        return false;
    }

    if (!keyFile || tagDirName.empty()) {
        tagDirName = "EXIF";
    }

    std::vector<const TagDirectory *> tagDirList;
    std::vector<Glib::ustring> tagDirPaths;

    FILE *f = nullptr;

    if (!keyFile) {
        // open the file in write mode
        f = g_fopen (commFName.c_str (), "wt");

        if (f == nullptr) {
            printf ("TagDirectory::keyFileDump(\"%s\") >>> Error: unable to open file with write access!\n", commFName.c_str());
            delete kf;
            return false;
        }

        try {

            kf->set_string ("RT General", "CachePath", options.cacheBaseDir);
            kf->set_string ("RT General", "AppVersion", RTVERSION);
            kf->set_integer ("RT General", "ProcParamsVersion", PPVERSION);
            kf->set_string ("RT General", "ImageFileName", imageFName);
            kf->set_string ("RT General", "OutputProfileFileName", profileFName);
            kf->set_string ("RT General", "DefaultProcParams", defaultPParams);
            kf->set_boolean ("RT General", "FlaggingMode", flagMode);

            kf->set_integer ("Common Data", "FrameCount", cfs->frameCount);
            kf->set_integer ("Common Data", "SampleFormat", cfs->sampleFormat);
            kf->set_boolean ("Common Data", "IsHDR", cfs->isHDR);
            kf->set_boolean ("Common Data", "IsPixelShift", cfs->isPixelShift);
            kf->set_double ("Common Data", "FNumber", cfs->fnumber);
            kf->set_double ("Common Data", "Shutter", cfs->shutter);
            kf->set_double ("Common Data", "FocalLength", cfs->focalLen);
            kf->set_integer ("Common Data", "ISO", cfs->iso);
            kf->set_string ("Common Data", "Lens", cfs->lens);
            kf->set_string ("Common Data", "Make", cfs->camMake);
            kf->set_string ("Common Data", "Model", cfs->camModel);

        } catch (Glib::KeyFileError&) {}
    }

    // recursively iterate over the tag list
    for (size_t i = 0; i < tags.size(); i++) {
        std::string tagName = tags[i]->nameToString ();

        if (tags[i]->isDirectory())
            for (int j = 0; tags[i]->getDirectory (j); j++) {
                // Accumulating the TagDirectories to dump later
                tagDirPaths.push_back ( Glib::ustring ( tagDirName + "/" + getDumpKey (tags[i]->getID(), tagName) ) );
                tagDirList.push_back (tags[i]->getDirectory (j));

                try {
                    kf->set_string (tagDirName, getDumpKey (tags[i]->getID(), tagName), "$subdir");
                } catch (Glib::KeyFileError&) {}
            } else {
            try {
                kf->set_string (tagDirName, getDumpKey (tags[i]->getID(), tagName), tags[i]->valueToString());
            } catch (Glib::KeyFileError&) {}
        }
    }

    // dumping the sub-directories
    for (size_t i = 0; i < tagDirList.size(); i++) {
        tagDirList.at (i)->CPBDump (commFName, imageFName, profileFName, defaultPParams, cfs, flagMode, kf, tagDirPaths.at (i));
    }

    if (!keyFile) {
        try {
            fprintf (f, "%s", kf->to_data().c_str());
        } catch (Glib::KeyFileError&) {}

        fclose (f);
        delete kf;
    }

    return true;
}

Glib::ustring TagDirectory::getDumpKey (int tagID, const Glib::ustring &tagName)
{
    Glib::ustring key;

    if (options.CPBKeys == CPBKT_TID || options.CPBKeys == CPBKT_TID_NAME) {
        key = Glib::ustring (Glib::ustring::format (std::fixed, std::hex, std::setfill (L'0'), std::setw (4), tagID));
    }

    if (options.CPBKeys == CPBKT_TID_NAME) {
        key += Glib::ustring ("_");
    }

    if (options.CPBKeys == CPBKT_TID_NAME || options.CPBKeys == CPBKT_NAME) {
        key += Glib::ustring (tagName);
    }

    return key;
}
void TagDirectory::addTag (Tag* tag)
{

    // look up if it already exists:
    if (getTag (tag->getID())) {
        delete tag;
    } else {
        tags.push_back (tag);
    }
}

void TagDirectory::addTagFront (Tag* tag)
{

    // look up if it already exists:
    if (getTag (tag->getID())) {
        delete tag;
    } else {
        tags.insert (tags.begin(), tag);
    }
}

void TagDirectory::replaceTag (Tag* tag)
{

    // look up if it already exists:
    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->getID() == tag->getID()) {
            delete tags[i];
            tags[i] = tag;
            return;
        }

    tags.push_back (tag);
}

Tag* TagDirectory::getTag (int ID) const
{

    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->getID() == ID) {
            return tags[i];
        }

    return nullptr;
}

Tag* TagDirectory::getTag (const char* name) const
{

    if (attribs) {
        for (int i = 0; attribs[i].ignore != -1; i++)
            if (!strcmp (attribs[i].name, name)) {
                return getTag (attribs[i].ID);
            }
    }

    return nullptr;
}

Tag* TagDirectory::getTagP (const char* name) const
{

    if (attribs)
        for (int i = 0; attribs[i].ignore != -1; i++) {
            // Yeah, self made comparison!
            const char *n = name;
            const char *a = attribs[i].name;

            while (*n && *a && *n == *a) {
                n++;
                a++;
            };

            if (!*a && (!*n || *n == '/')) {
                // we reached the end of the subpart of name and the end of attribs->name, so they match
                if (*n == '/') {
                    Tag* tag = getTag (attribs[i].ID);
                    TagDirectory *tagDir;

                    if (attribs[i].subdirAttribs && tag && (tagDir = tag->getDirectory())) {
                        return tagDir->getTagP (n + 1);
                    } else {
                        return nullptr;
                    }
                } else {
                    return getTag (attribs[i].ID);
                }
            }
        }

    return nullptr;
}

Tag* TagDirectory::findTag (const char* name, bool lookUpward) const
{
    Tag* t = getTag(name);
    if (t) {
        return t;
    }

    for (auto tag : tags) {
        if (tag->isDirectory()) {
            TagDirectory *dir;
            int i = 0;
            while ((dir = tag->getDirectory(i)) != nullptr) {
                TagDirectory *dir = tag->getDirectory();
                Tag* t = dir->findTag (name);

                if (t) {
                    return t;
                }
                ++i;
            }
        }
    }

    if (lookUpward && parent) {
        Tag* t = parent->findTagUpward(name);

        if (t) {
            return t;
        }
    }

    return nullptr;
}

std::vector<const Tag*> TagDirectory::findTags (int ID)
{

    std::vector<const Tag*> tagList;

    //assuming that an entry can only exist once
    Tag* t = getTag(ID);
    if (t) {
        tagList.push_back(t);
    }

    for (auto tag : tags) {
        if (tag->isDirectory()) {
            TagDirectory *dir;
            int i = 0;
            while ((dir = tag->getDirectory(i)) != nullptr) {
                std::vector<const Tag*> subTagList = dir->findTags (ID);

                if (!subTagList.empty()) {
                    // concatenating the 2 vectors
                    // not really optimal in a memory efficiency pov
                    for (auto tag2 : subTagList) {
                        tagList.push_back(tag2);
                    }
                }
                ++i;
            }
        }
    }

    return tagList;
}

std::vector<const Tag*> TagDirectory::findTags (const char* name)
{

    std::vector<const Tag*> tagList;

    //assuming that an entry can only exist once
    Tag* t = getTag(name);
    if (t) {
        tagList.push_back(t);
    }

    for (auto tag : tags) {
        if (tag->isDirectory()) {
            TagDirectory *dir;
            int i = 0;
            while ((dir = tag->getDirectory(i)) != nullptr) {
                std::vector<const Tag*> subTagList = dir->findTags (name);

                if (!subTagList.empty()) {
                    // concatenating the 2 vectors
                    // not really optimal in a memory efficiency pov, but adding 10 items should be a maximum
                    for (auto tag2 : subTagList) {
                        tagList.push_back(tag2);
                    }
                }
                ++i;
            }
        }
    }

    return tagList;
}


Tag* TagDirectory::findTagUpward (const char* name) const
{
    Tag* t = findTag(name);
    if (t) {
        return t;
    }

    if (parent) {
        Tag* t = parent->findTagUpward(name);

        if (t) {
            return t;
        }
    }

    return nullptr;
}


// Searches a simple value, as either attribute or element
// only for simple values, not for entries with special chars or free text
bool TagDirectory::getXMPTagValue (const char* name, char* value) const
{
    *value = 0;

    if (!getTag ("ApplicationNotes")) {
        return false;
    }

    char *sXMP = (char*)getTag ("ApplicationNotes")->getValue();

    // Check for full word
    char *pos = sXMP;

    bool found = false;

    do {
        pos = strstr (pos, name);

        if (pos) {
            char nextChar = * (pos + strlen (name));

            if (nextChar == ' ' || nextChar == '>' || nextChar == '=') {
                found = true;
            } else {
                pos += strlen (name);
            }
        }
    } while (pos && !found);

    if (!found) {
        return false;
    }

    char *posTag = strchr (pos, '>');
    char *posAttr = strchr (pos, '"');

    if (!posTag && !posAttr) {
        return false;
    }

    if (posTag && (!posAttr || posTag < posAttr)) {
        // Tag
        pos = strchr (posTag + 1, '<');
        strncpy (value, posTag + 1, pos - posTag - 1);
        value[pos - posTag - 1] = 0;
        return true;
    } else if (posAttr && (!posTag || posAttr < posTag)) {
        // Attribute
        pos = strchr (posAttr + 1, '"');
        strncpy (value, posAttr + 1, pos - posAttr - 1);
        value[pos - posAttr - 1] = 0;
        return true;
    } else {
        return false;
    }
}

void TagDirectory::keepTag (int ID)
{
    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->getID() == ID) {
            tags[i]->setKeep (true);
        }
}

int TagDirectory::calculateSize ()
{

    int size = 2; // space to store the number of tags

    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->getKeep()) {
            size += 12 + tags[i]->calculateSize ();
        }

    size += 4; // next ifd pointer
    return size;
}

TagDirectory* TagDirectory::clone (TagDirectory* parent)
{

    TagDirectory* td = new TagDirectory (parent, attribs, order);

    for (size_t i = 0; i < tags.size(); i++) {
        td->tags.push_back (tags[i]->clone (td));
    }

    return td;
}

int TagDirectory::write (int start, unsigned char* buffer)
{

    int size = calculateSize ();
    int tagnum = 0;
    int nondirspace = 0;

    for (size_t i = 0; i < tags.size(); i++)
        if (tags[i]->getKeep()) {
            tagnum++;

            if (!tags[i]->isDirectory()) {
                nondirspace += tags[i]->calculateSize();
            }
        }

    int nextValOffs = start + 2 + tagnum * 12 + 4;
    int nextDirOffs = nextValOffs + nondirspace;
    int pos = start;
    sset2 (tagnum, buffer + start, order);
    pos += 2;
    int maxPos = start + size;

    for (size_t i = 0; i < tags.size(); i++) {
        if (tags[i]->getKeep()) {
            if (!tags[i]->isDirectory()) {
                nextValOffs = tags[i]->write (pos, nextValOffs, buffer);    // pos: where to put the tag, dataoffset: the place where the value can be put. return: next data offset
            } else {
                nextDirOffs = tags[i]->write (pos, nextDirOffs, buffer);    // pos: where to put the tag, dataoffset: the place where the value can be put. return: next data offset
            }

            pos += 12;
        }
    }

    sset4 (0, buffer + pos, order);
    return maxPos;
}

void TagDirectory::applyChange (std::string name, Glib::ustring value)
{

    std::string::size_type dp = name.find_first_of ('.');
    std::string fseg = name.substr (0, dp);

    // this is a final segment: apply change
    if (dp == std::string::npos) {

        Tag* t = nullptr;

        for (size_t i = 0; i < tags.size(); i++)
            if (tags[i]->nameToString() == fseg) {
                t = tags[i];
                break;
            }

        if (value == "#keep" && t) {
            t->setKeep (true);
        } else if (value == "#delete" && t) {
            t->setKeep (false);
        } else if (t && !t->isDirectory()) {
            if (name == "UserComment") {
                // UserComment can be Unicode
                t->userCommentFromString (value);
            } else {
                t->valueFromString (value);
            }
        } else {
            const TagAttrib* attrib = nullptr;

            for (int i = 0; attribs[i].ignore != -1; i++)
                if (!strcmp (attribs[i].name, fseg.c_str())) {
                    attrib = &attribs[i];
                    break;
                }

            if (attrib) {
                Tag* nt = new Tag (this, attrib);
                if (name == "UserComment") {
                    // UserComment can be Unicode
                    nt->initUserComment (value);
                } else {
                    nt->initString (value.c_str());
                }
                addTag (nt);
            }
        }
    }
    // this is a subdirectory
    else {
        // try to find it
        std::string::size_type dp1 = fseg.find_first_of ('[');
        std::string basename = fseg.substr (0, dp1);
        Tag* t = nullptr;
        int dirnum = -1;

        for (size_t i = 0; i < tags.size(); i++)
            if (tags[i]->isDirectory()) {
                for (int j = 0; tags[i]->getDirectory (j); j++) {
                    if (tags[i]->nameToString (j) == fseg) {
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

        if (!t && value != "#keep" && value != "#delete") {
            const TagAttrib* attrib = nullptr;

            for (int i = 0; attribs[i].ignore != -1; i++)
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

        if (t && dirnum >= 0) {
            t->getDirectory (dirnum)->applyChange (name.substr (dp + 1, std::string::npos), value);
        }
    }
}

TagDirectoryTable::TagDirectoryTable ()
    : values (nullptr), zeroOffset (0), valuesSize (0), defaultType (INVALID)
{
}

TagDirectoryTable::TagDirectoryTable (TagDirectory* p, unsigned char *v, int memsize, int offs, TagType type, const TagAttrib* ta, ByteOrder border)
    : TagDirectory (p, ta, border), zeroOffset (offs), valuesSize (memsize), defaultType ( type )
{
    values = new unsigned char[valuesSize];
    memcpy (values, v, valuesSize);

    // Security ; will avoid to read above the buffer limit if the RT's tagDirectoryTable is longer that what's in the file
    int count = valuesSize / getTypeSize (type);

    for (const TagAttrib* tattr = ta; tattr->ignore != -1 && tattr->ID < count; ++tattr) {
        Tag* newTag = new Tag (this, tattr, (values + zeroOffset + tattr->ID * getTypeSize (type)), tattr->type == AUTO ? type : tattr->type);
        tags.push_back (newTag); // Here we can insert more tag in the same offset because of bitfield meaning
    }
}

TagDirectoryTable::TagDirectoryTable (TagDirectory* p, FILE* f, int memsize, int offs, TagType type, const TagAttrib* ta, ByteOrder border)
    : TagDirectory (p, ta, border), zeroOffset (offs), valuesSize (memsize), defaultType ( type )
{
    values = new unsigned char[valuesSize];
    fread (values, 1, valuesSize, f);

    // Security ; will avoid to read above the buffer limit if the RT's tagDirectoryTable is longer that what's in the file
    int count = valuesSize / getTypeSize (type);

    for (const TagAttrib* tattr = ta; tattr->ignore != -1 && tattr->ID < count; ++tattr) {
        Tag* newTag = new Tag (this, tattr, (values + zeroOffset + tattr->ID * getTypeSize (type)), tattr->type == AUTO ? type : tattr->type);
        tags.push_back (newTag); // Here we can insert more tag in the same offset because of bitfield meaning
    }
}
TagDirectory* TagDirectoryTable::clone (TagDirectory* parent)
{

    TagDirectory* td = new TagDirectoryTable (parent, values, valuesSize, zeroOffset, defaultType, attribs, order);
    return td;
}

TagDirectoryTable::~TagDirectoryTable()
{
    if (values) {
        delete [] values;
    }
}
int TagDirectoryTable::calculateSize ()
{
    return valuesSize;
}

int TagDirectoryTable::write (int start, unsigned char* buffer)
{
    if ( values && valuesSize) {
        memcpy (buffer + start, values, valuesSize);
        return start + valuesSize;
    } else {
        return start;
    }
}

//--------------- class Tag ---------------------------------------------------
// this class represents a tag stored in the directory
//-----------------------------------------------------------------------------

Tag::Tag (TagDirectory* p, FILE* f, int base)
    : type (INVALID), count (0), value (nullptr), allocOwnMemory (true), attrib (nullptr), parent (p), directory (nullptr)
{

    ByteOrder order = getOrder();

    tag   = get2 (f, order);
    type  = (TagType)get2 (f, order);
    count = get4 (f, order);

    if (!count) {
        count = 1;
    }

    makerNoteKind = NOMK;
    keep = false;

    // filter out invalid tags
    // note the large count is to be able to pass LeafData ASCII tag which can be up to almost 10 megabytes,
    // (only a small part of it will actually be parsed though)
    if ((int)type < 1 || (int)type > 14 || count > 10 * 1024 * 1024) {
        type = INVALID;
        valuesize = 0;
        return;
    }

    // store next Tag's position in file
    int save = ftell (f) + 4;

    // load value field (possibly seek before)
    valuesize = count * getTypeSize (type);

    if (valuesize > 4) {
        fseek (f, get4 (f, getOrder()) + base, SEEK_SET);
    }

    attrib = parent->getAttrib (tag);

    if (attrib && (attrib->action == AC_WRITE || attrib->action == AC_NEW)) {
        keep = true;
    }

    if ( tag == 0xc634 ) { // DNGPrivateData
        int currPos = ftell (f);
        const int buffersize = 32;
        char buffer[buffersize], *p = buffer;

        while ( fread (p, 1, 1, f ) && *p != 0 && p - buffer < buffersize - 1 ) {
            p++;
        }

        *p = 0;

        if ( !strncmp (buffer, "Adobe", 5) ) {
            fread (buffer, 1, 14, f );

            if ( !strncmp ( buffer, "MakN", 4) ) {
                ByteOrder bom = ((buffer[8] == 'M' && buffer[9] == 'M') ? MOTOROLA : INTEL) ;
                Tag* tmake = parent->getRoot()->findTag ("Make");
                std::string make ( tmake ? tmake->valueToString() : "");
                int save = ftell (f);
                int originalOffset = sget4 ( (unsigned char*)&buffer[10], ( make.find ("SONY") != std::string::npos ) || ( make.find ("Canon") != std::string::npos ) || ( make.find ("OLYMPUS") != std::string::npos ) ? MOTOROLA : bom );

                if ( !parseMakerNote (f, save - originalOffset, bom )) {
                    type = INVALID;
                }
            }
        } else if ( !strncmp (buffer, "PENTAX", 6) ) {
            makerNoteKind = HEADERIFD;
            fread (buffer, 1, 2, f);
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, currPos, pentaxAttribs, strncmp (buffer, "MM", 2) ? INTEL : MOTOROLA);
            directory[1] = nullptr;
        } else
            /* SONY uses this tag to write hidden info and pointer to private encrypted tags
            {
             unsigned offset =sget4((unsigned char*)buffer, order);
             fseek(f,offset,SEEK_SET);
             makerNoteKind = TABLESUBDIR;
             directory = new TagDirectory*[2];
             directory[0] = new TagDirectory (parent, f, base, sonyDNGMakerNote, order);
             directory[1] = NULL;
             fseek (f, save, SEEK_SET);
             return;
            }*/
        {
            type = INVALID;
        }
    }

    if (tag == 0x002e) { // location of the embedded preview image in raw files of Panasonic cameras
        ExifManager eManager(f, nullptr, true);
        eManager.parseJPEG(ftell(f)); // try to parse the exif data from the preview image

        if (eManager.roots.size()) {
            const TagDirectory* const previewdir = eManager.roots.at(0);
            if (previewdir->getTag ("Exif")) {
                if (previewdir->getTag ("Make")) {
                    if (previewdir->getTag ("Make")->valueToString() == "Panasonic") { // "make" is not yet available here, so get it from the preview tags to assure we're doing the right thing
                        Tag* t = new Tag (parent->getRoot(), lookupAttrib (ifdAttribs, "Exif")); // replace raw exif with preview exif assuming there are the same
                        t->initSubDir (previewdir->getTag ("Exif")->getDirectory());
                        parent->getRoot()->addTag (t);
                    }
                }
            }
        }
    }

    // if this tag is the makernote, it needs special treatment (brand specific parsing)
    if (tag == 0x927C && attrib && !strcmp (attrib->name, "MakerNote") ) {
        if ( !parseMakerNote (f, base, order )) {
            type = INVALID;
            fseek (f, save, SEEK_SET);
            return;
        }
    } else if (attrib && attrib->subdirAttribs) {
        // Some subdirs are specific of maker and model
        char make[128], model[128];
        make[0] = 0;
        model[0] = 0;
        Tag* tmake = parent->getRoot()->getTag ("Make");

        if (tmake) {
            tmake->toString (make);
        }

        Tag* tmodel = parent->getRoot()->getTag ("Model");

        if (tmodel) {
            tmodel->toString (model);
        }

        if (!strncmp (make, "SONY", 4)) {
            switch ( tag ) {
                case 0x0010:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;

                    if (count == 15360) {
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, sonyCameraInfoAttribs, order);
                    } else {
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, sonyCameraInfo2Attribs, order);
                    }

                    break;

                case 0x0114:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;

                    if (count == 280 || count == 364) {
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, SHORT, sonyCameraSettingsAttribs, MOTOROLA);
                    } else if (count == 332) {
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, SHORT, sonyCameraSettingsAttribs2, MOTOROLA);
                    } else if (count == 1536 || count == 2048) {
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, sonyCameraSettingsAttribs3, INTEL);
                    } else {
                        // Unknown CameraSettings
                        delete [] directory;
                        directory = nullptr;
                        type = INVALID;
                    }

                    makerNoteKind = directory ? TABLESUBDIR : NOMK;
                    break;

                case 0x9405:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, SHORT, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                default:
                    goto defsubdirs;
            }
        } else if ((!strncmp (make, "PENTAX", 6)) || (!strncmp (make, "RICOH", 5) && !strncmp (model, "PENTAX", 6))) { // Either the former Pentax brand or the RICOH brand + PENTAX model"
            switch ( tag ) {
                case 0x007d:
                case 0x0205:
                case 0x0208:
                case 0x0216:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                case 0x0215:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, LONG, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                case 0x005c:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;

                    if (count == 4) {     // SRInfo
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, pentaxSRInfoAttribs, order);
                    } else if (count == 2) { // SRInfo2
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, pentaxSRInfo2Attribs, order);
                    } else {
                        // Unknown SRInfo
                        delete [] directory;
                        directory = nullptr;
                        type = INVALID;
                    }

                    makerNoteKind = directory ? TABLESUBDIR : NOMK;
                    break;

                case 0x0206:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;

                    if (count == 21) {     // AEInfo2
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, pentaxAEInfo2Attribs, order);
                    } else if (count == 48) { // AEInfo3
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, pentaxAEInfo3Attribs, order);
                    } else if (count <= 25) { // AEInfo
                        directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, pentaxAEInfoAttribs, order);
                    } else {
                        // Unknown AEInfo
                        delete [] directory;
                        directory = nullptr;
                        type = INVALID;
                    }

                    makerNoteKind = directory ? TABLESUBDIR : NOMK;
                    break;

                case 0x0207: {
                    // There are 2 format pentaxLensDataAttribs
                    int offsetFirst = 4;  // LensInfo2

                    if ( strstr (model, "*ist") || strstr (model, "GX-1") || strstr (model, "K200D") || (strstr (model, "K100D") && !strstr (model, "K100D Super")) || strstr (model, "K110D") || strstr (model, "645Z")) {
                        offsetFirst = 3;    // LensInfo
                    } else if ( strstr (model, "645D") ) {
                        offsetFirst = 13;    // LensInfo3
                    } else if ( strstr (model, "K-01") || strstr (model, "K-30") || strstr (model, "K-50")) {
                        offsetFirst = 15;    // LensInfo5
                    } else if ( strstr (model, "K-5") || strstr (model, "K-r") ) {
                        offsetFirst = 12;    // LensInfo4
                    } else if (!strncmp (make, "RICOH", 5)) { // all PENTAX camera model produced under the RICOH era uses LensInfo5, for now...
                        offsetFirst = 15;  // LensInfo5 too
                    }

                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, offsetFirst, BYTE, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                }
                break;

                case 0x0239:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                default:
                    goto defsubdirs;
            }
        } else if (!strncmp (make, "Canon", 5)) {
            switch ( tag ) {
                case 0x0001:
                case 0x0002:
                case 0x0004:
                case 0x0005:
                case 0x0093:
                case 0x0098:
                case 0x00a0:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, SSHORT, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                case 0x009a:
                case 0x4013:
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, LONG, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;

                default:
                    goto defsubdirs;
            }
        } else if (!strncmp (make, "NIKON", 5)) {
            switch (tag) {
                case 0x0025: {
                    directory = new TagDirectory*[2];
                    directory[1] = nullptr;
                    directory[0] = new TagDirectoryTable (parent, f, valuesize, 0, BYTE, attrib->subdirAttribs, order);
                    makerNoteKind = TABLESUBDIR;
                    break;
                }

                default:
                    goto defsubdirs;
            }
        } else if (type == UNDEFINED) {
            count = 1;
            type = LONG;
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, base, attrib->subdirAttribs, order);
            directory[1] = nullptr;
        } else {
            goto defsubdirs;
        }
    } else {
        // read value
        value = new unsigned char [valuesize + 1];
        fread (value, 1, valuesize, f);
        value[valuesize] = '\0';
    }

    // seek back to the saved position
    fseek (f, save, SEEK_SET);
    return;

defsubdirs:
    // read value
    value = new unsigned char [valuesize];
    fread (value, 1, valuesize, f);

    // count the number of valid subdirs
    int sdcount = count;

    if (sdcount > 0) {
        if (parent->getAttribTable() == olympusAttribs) {
            sdcount = 1;
        }

        // allocate space
        directory = new TagDirectory*[sdcount + 1];

        // load directories
        for (size_t j = 0, i = 0; j < count; j++, i++) {
            int newpos = base + toInt (j * 4, LONG);
            fseek (f, newpos, SEEK_SET);
            directory[i] = new TagDirectory (parent, f, base, attrib->subdirAttribs, order);
        }

        // set the terminating NULL
        directory[sdcount] = nullptr;
    } else {
        type = INVALID;
    }

    // seek back to the saved position
    fseek (f, save, SEEK_SET);
    return;

}

bool Tag::parseMakerNote (FILE* f, int base, ByteOrder bom )
{
    value = nullptr;
    Tag* tmake = parent->getRoot()->findTag ("Make");
    std::string make ( tmake ? tmake->valueToString() : "");

    Tag* tmodel = parent->getRoot()->findTag ("Model");
    std::string model ( tmodel ? tmodel->valueToString() : "");

    if ( make.find ( "NIKON" ) != std::string::npos ) {
        if ( model.find ("NIKON E700") != std::string::npos ||
                model.find ("NIKON E800") != std::string::npos ||
                model.find ("NIKON E900") != std::string::npos ||
                model.find ("NIKON E900S") != std::string::npos ||
                model.find ("NIKON E910") != std::string::npos ||
                model.find ("NIKON E950") != std::string::npos ) {
            makerNoteKind = HEADERIFD;
            valuesize = 8;
            value = new unsigned char[8];
            fread (value, 1, 8, f);
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, base, nikon2Attribs, bom);
            directory[1] = nullptr;
        } else if ( model.find ("NIKON E990") != std::string::npos ||
                    (model.find ("NIKON D1") != std::string::npos && model.size() > 8 && model.at (8) != '0')) {
            makerNoteKind = IFD;
            directory = new TagDirectory*[2];
            directory[0] = new TagDirectory (parent, f, base, nikon3Attribs, bom);
            directory[1] = nullptr;
        } else {
            // needs refinement! (embedded tiff header parsing)
            makerNoteKind = NIKON3;
            valuesize = 18;
            value = new unsigned char[18];
            int basepos = ftell (f);
            fread (value, 1, 18, f);
            directory = new TagDirectory*[2];
            // byte order for makernotes can be different from exif byte order. We have to get it from makernotes header
            ByteOrder MakerNoteOrder;

            if (value[10] == 'M' && value[11] == 'M') {
                MakerNoteOrder = rtexif::MOTOROLA;
            } else {
                MakerNoteOrder = rtexif::INTEL;
            }

            directory[0] = new TagDirectory (parent, f, basepos + 10, nikon3Attribs, MakerNoteOrder);
            directory[1] = nullptr;
        }
    } else if ( make.find ( "Canon" ) != std::string::npos  ) {
        makerNoteKind = IFD;
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, canonAttribs, bom);
        directory[1] = nullptr;
    } else if ( make.find ( "PENTAX" ) != std::string::npos ) {
        makerNoteKind = HEADERIFD;
        valuesize = 6;
        value = new unsigned char[6];
        fread (value, 1, 6, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, pentaxAttribs, bom);
        directory[1] = nullptr;
    } else if ( (make.find ( "RICOH" ) != std::string::npos ) && (model.find ("PENTAX") != std::string::npos) ) {
        makerNoteKind = HEADERIFD;
        valuesize = 10;
        value = new unsigned char[10];
        fread (value, 1, 10, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, ftell (f) - 10, pentaxAttribs, bom);
        directory[1] = nullptr;
    } else if ( make.find ( "FUJIFILM" ) != std::string::npos ) {
        makerNoteKind = FUJI;
        valuesize = 12;
        value = new unsigned char[12];
        fread (value, 1, 12, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, ftell (f) - 12, fujiAttribs, INTEL);
        directory[1] = nullptr;
    } else if ( make.find ( "KONICA MINOLTA" ) != std::string::npos || make.find ( "Minolta" ) != std::string::npos ) {
        makerNoteKind = IFD;
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, minoltaAttribs, bom);
        directory[1] = nullptr;
    } else if ( make.find ( "SONY" ) != std::string::npos ) {
        valuesize = 12;
        value = new unsigned char[12];
        fread (value, 1, 12, f);

        if (!strncmp ((char*)value, "SONY DSC", 8)) {
            makerNoteKind = HEADERIFD;
        } else {
            makerNoteKind = IFD;
            fseek (f, -12, SEEK_CUR);
        }

        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, sonyAttribs, bom );
        directory[1] = nullptr;
    } else if ( make.find ( "OLYMPUS" ) != std::string::npos ) {
        makerNoteKind = HEADERIFD;
        valuesize = 8;
        value = new unsigned char[12];
        fread (value, 1, 8, f);
        directory = new TagDirectory*[2];
        directory[1] = nullptr;

        if (!strncmp ((char*)value, "OLYMPUS", 7)) {
            makerNoteKind = OLYMPUS2;
            fread (value + 8, 1, 4, f);
            valuesize = 12;
            directory[0] = new TagDirectory (parent, f, ftell (f) - 12, olympusAttribs, value[8] == 'I' ? INTEL : MOTOROLA);
        } else {
            directory[0] = new TagDirectory (parent, f, base, olympusAttribs, bom);
        }
    } else if ( make.find ( "Panasonic" ) != std::string::npos) {
        makerNoteKind = HEADERIFD;
        valuesize = 12;
        value = new unsigned char[12];
        fread (value, 1, 12, f);
        directory = new TagDirectory*[2];
        directory[0] = new TagDirectory (parent, f, base, panasonicAttribs, bom);
        directory[1] = nullptr;
    } else {
        return false;
    }

    return true;
}

Tag* Tag::clone (TagDirectory* parent)
{

    Tag* t = new Tag (parent, attrib);

    t->tag = tag;
    t->type = type;
    t->count = count;
    t->keep = keep;
    t->valuesize = valuesize;

    if (value) {
        t->value = new unsigned char [valuesize];
        memcpy (t->value, value, valuesize);
    } else {
        value = nullptr;
    }

    t->makerNoteKind = makerNoteKind;

    if (directory) {
        int ds = 0;

        for (; directory[ds]; ds++);

        t->directory = new TagDirectory*[ds + 1];

        for (int i = 0; i < ds; i++) {
            t->directory[i] = directory[i]->clone (parent);
        }

        t->directory[ds] = nullptr;
    } else {
        t->directory = nullptr;
    }

    return t;
}

Tag::~Tag ()
{

    // delete value
    if (value && allocOwnMemory) {
        delete [] value;
    }

    // if there are directories behind the tag, delete them
    if (directory) {
        int i = 0;

        while (directory[i]) {
            delete directory[i++];
        }

        delete [] directory;
    }
}

void Tag::setInt (int v, int ofs, TagType astype)
{

    if (astype == SHORT) {
        sset2 (v, value + ofs, getOrder());
    } else if (astype == RATIONAL) {
        sset4 (v, value + ofs, getOrder());
        sset4 (1, value + ofs + 4, getOrder());
    } else {
        sset4 (v, value + ofs, getOrder());
    }
}

void Tag::fromInt (int v)
{

    if (type == SHORT) {
        sset2 (v, value, getOrder());
    } else {
        sset4 (v, value, getOrder());
    }
}

void Tag::fromString (const char* v, int size)
{

    if ( value && allocOwnMemory) {
        delete [] value;
    }

    if (size < 0) {
        valuesize = strlen (v) + 1;
    } else {
        valuesize = size;
    }

    count = valuesize;

    if ( allocOwnMemory ) {
        value = new unsigned char [valuesize];
    }

    if(value) {
        memcpy ((char*)value, v, valuesize);
    }
}

int Tag::toInt (int ofs, TagType astype) const
{
    if (attrib) {
        return attrib->interpreter->toInt (this, ofs, astype);
    }

    int a;

    if (astype == INVALID) {
        astype = type;
    }

    switch (astype) {
        //case SBYTE: return (signed char)(value[ofs]);
        case SBYTE:
            return int ((reinterpret_cast<signed char*> (value))[ofs]);

        case BYTE:
            return value[ofs];

        case ASCII:
            return 0;

        case SSHORT:
            return (int)int2_to_signed (sget2 (value + ofs, getOrder()));

        case SHORT:
            return (int)sget2 (value + ofs, getOrder());

        case SLONG:
        case LONG:
            return (int)sget4 (value + ofs, getOrder());

        case SRATIONAL:
        case RATIONAL:
            a = (int)sget4 (value + ofs + 4, getOrder());
            return a == 0 ? 0 : (int)sget4 (value + ofs, getOrder()) / a;

        case FLOAT:
            return (int)toDouble (ofs);

        case UNDEFINED:
            return 0;

        default:
            return 0; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
    }

    return 0;
}

double Tag::toDouble (int ofs) const
{
    if (attrib) {
        return attrib->interpreter->toDouble (this, ofs);
    }

    union IntFloat {
        uint32_t i;
        float f;
    } conv;

    double ud, dd;

    switch (type) {
        case SBYTE:
            return (double) (int ((reinterpret_cast<signed char*> (value))[ofs]));

        case BYTE:
            return (double) ((int)value[ofs]);

        case ASCII:
            return 0.0;

        case SSHORT:
            return (double)int2_to_signed (sget2 (value + ofs, getOrder()));

        case SHORT:
            return (double) ((int)sget2 (value + ofs, getOrder()));

        case SLONG:
        case LONG:
            return (double) ((int)sget4 (value + ofs, getOrder()));

        case SRATIONAL:
        case RATIONAL:
            ud = (int)sget4 (value + ofs, getOrder());
            dd = (int)sget4 (value + ofs + 4, getOrder());
            return dd == 0. ? 0. : (double)ud / (double)dd;

        case FLOAT:
            conv.i = sget4 (value + ofs, getOrder());
            return conv.f;  // IEEE FLOATs are already C format, they just need a recast

        case UNDEFINED:
            return 0.;

        default:
            return 0.; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
    }

}

/**
 * @brief Create an array of the elements
 */
double *Tag::toDoubleArray (int ofs)
{
    double *values = new double[count];

    for (unsigned int i = 0; i < count; ++i) {
        values[i] = toDouble (ofs + i * getTypeSize (type));
    }

    return values;
}

void Tag::toRational (int& num, int& denom, int ofs)
{

    switch (type) {
        case BYTE:
            num = (int)value[ofs];
            denom = 1;
            break;

        case ASCII:
            num = 0;
            denom = 0;
            break;

        case SSHORT:
        case SHORT:
            num = (int)sget2 (value + ofs, getOrder());
            denom = 1;
            break;

        case SLONG:
        case LONG:
            num = (int)sget4 (value + ofs, getOrder());
            denom = 1;
            break;

        case SRATIONAL:
        case RATIONAL:
            num = (int)sget4 (value + ofs, getOrder());
            denom = (int)sget4 (value + ofs + 4, getOrder());
            break;

        case FLOAT:
            num = 0;
            denom = 0;
            break;

        case UNDEFINED:
            num = 0;
            denom = 0;
            break;

        default:
            num = 0;
            denom = 0;
            break; // Quick fix for missing cases (INVALID, DOUBLE, OLYUNDEF, SUBDIR)
    }
}

void Tag::toString (char* buffer, int ofs)
{

    if (type == UNDEFINED && !directory) {
        bool isstring = true;
        unsigned int i = 0;

        for (i = 0; i + ofs < count && i < 64 && value[i + ofs]; i++)
            if (value[i + ofs] < 32 || value[i + ofs] > 126) {
                isstring = false;
            }

        if (isstring) {
            int j = 0;

            for (i = 0; i + ofs < count && i < 64 && value[i + ofs]; i++) {
                if (value[i + ofs] == '<' || value[i + ofs] == '>') {
                    buffer[j++] = '\\';
                }

                buffer[j++] = value[i + ofs];
            }

            buffer[j++] = 0;
            return;
        }
    } else if (type == ASCII) {
        sprintf (buffer, "%.64s", value + ofs);
        return;
    }

    size_t maxcount = 4;

    if (count < 4) {
        maxcount = count;
    }

    strcpy (buffer, "");

    for (ssize_t i = 0; i < std::min<int>(maxcount, valuesize - ofs); i++) {
        if (i > 0) {
            strcat (buffer, ", ");
        }

        char* b = buffer + strlen (buffer);

        switch (type) {
            case UNDEFINED:
            case BYTE:
                sprintf (b, "%d", value[i + ofs]);
                break;

            case SSHORT:
                sprintf (b, "%d", toInt (2 * i + ofs));
                break;

            case SHORT:
                sprintf (b, "%u", toInt (2 * i + ofs));
                break;

            case SLONG:
                sprintf (b, "%d", toInt (4 * i + ofs));
                break;

            case LONG:
                sprintf (b, "%u", toInt (4 * i + ofs));
                break;

            case SRATIONAL:
            case RATIONAL:
                sprintf (b, "%d/%d", (int)sget4 (value + 8 * i + ofs, getOrder()), (int)sget4 (value + 8 * i + ofs + 4, getOrder()));
                break;

            case FLOAT:
                sprintf (b, "%g", toDouble (8 * i + ofs));
                break;

            default:
                break;
        }
    }

    if (count > maxcount) {
        strcat (buffer, "...");
    }
}

std::string Tag::nameToString (int i)
{

    char buffer[1025];

    if (attrib) {
        strncpy (buffer, attrib->name, 1024);
    } else {
        sprintf (buffer, "0x%x", tag);
    }

    if (i > 0) {
        sprintf (buffer + strlen (buffer) - 1, "[%d]", i);
    }

    return buffer;
}

std::string Tag::valueToString ()
{

    if (attrib && attrib->interpreter) {
        return attrib->interpreter->toString (this);
    } else {
        char buffer[1024];
        toString (buffer);
        return buffer;
    }
}

void Tag::valueFromString (const std::string& value)
{

    if (attrib && attrib->interpreter) {
        attrib->interpreter->fromString (this, value);
    }
}

void Tag::userCommentFromString (const Glib::ustring& text)
{

    if (!allocOwnMemory) {
        return;
    }
    if (value) {
        delete [] value;
        value = nullptr;
    }
    initUserComment(text);
}

int Tag::calculateSize ()
{
    int size = 0;

    if (directory) {
        int j;

        for (j = 0; directory[j]; j++) {
            size += directory[j]->calculateSize ();
        }

        if (j > 1) {
            size += 4 * j;
        }
        if (makerNoteKind != NOMK) {
            count = directory[0]->calculateSize () / getTypeSize (type);
        }
    } else if (valuesize > 4) {
        size += valuesize + (valuesize % 2);    // we align tags to even byte positions
    }


    if (makerNoteKind == NIKON3 || makerNoteKind == OLYMPUS2 || makerNoteKind == FUJI || makerNoteKind == HEADERIFD) {
        size += valuesize;
    }

    return size;
}

int Tag::write (int offs, int dataOffs, unsigned char* buffer)
{

    if ((int)type == 0 || offs > 65500) {
        return dataOffs;
    }

    sset2 (tag, buffer + offs, parent->getOrder());
    offs += 2;
    unsigned short typ = type;
    sset2 (typ, buffer + offs, parent->getOrder());
    offs += 2;
    sset4 (count, buffer + offs, parent->getOrder());
    offs += 4;

    if (!directory) {
        if (valuesize > 4) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            memcpy (buffer + dataOffs, value, valuesize);

            if (valuesize % 2) {
                buffer[dataOffs + valuesize] = 0;    // zero padding required by the exif standard
            }

            return dataOffs + valuesize + (valuesize % 2);
        } else {
            memcpy (buffer + offs, value, valuesize);
            return dataOffs;
        }
    } else {
        if (makerNoteKind == NIKON3) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            memcpy (buffer + dataOffs, value, 18);
            dataOffs += 10;
            dataOffs += directory[0]->write (8, buffer + dataOffs);
            return dataOffs;
        } else if (makerNoteKind == OLYMPUS2 || makerNoteKind == FUJI) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            memcpy (buffer + dataOffs, value, valuesize);
            dataOffs += valuesize + directory[0]->write (valuesize, buffer + dataOffs);
            return dataOffs;
        } else if (makerNoteKind == HEADERIFD) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            memcpy (buffer + dataOffs, value, valuesize);
            dataOffs += valuesize;
            dataOffs += directory[0]->write (dataOffs, buffer);
            return dataOffs;
        } else if ( makerNoteKind == TABLESUBDIR) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            dataOffs = directory[0]->write (dataOffs, buffer);
            return dataOffs;
        } else if (!directory[1]) {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            return directory[0]->write (dataOffs, buffer);
        } else {
            sset4 (dataOffs, buffer + offs, parent->getOrder());
            int linkOffs = dataOffs;

            for (int i = 0; directory[i]; i++) {
                dataOffs += 4;
            }

            for (int i = 0; directory[i]; i++) {
                sset4 (dataOffs, buffer + linkOffs, parent->getOrder());
                linkOffs += 4;
                dataOffs = directory[i]->write (dataOffs, buffer);
            }

            return dataOffs;
        }
    }
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr)
    : tag (attr ? attr->ID : -1), type (INVALID), count (0), value (nullptr), valuesize (0), keep (true), allocOwnMemory (true), attrib (attr), parent (p), directory (nullptr), makerNoteKind (NOMK)
{
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, int data, TagType t)
    : tag (attr ? attr->ID : -1), type (t), count (1), value (nullptr), valuesize (0), keep (true), allocOwnMemory (true), attrib (attr), parent (p), directory (nullptr), makerNoteKind (NOMK)
{

    initInt (data, t);
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, unsigned char *data, TagType t)
    : tag (attr ? attr->ID : -1), type (t), count (1), value (nullptr), valuesize (0), keep (true), allocOwnMemory (false), attrib (attr), parent (p), directory (nullptr), makerNoteKind (NOMK)
{

    initType (data, t);
}

Tag::Tag (TagDirectory* p, const TagAttrib* attr, const char* text)
    : tag (attr ? attr->ID : -1), type (ASCII), count (1), value (nullptr), valuesize (0), keep (true), allocOwnMemory (true), attrib (attr), parent (p), directory (nullptr), makerNoteKind (NOMK)
{

    initString (text);
}

void Tag::initType (unsigned char *data, TagType type)
{
    valuesize = getTypeSize (type);

    if ( allocOwnMemory ) {
        value = new unsigned char[valuesize];
        memcpy ((char*)value, data, valuesize);
    } else {
        value = data;
    }
}

void Tag::initInt (int data, TagType t, int cnt)
{

    type = t;

    if (t == LONG) {
        valuesize = 4;
    } else if (t == SHORT) {
        valuesize = 2;
    } else if (t == BYTE) {
        valuesize = 1;
    } else if (t == RATIONAL) {
        valuesize = 8;
    }

    count = cnt;
    valuesize *= count;
    value = new unsigned char[valuesize];
    setInt (data, 0, t);
}

void Tag::initUserComment (const Glib::ustring &text)
{
    type = UNDEFINED;
    if (text.is_ascii()) {
        count = 8 + strlen (text.c_str());
        valuesize = count;
        value = new unsigned char[valuesize];
        memcpy((char*)value, "ASCII\0\0\0", 8);
        memcpy((char*)value + 8, text.c_str(), valuesize - 8);
    } else {
        wchar_t *commentStr = (wchar_t*)g_utf8_to_utf16 (text.c_str(), -1, NULL, NULL, NULL);
        count = 8 + wcslen(commentStr)*2;
        valuesize = count;
        value = (unsigned char*)new char[valuesize];
        memcpy((char*)value, "UNICODE\0", 8);
        memcpy((char*)value + 8, (char*)commentStr, valuesize - 8);
        g_free(commentStr);
    }
}

void Tag::initString (const char* text)
{

    type = ASCII;
    count = strlen (text) + 1;
    valuesize = count;
    value = new unsigned char[valuesize];
    strcpy ((char*)value, text);
}

void Tag::initSubDir ()
{
    type = LONG;
    valuesize = 4;
    count = 1;
    value = new unsigned char[4];
    setInt (0);
    directory = new TagDirectory*[2];
    directory[0] = new TagDirectory (parent, attrib ? attrib->subdirAttribs : nullptr, parent->getOrder());
    directory[1] = nullptr;
}

void Tag::initSubDir (TagDirectory* dir)
{
    type = LONG;
    valuesize = 4;
    count = 1;
    value = new unsigned char[4];
    setInt (0);
    directory = new TagDirectory*[2];
    directory[0] = dir;
    directory[1] = nullptr;
}

void Tag::initMakerNote (MNKind mnk, const TagAttrib* ta)
{
    type = UNDEFINED;
    valuesize = 4;
    count = 1;
    value = new unsigned char[4];
    setInt (0);
    directory = new TagDirectory*[2];
    directory[0] = new TagDirectory (parent, ta, parent->getOrder());
    directory[1] = nullptr;
    makerNoteKind = mnk;
}

void Tag::initUndefArray (const char* data, int len)
{
    type = UNDEFINED;
    count = valuesize = len;
    value = new unsigned char[valuesize];
    memcpy (value, data, len);
}

void Tag::initLongArray (const char* data, int len)
{
    type = LONG;
    count = (len + 3) / 4;
    valuesize = count * 4;
    value = new unsigned char[valuesize];
    memcpy (value, data, len);
}

void Tag::initRational (int num, int den)
{
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


const TagAttrib* lookupAttrib (const TagAttrib* dir, const char* field)
{

    for (int i = 0; dir[i].ignore != -1; i++)
        if (!strcmp (dir[i].name, field)) {
            return &dir[i];
        }

    return nullptr;
}

void ExifManager::setIFDOffset(unsigned int offset)
{
    IFDOffset = offset;
}

void ExifManager::parseCIFF ()
{

    TagDirectory* root = new TagDirectory (nullptr, ifdAttribs, INTEL);
    Tag* exif = new Tag (root, lookupAttrib (ifdAttribs, "Exif"));
    exif->initSubDir ();
    Tag* mn = new Tag (exif->getDirectory(), lookupAttrib (exifAttribs, "MakerNote"));
    mn->initMakerNote (IFD, canonAttribs);
    root->addTag (exif);
    exif->getDirectory()->addTag (mn);
    parseCIFF (rml->ciffLength, root);
    root->sort ();
}

Tag* ExifManager::saveCIFFMNTag (TagDirectory* root, int len, const char* name)
{
    int s = ftell (f);
    if(s >= 0) {
        char* data = new char [len];
        fread (data, len, 1, f);
        TagDirectory* mn = root->getTag ("Exif")->getDirectory()->getTag ("MakerNote")->getDirectory();
        Tag* cs = new Tag (mn, lookupAttrib (canonAttribs, name));
        cs->initUndefArray (data, len);
        mn->addTag (cs);
        fseek (f, s, SEEK_SET);
        delete [] data;
        return cs;
    } else {
        return nullptr;
    }
}

void ExifManager::parseCIFF (int length, TagDirectory* root)
{

    if (!f) {
        #ifndef NDEBUG
        std::cerr << "ERROR : no file opened !" << std::endl;
        #endif
        return;
    }

    char buffer[1024];
    Tag* t;

    fseek (f, rml->ciffBase + length - 4, SEEK_SET);

    int dirStart = get4 (f, INTEL) + rml->ciffBase;
    fseek (f, dirStart, SEEK_SET);

    int numOfTags = get2 (f, INTEL);

    if (numOfTags > 100) {
        return;
    }

    float exptime, shutter, aperture, fnumber, ev;
    exptime = fnumber = shutter = aperture = ev = -1000.f;
    int focal_len, iso;
    focal_len = iso = -1;

    TagDirectory* exif = root->getTag ("Exif")->getDirectory();

    time_t timestamp = time (nullptr);

    for (int i = 0; i < numOfTags; i++) {

        int type = get2 (f, INTEL);
        int len  = get4 (f, INTEL);
        int nextPos = ftell (f) + 4;

        // seek to the location of the value
        fseek (f, rml->ciffBase + get4 (f, INTEL), SEEK_SET);

        if ((((type >> 8) + 8) | 8) == 0x38) {
            ExifManager(
                f,
                std::unique_ptr<rtengine::RawMetaDataLocation>(
                    new rtengine::RawMetaDataLocation(
                        ftell(f),
                        len
                    )
                ),
                true
            ).parseCIFF(len, root); // Parse a sub-table
        }

        if (type == 0x0810) {
            fread (buffer, 64, 1, f);
            t = new Tag (root, lookupAttrib (ifdAttribs, "Artist"));
            t->initString (buffer);
            root->addTag (t);
        }

        if (type == 0x080a) {
            fread (buffer, 64, 1, f);
            t = new Tag (root, lookupAttrib (ifdAttribs, "Make"));
            t->initString (buffer);
            root->addTag (t);
            fseek (f, strlen (buffer) - 63, SEEK_CUR);
            fread (buffer, 64, 1, f);
            t = new Tag (root, lookupAttrib (ifdAttribs, "Model"));
            t->initString (buffer);
            root->addTag (t);
        }

        if (type == 0x1818) {
            ev = int_to_float (get4 (f, INTEL));
            shutter = int_to_float (get4 (f, INTEL));
            exptime = pow (2, -shutter);
            aperture = int_to_float (get4 (f, INTEL));
            fnumber = pow (2, aperture / 2);

        }

        ExifManager exifManager(f, nullptr, true);
        if (type == 0x102d) {
            Tag* t = exifManager.saveCIFFMNTag (root, len, "CanonCameraSettings");
            int mm = t->toInt (34, SHORT);
            Tag* nt = new Tag (exif, lookupAttrib (exifAttribs, "MeteringMode"));

            switch (mm) {
                case 0:
                    nt->initInt (5, SHORT);
                    break;

                case 1:
                    nt->initInt (3, SHORT);
                    break;

                case 2:
                    nt->initInt (1, SHORT);
                    break;

                case 3:
                    nt->initInt (5, SHORT);
                    break;

                case 4:
                    nt->initInt (6, SHORT);
                    break;

                case 5:
                    nt->initInt (2, SHORT);
                    break;
            }

            exif->addTag (nt);
            nt = new Tag (exif, lookupAttrib (exifAttribs, "MaxApertureValue"));
            nt->initRational (t->toInt (52, SHORT), 32);
            exif->addTag (nt);
            int em = t->toInt (40, SHORT);
            nt = new Tag (exif, lookupAttrib (exifAttribs, "ExposureProgram"));

            switch (em) {
                case 0:
                    nt->initInt (2, SHORT);
                    break;

                case 1:
                    nt->initInt (2, SHORT);
                    break;

                case 2:
                    nt->initInt (4, SHORT);
                    break;

                case 3:
                    nt->initInt (3, SHORT);
                    break;

                case 4:
                    nt->initInt (1, SHORT);
                    break;

                default:
                    nt->initInt (0, SHORT);
                    break;
            }

            exif->addTag (nt);
            nt = new Tag (exif, lookupAttrib (exifAttribs, "Flash"));

            if (t->toInt (8, SHORT) == 0) {
                nt->initInt (0, SHORT);
            } else {
                nt->initInt (1, SHORT);
            }

            exif->addTag (nt);
            nt = new Tag (exif, lookupAttrib (exifAttribs, "MaxApertureValue"));
            nt->initRational (t->toInt (52, SHORT), 32);
            exif->addTag (nt);
        }

        if (type == 0x1029) {
            exifManager.saveCIFFMNTag (root, len, "CanonFocalLength");
        }

        if (type == 0x1031) {
            exifManager.saveCIFFMNTag (root, len, "SensorInfo");
        }

        if (type == 0x1033) {
            exifManager.saveCIFFMNTag (root, len, "CustomFunctions");
        }

        if (type == 0x1038) {
            exifManager.saveCIFFMNTag (root, len, "CanonAFInfo");
        }

        if (type == 0x1093) {
            exifManager.saveCIFFMNTag (root, len, "CanonFileInfo");
        }

        if (type == 0x10a9) {
            exifManager.saveCIFFMNTag (root, len, "ColorBalance");
        }

        if (type == 0x102a) {
            exifManager.saveCIFFMNTag (root, len, "CanonShotInfo");

            iso = pow (2, (get4 (f, INTEL), get2 (f, INTEL)) / 32.0 - 4) * 50;
            aperture  = (get2 (f, INTEL), (short)get2 (f, INTEL)) / 32.0f;
            fnumber = pow (2, aperture / 2);
            shutter = ((short)get2 (f, INTEL)) / 32.0f;
            ev = ((short)get2 (f, INTEL)) / 32.0f;
            fseek (f, 34, SEEK_CUR);

            if (shutter > 1e6) {
                shutter = get2 (f, INTEL) / 10.0f;
            }

            exptime   = pow (2, -shutter);
        }

        if (type == 0x5029) {
            focal_len = len >> 16;

            if ((len & 0xffff) == 2) {
                focal_len /= 32;
            }
        }

//    if (type == 0x5813) flash_used = int_to_float(len);
        if (type == 0x580e) {
            timestamp  = len;
        }

        if (type == 0x180e) {
            timestamp  = get4 (f, INTEL);
        }

        if ((type | 0x4000) == 0x580e) {
            timestamp = mktime (gmtime (&timestamp));
        }

        fseek (f, nextPos, SEEK_SET);
    }

    if (shutter > -999) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "ShutterSpeedValue"));
        t->initRational ((int) (shutter * 10000), 10000);
        exif->addTag (t);
    }

    if (exptime > -999) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "ExposureTime"));
        t->initRational ((int) (exptime * 10000), 10000);
        exif->addTag (t);
    }

    if (aperture > -999) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "ApertureValue"));
        t->initRational ((int) (aperture * 10), 10);
        exif->addTag (t);
    }

    if (fnumber > -999) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "FNumber"));
        t->initRational ((int) (fnumber * 10), 10);
        exif->addTag (t);
    }

    if (ev > -999) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "ExposureBiasValue"));
        t->initRational ((int) (ev * 1000), 1000);
        exif->addTag (t);
    }

    if (iso > 0) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "ISOSpeedRatings"));
        t->initInt (iso, LONG);
        exif->addTag (t);
    }

    if (focal_len > 0) {
        t = new Tag (exif, lookupAttrib (exifAttribs, "FocalLength"));
        t->initRational (focal_len * 32, 32);
        exif->addTag (t);
    }

    if (timestamp != time (nullptr)) {
        struct tm* tim = localtime (&timestamp);
        strftime (buffer, 20, "%Y:%m:%d %H:%M:%S", tim);
        t = new Tag (exif, lookupAttrib (exifAttribs, "DateTimeOriginal"));
        t->initString (buffer);
        exif->addTag (t);
        t = new Tag (exif, lookupAttrib (exifAttribs, "DateTimeDigitized"));
        t->initString (buffer);
        exif->addTag (t);
        t = new Tag (root, lookupAttrib (ifdAttribs, "DateTime"));
        t->initString (buffer);
        root->addTag (t);
    }

    roots.push_back(root);

}

static void
parse_leafdata (TagDirectory* root, ByteOrder order)
{

    Tag *leafdata = root->getTag ("LeafData");

    if (!leafdata) {
        return;
    }

    unsigned char *value = leafdata->getValue();
    int valuesize = leafdata->getValueSize();

    // parse LeafData tag, a tag specific to Leaf digital backs, and has a custom
    // format with 52 byte tag headers starting with "PKTS"
    const char *PKTS_tag = (order == MOTOROLA) ? "PKTS" : "STKP";
    char *hdr;
    int pos = 0;

    // There are lots of sub-tags in here, but for now we only care about those directly
    // useful to RT, which is ISO and rotation. Shutter speed and aperture is not
    // available here.
    int iso_speed = 0;
    int rotation_angle = 0;
    int found_count = 0;

    while (pos + (int)sizeof (hdr) <= valuesize && found_count < 2) {
        hdr = (char *)&value[pos];

        if (strncmp (hdr, PKTS_tag, 4) != 0) {
            // in a few cases the header can be offset a few bytes, don't know why
            // it does not seem to be some sort of alignment, it appears random,
            // this check takes care of it, restart if we find an offset match.
            int offset = 1;

            for (; offset <= 3; offset++) {
                if (strncmp (&hdr[offset], PKTS_tag, 4) == 0) {
                    pos += offset;
                    break;
                }
            }

            if (offset <= 3) {
                continue;
            }

            break;
        }

        int size = sget4 ((unsigned char *)&hdr[48], order);

        if (pos + size > valuesize) {
            break;
        }

        pos += 52;
        char *val = (char *)&value[pos];

        if (strncmp (&hdr[8], "CameraObj_ISO_speed", 19) == 0) {
            iso_speed = 25 * (1 << (atoi (val) - 1));
            found_count++;
        } else if (strncmp (&hdr[8], "ImgProf_rotation_angle", 22) == 0) {
            rotation_angle = atoi (val);
            found_count++;
        } else {
            // check if this is a sub-directory, include test for that strange offset of next header
            if (size >= 8 &&
                    (strncmp (val, PKTS_tag, 4) == 0 ||
                     strncmp (&val[1], PKTS_tag, 4) == 0 ||
                     strncmp (&val[2], PKTS_tag, 4) == 0 ||
                     strncmp (&val[3], PKTS_tag, 4) == 0)) {
                // start of next hdr, this is a sub-directory, we skip those for now.
                size = 0;
            }
        }

        pos += size;
    }

    // create standard tags from the custom Leaf tags
    Tag* exif = root->getTag ("Exif");

    if (!exif) {
        exif = new Tag (root, root->getAttrib ("Exif"));
        exif->initSubDir();
        root->addTagFront (exif);
    }

    if (!exif->getDirectory()->getTag ("ISOSpeedRatings")) {
        Tag *t = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
        t->initInt (iso_speed, LONG);
        exif->getDirectory()->addTagFront (t);
    }

    if (!root->getTag ("Orientation")) {
        int orientation;

        switch (rotation_angle) {
            case 0:
                orientation = 1;
                break;

            case 90:
                orientation = 6;
                break;

            case 180:
                orientation = 3;
                break;

            case 270:
                orientation = 8;
                break;

            default:
                orientation = 1;
                break;
        }

        Tag *t = new Tag (root, root->getAttrib ("Orientation"));
        t->initInt (orientation, SHORT);
        root->addTagFront (t);
    }

    // now look in ApplicationNotes tag for additional information
    Tag *appnotes = root->getTag ("ApplicationNotes");

    if (!appnotes) {
        return;
    }

    char *xmp = (char *)appnotes->getValue();
    char *end, *p;

    // Quick-and-dirty value extractor, no real xml parsing.
    // We could make it more generic, but we just get most important
    // values we know use to be in there.
    if ((p = strstr (xmp, "xmlns:tiff")) != nullptr &&
            (end = strstr (p, "</rdf:Description>")) != nullptr) {
        *end = '\0';

        while ((p = strstr (p, "<tiff:")) != nullptr) {
            char *tag = &p[6], *tagend;

            if ((tagend = strchr (tag, '>')) == nullptr) {
                break;
            }

            *tagend = '\0';
            char *val = &tagend[1];

            if ((p = strstr (val, "</tiff:")) == nullptr) {
                *tagend = '>';
                break;
            }

            *p = '\0';

            if (root->getAttrib (tag) && !root->getTag (tag)) {
                Tag *t = new Tag (root, root->getAttrib (tag));

                if (strcmp (tag, "Make") == 0 ||
                        strcmp (tag, "Model") == 0) {
                    if (strcmp (tag, "Model") == 0) {
                        // Leaf adds back serial number and camera model to the 'Model'
                        // tag, we strip that away here so the back can be recognized
                        // and matched against DCP profile
                        char *p1 = strchr (val, '(');

                        if (p1 != nullptr) {
                            *p1 = '\0';
                        }

                        // Model name also contains a leading "Leaf " which we already
                        // have in the Make name, remove that.
                        if (strstr (val, "Leaf ") == val) {
                            t->initString (&val[5]);
                        } else {
                            t->initString (val);
                        }

                        if (p1 != nullptr) {
                            *p1 = '(';
                        }
                    } else {
                        t->initString (val);
                    }

                    root->addTagFront (t);
                } else {
                    delete t;
                }
            }

            *p = '<';
            *tagend = '>';
        }

        *end = '<';
    }

    if ((p = strstr (xmp, "xmlns:exif")) != nullptr &&
            (end = strstr (p, "</rdf:Description>")) != nullptr) {
        *end = '\0';

        while ((p = strstr (p, "<exif:")) != nullptr) {
            char *tag = &p[6], *tagend;

            if ((tagend = strchr (tag, '>')) == nullptr) {
                break;
            }

            *tagend = '\0';
            char *val = &tagend[1];

            if ((p = strstr (val, "</exif:")) == nullptr) {
                *tagend = '>';
                break;
            }

            *p = '\0';

            if (exif->getDirectory()->getAttrib (tag) && !exif->getDirectory()->getTag (tag)) {
                Tag *t = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib (tag));
                int num, denom;
                struct tm tm;

                if (strcmp (tag, "ApertureValue") == 0 && sscanf (val, "%d/%d", &num, &denom) == 2) {
                    t->initRational (num, denom);
                    exif->getDirectory()->addTagFront (t);
                    // we also make an "FNumber" tag since many tools don't interpret ApertureValue
                    // according to Exif standard
                    t = new Tag (exif->getDirectory(), lookupAttrib (exifAttribs, "FNumber"));
                    double f = pow (sqrt (2.0), ((double)num / denom));

                    if (f > 10.0) {
                        t->initRational ((int)floor (f), 1);
                    } else {
                        t->initRational ((int)floor (f * 10.0), 10);
                    }

                    exif->getDirectory()->addTagFront (t);
                } else if (strcmp (tag, "ShutterSpeedValue") == 0 && sscanf (val, "%d/%d", &num, &denom) == 2) {
                    t->initRational (num, denom);
                    exif->getDirectory()->addTagFront (t);
                    // we also make an "ExposureTime" tag since many tools don't interpret ShutterSpeedValue
                    // according to Exif standard
                    t = new Tag (exif->getDirectory(), lookupAttrib (exifAttribs, "ExposureTime"));
                    double f = 1.0 / pow (2.0, ((double)num / denom));

                    if (f > 10.0) {
                        t->initRational ((int)floor (f), 1);
                    } else if (f > 1.0) {
                        t->initRational ((int)floor (f * 10.0), 10);
                    } else if (f == 1.0) {
                        t->initRational (1, 1);
                    } else {
                        f = 1.0 / f;
                        static const double etimes[] = {
                            10000, 8000, 6400, 6000, 5000,
                            4000, 3200, 3000, 2500,
                            2000, 1600, 1500, 1250,
                            1000, 800, 750, 640,
                            500, 400, 350, 320,
                            250, 200, 180, 160,
                            125, 100, 90, 80,
                            60, 50, 45, 40,
                            30, 25, 22, 20,
                            15, 13, 11, 10,
                            8, 6, 5,
                            4, 3, 2.5,
                            2, 1.6, 1.5, 1.3,
                            1, -1
                        };
                        double diff = etimes[0];
                        int idx = -1;

                        for (int i = 1; etimes[i] > 0; i++) {
                            if (abs (etimes[i] - f) < diff) {
                                idx = i;
                                diff = abs (etimes[i] - f);
                            }
                        }

                        if (idx != -1 && f < etimes[0]) {
                            f = etimes[idx];
                        }

                        if (f < 2) {
                            t->initRational (10, (int) (10 * f));
                        } else {
                            t->initRational (1, (int)f);
                        }
                    }

                    exif->getDirectory()->addTagFront (t);
                } else if (strcmp (tag, "FocalLength") == 0 && sscanf (val, "%d/%d", &num, &denom) == 2) {
                    t->initRational (num, denom);
                    exif->getDirectory()->addTagFront (t);
                } else if (strcmp (tag, "ISOSpeedRatings") == 0) {
                    char *p1 = val;

                    while (*p1 != '\0' && !isdigit (*p1)) {
                        p1++;
                    }

                    if (*p1 != '\0') {
                        t->initInt (atoi (p1), LONG);
                        exif->getDirectory()->addTagFront (t);
                    } else {
                        delete t;
                    }
                } else if (strcmp (tag, "DateTimeOriginal") == 0 &&
                           sscanf (val, "%d-%d-%dT%d:%d:%dZ",
                                   &tm.tm_year, &tm.tm_mon,
                                   &tm.tm_mday, &tm.tm_hour,
                                   &tm.tm_min, &tm.tm_sec) == 6) {
                    char tstr[64];
                    sprintf (tstr, "%04d:%02d:%02d %02d:%02d:%02d", tm.tm_year, tm.tm_mon,
                             tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
                    t->initString (tstr);
                    exif->getDirectory()->addTagFront (t);
                } else {
                    delete t;
                }
            }

            *p = '<';
            *tagend = '>';
        }

        *end = '<';
    }
}

void ExifManager::parseRaw (bool skipIgnored) {
    parse(true, skipIgnored);
}

void ExifManager::parseStd (bool skipIgnored) {
    parse(false, skipIgnored);
}

void ExifManager::parse (bool isRaw, bool skipIgnored)
{
    int ifdOffset = IFDOffset;

    if (!f) {
        #ifndef NDEBUG
        std::cerr << "ERROR : no file opened !" << std::endl;
        #endif
        return;
    }
    setlocale (LC_NUMERIC, "C"); // to set decimal point in sscanf

    if (order == ByteOrder::UNKNOWN) {
        // read tiff header
        fseek (f, rml->exifBase, SEEK_SET);
        unsigned short bo;
        fread (&bo, 1, 2, f);
        order = (ByteOrder) ((int)bo);
        get2 (f, order);
        if (!ifdOffset) {
            ifdOffset = get4 (f, order);
        }
    }

    do {
        // seek to IFD
        fseek (f, rml->exifBase + ifdOffset, SEEK_SET);

        // first read the IFD directory
        TagDirectory* root =  new TagDirectory (nullptr, f, rml->exifBase, ifdAttribs, order, skipIgnored);

        // fix ISO issue with nikon and panasonic cameras
        Tag* make = root->getTag ("Make");
        Tag* exif = root->getTag ("Exif");

        if (exif && !exif->getDirectory()->getTag ("ISOSpeedRatings")) {
            if (make && !strncmp ((char*)make->getValue(), "NIKON", 5)) {
                Tag* mn   = exif->getDirectory()->getTag ("MakerNote");

                if (mn) {
                    Tag* iso = mn->getDirectory()->getTag ("ISOSpeed");

                    if (iso) {
                        std::string isov = iso->valueToString ();
                        Tag* niso = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
                        niso->initInt (atoi (isov.c_str()), SHORT);
                        exif->getDirectory()->addTagFront (niso);
                    }
                }
            } else if (make && (!strncmp ((char*)make->getValue(), "Panasonic", 9) || !strncmp ((char*)make->getValue(), "LEICA", 5))) {
                Tag* iso = root->getTag ("PanaISO");

                if (iso) {
                    std::string isov = iso->valueToString ();
                    Tag* niso = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
                    niso->initInt (atoi (isov.c_str()), SHORT);
                    exif->getDirectory()->addTagFront (niso);
                }
            }
        }

        if (make && !strncmp ((char*)make->getValue(), "Kodak", 5)) {
            if (!exif) {
                // old Kodak cameras may have exif tags in IFD0, reparse and create an exif subdir
                fseek (f, rml->exifBase + ifdOffset, SEEK_SET);
                TagDirectory* exifdir =  new TagDirectory (nullptr, f, rml->exifBase, exifAttribs, order, true);

                exif = new Tag (root, root->getAttrib ("Exif"));
                exif->initSubDir (exifdir);
                root->addTagFront (exif);

                if (!exif->getDirectory()->getTag ("ISOSpeedRatings") && exif->getDirectory()->getTag ("ExposureIndex")) {
                    Tag* niso = new Tag (exif->getDirectory(), exif->getDirectory()->getAttrib ("ISOSpeedRatings"));
                    niso->initInt (exif->getDirectory()->getTag ("ExposureIndex")->toInt(), SHORT);
                    exif->getDirectory()->addTagFront (niso);
                }
            }

            Tag *kodakIFD = root->getTag ("KodakIFD");

            if (kodakIFD && kodakIFD->getDirectory()->getTag ("TextualInfo")) {
                parseKodakIfdTextualInfo (kodakIFD->getDirectory()->getTag ("TextualInfo"), exif);
            }
        }

        parse_leafdata (root, order);

        if (make && !strncmp ((char*)make->getValue(), "Hasselblad", 10)) {
            /*
            Figuring out the Hasselblad model is a mess. Hasselblad raw data comes in four slightly
            different containers, 3FR (directly from CF card), FFF (same as 3FR but filtered through
            Phocus, calibration data applied and a bit different tags), Adobe-generated DNGs and
            Phocus-generated DNGs.

            FFF usually has a sane model name in Model (and is used as reference for what we shall
            call the different Hasselblad models), but 3FR only says like "Hasselblad H3D" for
            all H3D models, or "Flash Sync" if the back has been used on a mechanical camera body.
            V-mount backs may have the model name of the V body instead of the back model. Etc...
            as said it's a mess.

            This code is supposed to handle all raw containers and end up with the same model
            regardless of container.

            We don't differ between single shot and multi-shot models, and probably there's no use
            of doing so. You need Hasselblad's own software to shoot multi-shot and can only do that
            tethered. In single-shot mode they should be exactly the same as the single-shot models.
                  */
            Tag *subd = root->getTag (0x14a);
            Tag *iw = (subd) ? subd->getDirectory()->getTag ("ImageWidth") : nullptr;
            int sensorWidth = (iw) ? iw->toInt() : 0;
            Tag* tmodel = root->getTag ("Model");
            const char *model = (tmodel) ? (const char *)tmodel->getValue() : "";

            if (strstr (model, "Hasselblad ") == model) {
                model += 11;
            } else {
                // if HxD is used in flash sync mode for example, we need to fetch model from this tag
                Tag* tmodel3 = root->getTag ("UniqueCameraModel");
                const char *model3 = (tmodel3) ? (const char *)tmodel3->getValue() : "";

                if (strstr (model3, "Hasselblad ") == model3) {
                    model = model3 + 11;
                }
            }

            // FIXME: due to lack of test files this Hasselblad model identification is not 100% complete
            // This needs checking out: CFV-39/CFV-50 3FR, H3DII vs H3D, old CF/CFH models

            if (!strcmp (model, "H3D")) {
                // We can't differ between H3D and H3DII for the 22, 31 and 39 models. There's was no H3D-50 so we know that is a
                // H3DII-50. At the time of writing I have no test files for the H3D vs H3DII models, so there still may be a chance
                // to differ between them. AFAIK Adobe's DNG converter don't differ between them, and actually call the H3DII-50
                // H3D-50 although Hasselblad never released such a model.
                switch (sensorWidth) {
                    case 4096:
                        tmodel->initString ("H3D-22");
                        break;

                    case 6542:
                        tmodel->initString ("H3D-31");
                        break;

                    case 7262:
                        tmodel->initString ("H3D-39");
                        break;

                    case 8282:
                        tmodel->initString ("H3DII-50");
                        break;
                }
            } else if (!strcmp (model, "H4D")) {
                switch (sensorWidth) {
                    case 6542:
                        tmodel->initString ("H4D-31");
                        break;

                    case 7410:
                        tmodel->initString ("H4D-40");
                        break;

                    case 8282:
                        tmodel->initString ("H4D-50");
                        break;

                    case 9044:
                        tmodel->initString ("H4D-60");
                        break;
                }
            } else if (!strcmp (model, "H5D")) {
                switch (sensorWidth) {
                    case 7410:
                        tmodel->initString ("H5D-40");
                        break;

                    case 8282:
                        tmodel->initString ("H5D-50");
                        break;

                    case 8374:
                        tmodel->initString ("H5D-50c");
                        break;

                    case 9044:
                        tmodel->initString ("H5D-60");
                        break;
                }
            } else if (!strcmp (model, "CFV")) {
                switch (sensorWidth) {
                    case 7262:
                        tmodel->initString ("CFV-39");
                        break;

                    case 8282:
                        tmodel->initString ("CFV-50");
                        break;

                    case 8374:
                        tmodel->initString ("CFV-50c");
                        break;
                }
            }

            // and a few special cases
            Tag* tmodel3 = root->getTag ("UniqueCameraModel");
            const char *model3 = (tmodel3) ? (const char *)tmodel3->getValue() : "";

            if (strstr (model3, "Hasselblad ") == model3) {
                model3 = model3 + 11;
            }

            if (!strcmp (model3, "ixpressCF132")) {
                tmodel->initString ("CF-22");
            } else if (!strcmp (model3, "Hasselblad96")) {
                tmodel->initString ("CFV"); // popularly called CFV-16, but the official name is CFV
            } else if (!strcmp (model3, "Hasselblad234")) {
                tmodel->initString ("CFV-39");
            } else if (sensorWidth == 4090) {
                tmodel->initString ("V96C");
            }

            // and yet some, this is for Adobe-generated DNG files
            Tag* tmodel4 = root->getTag ("LocalizedCameraModel");

            if (tmodel4) {
                const char *model4 = (const char *)tmodel4->getValue();

                if (strstr (model4, "Hasselblad ") == model4) {
                    model4 = model4 + 11;
                }

                if (!strcmp (model4, "ixpressCF132-22")) {
                    tmodel->initString ("CF-22");
                } else if (!strcmp (model4, "Hasselblad96-16")) {
                    tmodel->initString ("CFV");
                } else if (!strcmp (model4, "Hasselblad234-39")) {
                    tmodel->initString ("CFV-39");
                } else if (!strcmp (model4, "H3D-50")) {
                    // Adobe names H3DII-50 incorrectly as H3D-50
                    tmodel->initString ("H3DII-50");
                } else if (strstr (model4, "H3D-") == model4 || strstr (model4, "H4D-") == model4 || strstr (model4, "H5D-") == model4) {
                    tmodel->initString (model4);
                }
            }
        }

        if (!root->getTag ("Orientation")) {
            if (make && !strncmp ((char*)make->getValue(), "Phase One", 9)) {
                int orientation = 0;
                Tag *iw = root->getTag ("ImageWidth");

                if (iw) {
                    // from dcraw, derive orientation from image width
                    orientation = "0653"[iw->toInt() & 3] - '0';
                }

                Tag *t = new Tag (root, root->getAttrib ("Orientation"));
                t->initInt (orientation, SHORT);
                root->addTagFront (t);
            }
        }

        // --- detecting image root IFD based on SubFileType, or if not provided, on PhotometricInterpretation

        bool frameRootDetected = false;

        if(!frameRootDetected) {
            std::vector<const Tag*> risTagList = root->findTags("RawImageSegmentation");
            if (!risTagList.empty()) {
                for (auto ris : risTagList) {
                    frames.push_back(ris->getParent());
                    frameRootDetected = true;

    #if PRINT_METADATA_TREE
                    printf("\n--------------- FRAME (RAWIMAGESEGMENTATION) ---------------\n\n");
                    ris->getParent()->printAll ();
    #endif
                }
            }
        }

        if(!frameRootDetected) {
            std::vector<const Tag*> sftTagList = root->findTags(TIFFTAG_SUBFILETYPE);
            if (!sftTagList.empty()) {
                for (auto sft : sftTagList) {
                    int sftVal = sft->toInt();
                    if (sftVal == (isRaw ? 0 : 2)) {
                        frames.push_back(sft->getParent());
                        frameRootDetected = true;

#if PRINT_METADATA_TREE
                        printf("\n--------------- FRAME (SUBFILETYPE) ---------------\n\n");
                        sft->getParent()->printAll ();
#endif
                    }
                }
            }
        }

        if(!frameRootDetected) {
            std::vector<const Tag*> sftTagList = root->findTags(TIFFTAG_OSUBFILETYPE);
            if (!sftTagList.empty()) {
                for (auto sft : sftTagList) {
                    int sftVal = sft->toInt();
                    if (sftVal == OFILETYPE_IMAGE) {
                        frames.push_back(sft->getParent());
                        frameRootDetected = true;

#if PRINT_METADATA_TREE
                        printf("\n--------------- FRAME (OSUBFILETYPE) ---------------\n\n");
                        sft->getParent()->printAll ();
#endif
                    }
                }
            }
        }

        if(!frameRootDetected) {
            std::vector<const Tag*> piTagList = root->findTags("PhotometricInterpretation");
            if (!piTagList.empty()) {
                for (auto pi : piTagList) {
                    int piVal = pi->toInt();
                    if (piVal == (isRaw ? 32803 : 2)) {
                        frames.push_back(pi->getParent());
                        //frameRootDetected = true;  not used afterward

#if PRINT_METADATA_TREE
                        printf("\n--------------- FRAME (PHOTOMETRIC) ---------------\n\n");
                        pi->getParent()->printAll ();
#endif
                    }
                }
            }
        }

        // --- getting next sibling root

        ifdOffset = get4 (f, order);

        roots.push_back(root);

#if PRINT_METADATA_TREE
        printf("\n~~~~~~~~~ ROOT ~~~~~~~~~~~~~~~~~~~~~~~~\n\n");
        root->printAll ();
#endif

    } while (ifdOffset && !onlyFirst);

    // Security check : if there's at least one root, there must be at least one image.
    // If the following occurs, then image detection above has failed or it's an unsupported file type.
    // Yet the result of this should be valid.
    if (!roots.empty() && frames.empty()) {
        frames.push_back(roots.at(0));
    }
}

void ExifManager::parseJPEG (int offset)
{
    if (!f) {
        #ifndef NDEBUG
        std::cerr << "ERROR : no file opened !" << std::endl;
        #endif
        return;
    }

    if(!fseek (f, offset, SEEK_SET)) {
        unsigned char c;
        if(fread (&c, 1, 1, f) == 1) {
            constexpr unsigned char markerl = 0xff;
            const char exifid[] = "Exif\0\0";
            char idbuff[8];
            int tiffbase = -1;

            while (fread (&c, 1, 1, f)) {
                if (c != markerl) {
                    continue;
                }

                if (fread (&c, 1, 1, f) && c == 0xe1) { // APP1 marker found
                    if (fread (idbuff, 1, 8, f) < 8) {
                        return;
                    }

                    if (!memcmp (idbuff + 2, exifid, 6)) {  // Exif info found
                        tiffbase = ftell (f);

                        // We need a RawMetaDataLocation to put the 'tiffbase' value
                        const bool rmlCreated = !rml;
                        if (rmlCreated) {
                            rml.reset(new rtengine::RawMetaDataLocation(0));
                        }
                        rml->exifBase = tiffbase;
                        parse (false);
                        if (rmlCreated) {
                            rml.reset();
                        }
                        return;
                    }
                }
            }
        }
    }
}

void ExifManager::parseTIFF (bool skipIgnored)
{
    if (!rml) {
        rml.reset(new rtengine::RawMetaDataLocation(0));
        parse(false, skipIgnored);
        rml.reset();
    } else {
        parse (false,skipIgnored);
    }
}

std::vector<Tag*> ExifManager::getDefaultTIFFTags (TagDirectory* forthis)
{

    std::vector<Tag*> defTags;

    defTags.reserve (12);
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "ImageWidth"), 0, LONG));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "ImageHeight"), 0, LONG));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "XResolution"), 300, RATIONAL));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "YResolution"), 300, RATIONAL));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "ResolutionUnit"), 2, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "Software"), "RawTherapee " RTVERSION));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "Orientation"), 1, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "SamplesPerPixel"), 3, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "BitsPerSample"), 8, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "PlanarConfiguration"), 1, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "PhotometricInterpretation"), 2, SHORT));
    defTags.push_back (new Tag (forthis, lookupAttrib (ifdAttribs, "Compression"), 1, SHORT));

    return defTags;
}



int ExifManager::createJPEGMarker (const TagDirectory* root, const rtengine::procparams::ExifPairs& changeList, int W, int H, unsigned char* buffer)
{

    // write tiff header
    int offs = 6;
    memcpy (buffer, "Exif\0\0", 6);
    ByteOrder order = INTEL;

    if (root) {
        order = root->getOrder ();
    }

    sset2 ((unsigned short)order, buffer + offs, order);
    offs += 2;
    sset2 (42, buffer + offs, order);
    offs += 2;
    sset4 (8, buffer + offs, order);

    TagDirectory* cl;

    if (root) {
        cl = (const_cast<TagDirectory*> (root))->clone (nullptr);
    } else {
        cl = new TagDirectory (nullptr, ifdAttribs, INTEL);
    }

    for (rtengine::procparams::ExifPairs::const_iterator i = changeList.begin(); i != changeList.end(); ++i) {
        cl->applyChange (i->first, i->second);
    }

    const std::vector<Tag*> defTags = getDefaultTIFFTags (cl);

    defTags[0]->setInt (W, 0, LONG);
    defTags[1]->setInt (H, 0, LONG);
    defTags[8]->setInt (8, 0, SHORT);

    for (int i = defTags.size() - 1; i >= 0; i--) {
        Tag* defTag = defTags[i];
        cl->replaceTag (defTag->clone (cl));
        delete defTag;
    }

    cl->sort ();
    int size = cl->write (8, buffer + 6);

    delete cl;

    return size + 6;
}

int ExifManager::createTIFFHeader (const TagDirectory* root, const rtengine::procparams::ExifPairs& changeList, int W, int H, int bps, const char* profiledata, int profilelen, const char* iptcdata, int iptclen, unsigned char *&buffer, unsigned &bufferSize)
{

// write tiff header
    int offs = 0;
    ByteOrder order = HOSTORDER;

    if (root) {
        order = root->getOrder ();
    }

    TagDirectory* cl;

    if (root) {
        cl = (const_cast<TagDirectory*> (root))->clone (nullptr);
        // remove some unknown top level tags which produce warnings when opening a tiff
        Tag *removeTag = cl->getTag (0x9003);

        if (removeTag) {
            removeTag->setKeep (false);
        }

        removeTag = cl->getTag (0x9211);

        if (removeTag) {
            removeTag->setKeep (false);
        }
    } else {
        cl = new TagDirectory (nullptr, ifdAttribs, HOSTORDER);
    }

// add tiff strip data
    int rps = 8;
    int strips = ceil ((double)H / rps);
    cl->replaceTag (new Tag (cl, lookupAttrib (ifdAttribs, "RowsPerStrip"), rps, LONG));
    Tag* stripBC   = new Tag (cl, lookupAttrib (ifdAttribs, "StripByteCounts"));
    stripBC->initInt (0, LONG, strips);
    cl->replaceTag (stripBC);
    Tag* stripOffs = new Tag (cl, lookupAttrib (ifdAttribs, "StripOffsets"));
    stripOffs->initInt (0, LONG, strips);
    cl->replaceTag (stripOffs);

    for (int i = 0; i < strips - 1; i++) {
        stripBC->setInt (rps * W * 3 * bps / 8, i * 4);
    }

    int remaining = (H - rps * floor ((double)H / rps)) * W * 3 * bps / 8;

    if (remaining) {
        stripBC->setInt (remaining, (strips - 1) * 4);
    } else {
        stripBC->setInt (rps * W * 3 * bps / 8, (strips - 1) * 4);
    }

    if (profiledata) {
        Tag* icc = new Tag (cl, lookupAttrib (ifdAttribs, "ICCProfile"));
        icc->initUndefArray (profiledata, profilelen);
        cl->replaceTag (icc);
    }

    if (iptcdata) {
        Tag* iptc = new Tag (cl, lookupAttrib (ifdAttribs, "IPTCData"));
        iptc->initLongArray (iptcdata, iptclen);
        cl->replaceTag (iptc);
    }

// apply list of changes
    for (rtengine::procparams::ExifPairs::const_iterator i = changeList.begin(); i != changeList.end(); ++i) {
        cl->applyChange (i->first, i->second);
    }

    // append default properties
    const std::vector<Tag*> defTags = getDefaultTIFFTags (cl);

    defTags[0]->setInt (W, 0, LONG);
    defTags[1]->setInt (H, 0, LONG);
    defTags[8]->initInt (0, SHORT, 3);

    for (int i = 0; i < 3; i++) {
        defTags[8]->setInt (bps, i * 2, SHORT);
    }

    for (int i = defTags.size() - 1; i >= 0; i--) {
        Tag* defTag = defTags[i];
        cl->replaceTag (defTag->clone (cl));
        delete defTag;
    }

// calculate strip offsets
    int size = cl->calculateSize ();
    int byps = bps / 8;

    for (int i = 0; i < strips; i++) {
        stripOffs->setInt (size + 8 + i * rps * W * 3 * byps, i * 4);
    }

    cl->sort ();
    bufferSize = cl->calculateSize() + 8;
    buffer = new unsigned char[bufferSize]; // this has to be deleted in caller
    sset2 ((unsigned short)order, buffer + offs, order);
    offs += 2;
    sset2 (42, buffer + offs, order);
    offs += 2;
    sset4 (8, buffer + offs, order);

    int endOffs = cl->write (8, buffer);

//  cl->printAll();
    delete cl;

    return endOffs;
}


int ExifManager::createPNGMarker(const TagDirectory* root, const rtengine::procparams::ExifPairs &changeList, int W, int H, int bps, const char* iptcdata, int iptclen, unsigned char *&buffer, unsigned &bufferSize)
{
// write tiff header
    int offs = 0;
    ByteOrder order = HOSTORDER;

    if (root) {
        order = root->getOrder ();
    }

    TagDirectory* cl;

    if (root) {
        cl = (const_cast<TagDirectory*> (root))->clone (nullptr);
        // remove some unknown top level tags which produce warnings when opening a tiff
        Tag *removeTag = cl->getTag (0x9003);

        if (removeTag) {
            removeTag->setKeep (false);
        }

        removeTag = cl->getTag (0x9211);

        if (removeTag) {
            removeTag->setKeep (false);
        }
    } else {
        cl = new TagDirectory (nullptr, ifdAttribs, HOSTORDER);
    }

    if (iptcdata) {
        Tag* iptc = new Tag (cl, lookupAttrib (ifdAttribs, "IPTCData"));
        iptc->initLongArray (iptcdata, iptclen);
        cl->replaceTag (iptc);
    }

// apply list of changes
    for (rtengine::procparams::ExifPairs::const_iterator i = changeList.begin(); i != changeList.end(); ++i) {
        cl->applyChange (i->first, i->second);
    }

    // append default properties
    const std::vector<Tag*> defTags = getDefaultTIFFTags (cl);

    defTags[0]->setInt (W, 0, LONG);
    defTags[1]->setInt (H, 0, LONG);
    defTags[8]->initInt (0, SHORT, 3);

    for (int i = 0; i < 3; i++) {
        defTags[8]->setInt (bps, i * 2, SHORT);
    }

    for (int i = defTags.size() - 1; i >= 0; i--) {
        Tag* defTag = defTags[i];
        cl->replaceTag (defTag->clone (cl));
        delete defTag;
    }

    cl->sort ();
    bufferSize = cl->calculateSize() + 8;
    buffer = new unsigned char[bufferSize]; // this has to be deleted in caller
    sset2 ((unsigned short)order, buffer + offs, order);
    offs += 2;
    sset2 (42, buffer + offs, order);
    offs += 2;
    sset4 (8, buffer + offs, order);

    int endOffs = cl->write (8, buffer);

//  cl->printAll();
    delete cl;

    return endOffs;
}


//-----------------------------------------------------------------------------
// global functions to read byteorder dependent data
//-----------------------------------------------------------------------------
unsigned short sget2 (unsigned char *s, rtexif::ByteOrder order)
{

    if (order == rtexif::INTEL) {
        return s[0] | s[1] << 8;
    } else {
        return s[0] << 8 | s[1];
    }
}

int sget4 (unsigned char *s, rtexif::ByteOrder order)
{

    if (order == rtexif::INTEL) {
        return s[0] | s[1] << 8 | s[2] << 16 | s[3] << 24;
    } else {
        return s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3];
    }
}

inline unsigned short get2 (FILE* f, rtexif::ByteOrder order)
{

    unsigned char str[2] = { 0xff, 0xff };
    fread (str, 1, 2, f);
    return rtexif::sget2 (str, order);
}

int get4 (FILE* f, rtexif::ByteOrder order)
{

    unsigned char str[4] = { 0xff, 0xff, 0xff, 0xff };
    fread (str, 1, 4, f);
    return rtexif::sget4 (str, order);
}

void sset2 (unsigned short v, unsigned char *s, rtexif::ByteOrder order)
{

    if (order == rtexif::INTEL) {
        s[0] = v & 0xff;
        v >>= 8;
        s[1] = v;
    } else {
        s[1] = v & 0xff;
        v >>= 8;
        s[0] = v;
    }
}

void sset4 (int v, unsigned char *s, rtexif::ByteOrder order)
{

    if (order == rtexif::INTEL) {
        s[0] = v & 0xff;
        v >>= 8;
        s[1] = v & 0xff;
        v >>= 8;
        s[2] = v & 0xff;
        v >>= 8;
        s[3] = v;
    } else {
        s[3] = v & 0xff;
        v >>= 8;
        s[2] = v & 0xff;
        v >>= 8;
        s[1] = v & 0xff;
        v >>= 8;
        s[0] = v;
    }
}

float int_to_float (int i)
{
    union {
        int i;
        float f;
    } u;
    u.i = i;
    return u.f;
}

short int int2_to_signed (short unsigned int i)
{
    union {
        short unsigned int i;
        short int s;
    } u;
    u.i = i;
    return u.s;
}

/* Function to parse and extract focal length and aperture information from description
 * @fullname must conform to the following formats
 * <focal>mm f/<aperture>
 * <focal>-<focal>mm f/<aperture>
 * <focal>-<focal>mm f/<aperture>-<aperture>
 * NB: no space between separator '-'; no space between focal length and 'mm'
 */
bool extractLensInfo (std::string &fullname, double &minFocal, double &maxFocal, double &maxApertureAtMinFocal, double &maxApertureAtMaxFocal)
{
    minFocal = 0.0;
    maxFocal = 0.0;
    maxApertureAtMinFocal = 0.0;
    maxApertureAtMaxFocal = 0.0;
    char buffer[1025];
    strncpy (buffer, fullname.c_str(), 1024);
    char *pF = strstr (buffer, "f/" );

    if ( pF ) {
        sscanf (pF + 2, "%lf-%lf", &maxApertureAtMinFocal, &maxApertureAtMaxFocal);

        if (maxApertureAtMinFocal > 0. && maxApertureAtMaxFocal == 0.) {
            maxApertureAtMaxFocal = maxApertureAtMinFocal;
        }

        char *pMM = pF - 3;

        while ( pMM[0] != 'm' && pMM[1] != 'm' && pMM > buffer) {
            pMM--;
        }

        if ( pMM[0] == 'm' && pMM[1] == 'm' ) {
            char *sp = pMM;

            while ( *sp != ' ' && sp > buffer ) {
                sp--;
            }

            sscanf (sp + 1, "%lf-%lf", &minFocal, &maxFocal);

            if (maxFocal == 0.) {
                maxFocal = minFocal;
            }

            return true;
        }
    }

    return false;
}

}

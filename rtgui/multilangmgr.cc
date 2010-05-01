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
#include <glib/gstdio.h>
#include <multilangmgr.h>
#include <string.h>
#include <safegtk.h>

MultiLangMgr langMgr;

Glib::ustring M (std::string key) { return langMgr.getStr (key); }

bool MultiLangMgr::load (Glib::ustring fname, MultiLangMgr* fb) {

    fallBack = fb;

    FILE *f = g_fopen (fname.c_str(), "rt");

    if (f==NULL)
        return false;
    
    transTable.clear ();
    
    char* buffer = new char[2048];
    
    while (buffer = fgets (buffer, 2048, f)) {
        // find separator
        int seppos = 0;
        while (buffer[seppos]!=0 && buffer[seppos]!=';')
            seppos++;
        // no separator found
        if (buffer[seppos]==0) 
            continue;
        // cut the last \n and \r characters
        int endpos = strlen(buffer)-1;
        while (buffer[endpos]=='\n' || buffer[endpos]=='\r')
            endpos--;
        buffer[endpos+1] = 0;
        // replace "\n" to '\n'
        int j = 0;
        for (int i=0; i<endpos+1; i++) 
            if (i<endpos && buffer[i]=='\\' && buffer[i+1]=='n') {
                buffer[j++] = '\n';
                i++;
            }
            else
                buffer[j++] = buffer[i];
        buffer[j] = 0;
        // cut to two parts
        buffer[seppos] = 0;
        transTable[buffer] = buffer + seppos + 1;
    }
    
    fclose (f);
    delete [] buffer;
    return true;
}

bool MultiLangMgr::save (Glib::ustring fname) {

    FILE *f = g_fopen (fname.c_str(), "wt");
    
    if (f==NULL)
        return false;

    std::map<std::string, Glib::ustring>::iterator r;
    for (r=transTable.begin (); r!=transTable.end(); r++)
        fprintf (f, "%s;%s\n", r->first.c_str(), safe_locale_from_utf8(r->second).c_str());

    fclose (f);
    return true;
}
        
Glib::ustring MultiLangMgr::getStr (std::string key) {

    std::map<std::string, Glib::ustring>::iterator r = transTable.find (key);
    if (r!=transTable.end()) 
        return r->second;
    else if (fallBack)
        return fallBack->getStr (key);
    else
        return "";
}

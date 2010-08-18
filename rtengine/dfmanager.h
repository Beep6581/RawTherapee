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
#include <string>
#include <glibmm/ustring.h>
#include <map>
#include <math.h>
#include <common.h>

namespace rtengine{

class dfInfo
{
public:

	Glib::ustring pathname; // filename of dark frame
	std::list<Glib::ustring> pathNames; // other similar dark frames, used to mediate
	std::string maker;  ///< manufacturer
	std::string model;  ///< model
	int iso;      ///< ISO (gain)
	double shutter;     ///< shutter or exposure time in sec
	time_t timestamp;   ///< seconds since 1 Jan 1970


	dfInfo(const Glib::ustring &name, const std::string &mak, const std::string &mod,int iso,double shut,time_t t)
	:pathname(name),maker(mak),model(mod),iso(iso),shutter(shut),timestamp(t),ri(NULL){}

	dfInfo( const dfInfo &o)
	:pathname(o.pathname),maker(o.maker),model(o.model),iso(o.iso),shutter(o.shutter),timestamp(o.timestamp),ri(NULL){}
	~dfInfo() { if( ri ) delete ri; }

	
	dfInfo &operator =(const dfInfo &o);
	bool operator <(const dfInfo &e2) const;

	// Calculate virtual distance between two shots; different model return infinite
	double distance(const std::string &mak, const std::string &mod, int iso, double shutter) const;

	static std::string key(const std::string &mak, const std::string &mod, int iso, double shut );
	std::string key(){ return key( maker,model,iso,shutter); }

	RawImage *getRawImage();
	std::list<badPix> &getBadPixels();

protected:
	RawImage *ri; ///< Dark Frame raw data
	std::list<badPix> badPixels; ///< Unreliable pixels

	void updateBadPixelList( RawImage *df );
	void updateRawImage();
};

class DFManager
{
public:
	void init( Glib::ustring pathname );
	Glib::ustring getPathname(){ return currentPath; };
	void getStat( int &totFiles, int &totTemplate);
	RawImage *searchDarkFrame( const std::string &mak, const std::string &mod, int iso, double shut, time_t t );
	RawImage *searchDarkFrame( Glib::ustring filename );
	std::list<badPix> searchBadPixels ( const std::string &mak, const std::string &mod, int iso, double shut, time_t t );

protected:
	typedef std::multimap<std::string,dfInfo> dfList_t;
	dfList_t dfList;
	bool initialized;
	Glib::ustring currentPath;
	bool addFileInfo(const Glib::ustring &filename );
	dfInfo &find( const std::string &mak, const std::string &mod, int isospeed, double shut, time_t t );
};

extern DFManager dfm;

}

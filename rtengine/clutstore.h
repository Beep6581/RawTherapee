#ifndef CLUT_STORE_INCLUDED
#define CLUT_STORE_INCLUDED

#include <gtkmm.h>
#include "../rtgui/threadutils.h"
#include "imagefloat.h"
#include <vector>
#include <map>

namespace rtengine { 

// simple CLUT interface
class CLUT
{
public:
    virtual void getRGB( float r, float g, float b, float &outR, float &outG, float &outB ) const = 0;
    virtual Glib::ustring profile() const  = 0;
protected:
    virtual ~CLUT() {};
};

class HaldCLUT : public CLUT
{
public:
    HaldCLUT();
    ~HaldCLUT();
    void load( Glib::ustring filename );
    bool isValid() const;

    void getRGB( float r, float g, float b, float &outR, float &outG, float &outB ) const;
    Glib::ustring profile() const;

    typedef std::vector<unsigned char> RawClut; // using 8 bit for reduce memory usage
    static void correct( const RawClut&, int level, float r, float g, float b, float &outR, float &outG, float &outB );
    static void correct( Imagefloat &clutImage, int level, float rr, float gg, float bb, float &outR, float &outG, float &outB );
    static Imagefloat* generateIdentImage( int level );
    static Imagefloat* loadFile( Glib::ustring filename, Glib::ustring workingColorSpace, int &outLevel );

private:
    
    void loadClut( Imagefloat *img, RawClut &outClut );
    
    Imagefloat *m_clutImage;
    int m_level;
    Glib::ustring m_filename;
    Glib::ustring m_profile;
};

// CLUT cache
class CLUTStore
{
public:
    CLUTStore();
    
    CLUT* getClut( const Glib::ustring& filename );
    void releaseClut( const CLUT* clut );

	void clearCache();

private:
	typedef std::map<Glib::ustring, std::pair<int, HaldCLUT*> > Cluts;
	
	Cluts m_cluts;
    MyMutex m_mutex;
};

void splitClutFilename( Glib::ustring filename, Glib::ustring &name, Glib::ustring &extension, Glib::ustring &profileName );

}; //namespace rtengine 

extern rtengine::CLUTStore clutStore;

namespace rtengine {

//support class for automate call of clutStore.releaseClut()
class ClutPtr
{
public:
    ClutPtr() : m_point( 0 ) {}
    explicit ClutPtr(CLUT *p) : m_point( p ) {}
    ~ClutPtr() { clutStore.releaseClut( m_point ); }
    const CLUT* operator-> () const { return m_point; }
    operator bool() const { return m_point != 0; }
    void set( CLUT *p ) { m_point = p; }

private:
    ClutPtr& operator=(ClutPtr const& cp );
    CLUT *m_point;
};

}; //namespace rtengine 

#endif

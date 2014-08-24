#include "clutstore.h"
#include "rt_math.h"
#include "stdimagesource.h"
#include "safegtk.h"

rtengine::CLUTStore clutStore;

using namespace rtengine;

const float MAXVAL8 = 255.;

CLUTStore::CLUTStore()
:   m_refCount( 0 )
{
}

CLUT* CLUTStore::getClut( Glib::ustring filename )
{
    //MyMutex::MyLock lock(m_mutex); // copypasted from iccstore
    CLUT *result = 0;
    if ( !m_lastHaldClut.isValid() || m_lastFilename != filename )
    {
        takeUpClut();
        m_lastHaldClut.load( filename );
    }
    if ( m_lastHaldClut.isValid() )
    {
        result = &m_lastHaldClut;
        m_lastFilename = filename;
        ++m_refCount;
    }
    else
    {
        m_mutex.unlock();
    }
    return result;
}

void CLUTStore::takeUpClut()
{
    m_mutex.lock();
}

void CLUTStore::releaseClut( CLUT *clut )
{
    if ( clut && --m_refCount == 0 )
    {
        m_mutex.unlock();
    }
}

void CLUTStore::clearCache()
{
    m_lastHaldClut.clear();
    m_lastFilename.clear();
}

void rtengine::splitClutFilename( Glib::ustring filename, Glib::ustring &name, Glib::ustring &extension, Glib::ustring &profileName )
{
    filename = Glib::path_get_basename( filename );
    name = filename;
    //remove dirs
    size_t lastSlashPos = filename.find_last_of( "/" );
    if ( lastSlashPos == Glib::ustring::npos )
    {
        lastSlashPos = filename.find_last_of( "\\" );
    }

    size_t lastDotPos = filename.find_last_of( '.' );
    if ( lastDotPos != Glib::ustring::npos )
    {
        name = filename.substr( 0, lastDotPos );
        extension = filename.substr( lastDotPos + 1, Glib::ustring::npos );
    }
    profileName = "sRGB"; // sRGB by default
    static std::vector<Glib::ustring> workingProfiles = rtengine::getWorkingProfiles();
    for ( std::vector<Glib::ustring>::iterator it = workingProfiles.begin();
        it != workingProfiles.end();
        ++it )
    {
        Glib::ustring &currentProfile = *it;
        if ( std::search( name.rbegin(), name.rend(), currentProfile.rbegin(), currentProfile.rend() ) == name.rbegin() )
        {
            profileName = currentProfile;
            name = name.substr( 0, name.size() - currentProfile.size() );
            break;
        }
    }
}

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

HaldCLUT::HaldCLUT()
:   m_clutImage( 0 ),
    m_level (0),
    m_profile( "sRGB" )
{}

HaldCLUT::~HaldCLUT()
{
}

void HaldCLUT::clear()
{
    if ( m_clutImage )
    {
        m_clutImage->free();
        m_clutImage = 0;
    }
    m_filename.clear();
}

void HaldCLUT::load( Glib::ustring filename )
{
    if ( m_filename != filename )
    {
        clear();
        m_clutImage = loadFile( filename, "", m_level );
        Glib::ustring name, ext;
        splitClutFilename( filename, name, ext, m_profile );
        if ( m_clutImage )
        {
            m_filename = filename;
        }
    }
}

Glib::ustring HaldCLUT::profile() const
{
    return m_profile;
}

Imagefloat* HaldCLUT::loadFile( Glib::ustring filename, Glib::ustring workingColorSpace, int &outLevel )
{
    Imagefloat *result = 0;
    StdImageSource imgSrc;
    if ( !safe_file_test( filename, Glib::FILE_TEST_EXISTS ) || imgSrc.load(filename) ) 
    {
        return result;
    }

    int fw, fh;
    imgSrc.getFullSize (fw, fh, TR_NONE);

    bool valid = false;
    //test on Hald format, copypasted from http://www.quelsolaar.com/technology/clut.html
    if ( fw == fh )
    {
        outLevel = 1;
        for(; outLevel * outLevel * outLevel < fw; outLevel++);
        if( !( outLevel * outLevel * outLevel > fw ) )
        {
            valid = true;
        }
    }

    if ( valid )
    {
        ColorTemp currWB = imgSrc.getWB();
        Imagefloat* baseImg = new Imagefloat (fw, fh);
        PreviewProps pp (0, 0, fw, fh, 1);

        procparams::ColorManagementParams icm;
        icm.working = workingColorSpace;

        imgSrc.getImage (currWB, TR_NONE, baseImg, pp, procparams::ToneCurveParams(), icm, procparams::RAWParams());
        if ( !workingColorSpace.empty() )
        {
            imgSrc.convertColorSpace(baseImg, icm, currWB, procparams::RAWParams());
        }
        result = baseImg;
    }
    return result;  
}

void HaldCLUT::loadClut( Imagefloat *img, RawClut &outClut )
{
    img->normalizeFloatTo1();
    int y_size = img->getH();
    int x_size = img->getW();
    outClut.resize( x_size * y_size * 3 );
    int clutIdx = 0;
    //int level = m_level * m_level;  (unused)
    for(int y = 0; y < y_size; y++)
    {
        for(int x = 0; x < x_size; x++)
        {
            outClut[ clutIdx * 3     ] = img->r( y, x ) * MAXVAL8;
            outClut[ clutIdx * 3 + 1 ] = img->g( y, x ) * MAXVAL8;
            outClut[ clutIdx * 3 + 2 ] = img->b( y, x ) * MAXVAL8;

            ++clutIdx;
        }
    }
}

Imagefloat* HaldCLUT::generateIdentImage( int level )
{
    int imageWidth = level * level * level;
    Imagefloat *resultImg = new Imagefloat( imageWidth, imageWidth );
    
    int cubeSideSize = level * level;
    float step = MAXVALF / (cubeSideSize - 1);
    int pos = 0;
    for( int b = 0; b < cubeSideSize; ++b )
    {
        for ( int g = 0; g < cubeSideSize; ++g )
        {
            for ( int r = 0; r < cubeSideSize; ++r )
            {
                int x = pos / imageWidth;
                int y = pos % imageWidth;
                resultImg->r( x, y ) = step * r;
                resultImg->g( x, y ) = step * g;
                resultImg->b( x, y ) = step * b;
                ++pos;
            }
        }
    }
    return resultImg;
}


bool HaldCLUT::isValid() const
{
    return m_clutImage != 0;
}

void HaldCLUT::getRGB( float rr, float gg, float bb, float &outR, float &outG, float &outB ) const
{
    rr /= MAXVALF;
    gg /= MAXVALF;
    bb /= MAXVALF;
    correct( *m_clutImage, m_level, rr, gg, bb, outR, outG, outB );
}

inline float valF( unsigned char val )
{
    return float( val ) / MAXVAL8;
}

// copypasted from http://www.quelsolaar.com/technology/clut.html
void HaldCLUT::correct( const HaldCLUT::RawClut& clut, int level, float rr, float gg, float bb, float &outR, float &outG, float &outB )
{
    int color, red, green, blue, i, j;
    float tmp[6], r, g, b;
    level =  level * level;

    red = rr * (float)(level - 1);
    if(red > level - 2)
        red = (float)level - 2;
    if(red < 0)
        red = 0;

    green = gg * (float)(level - 1);
    if(green > level - 2)
        green = (float)level - 2;
    if(green < 0)
        green = 0;

    blue = bb * (float)(level - 1);
    if(blue > level - 2)
        blue = (float)level - 2;
    if(blue < 0)
        blue = 0;

    r = rr * (float)(level - 1) - red;
    g = gg * (float)(level - 1) - green;
    b = bb * (float)(level - 1) - blue;

    color = red + green * level + blue * level * level;

    i = color * 3;
    j = (color + 1) * 3;

    tmp[0] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[1] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[2] = valF( clut[i] ) * (1 - r) + valF( clut[j] ) * r;

    i = (color + level) * 3;
    j = (color + level + 1) * 3;

    tmp[3] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[4] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[5] = valF( clut[i] ) * (1 - r) + valF( clut[j] ) * r;

    outR = tmp[0] * (1 - g) + tmp[3] * g;
    outG = tmp[1] * (1 - g) + tmp[4] * g;
    outB = tmp[2] * (1 - g) + tmp[5] * g;

    i = (color + level * level) * 3;
    j = (color + level * level + 1) * 3;

    tmp[0] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[1] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[2] = valF( clut[i] ) * (1 - r) + valF( clut[j] ) * r;

    i = (color + level + level * level) * 3;
    j = (color + level + level * level + 1) * 3;

    tmp[3] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[4] = valF( clut[i++] ) * (1 - r) + valF( clut[j++] ) * r;
    tmp[5] = valF( clut[i] ) * (1 - r) + valF( clut[j] ) * r;

    tmp[0] = tmp[0] * (1 - g) + tmp[3] * g;
    tmp[1] = tmp[1] * (1 - g) + tmp[4] * g;
    tmp[2] = tmp[2] * (1 - g) + tmp[5] * g;

    outR = outR * (1 - b) + tmp[0] * b;
    outG = outG * (1 - b) + tmp[1] * b;
    outB = outB * (1 - b) + tmp[2] * b;
}

inline void pos2xy( int pos, int imageSideSize, int &outX, int &outY )
{
    outX = pos / imageSideSize;
    outY = pos % imageSideSize;
}

void HaldCLUT::correct( Imagefloat &clutImage, int level, float rr, float gg, float bb, float &outR, float &outG, float &outB )
{
    int color, red, green, blue, i, j;
    float tmp[6], r, g, b;
    level =  level * level;
    int imageSideSize = clutImage.getW();

    red = rr * (float)(level - 1);
    if(red > level - 2)
        red = (float)level - 2;
    if(red < 0)
        red = 0;

    green = gg * (float)(level - 1);
    if(green > level - 2)
        green = (float)level - 2;
    if(green < 0)
        green = 0;

    blue = bb * (float)(level - 1);
    if(blue > level - 2)
        blue = (float)level - 2;
    if(blue < 0)
        blue = 0;

    r = rr * (float)(level - 1) - red;
    g = gg * (float)(level - 1) - green;
    b = bb * (float)(level - 1) - blue;

    color = red + green * level + blue * level * level;


    i = color;
    j = color + 1;
    int xi, yi, xj, yj;
    pos2xy( i, imageSideSize, xi, yi );
    pos2xy( j, imageSideSize, xj, yj );

    tmp[0] = clutImage.r( xi, yi ) * (1 - r) + clutImage.r( xj, yj ) * r;
    tmp[1] = clutImage.g( xi, yi ) * (1 - r) + clutImage.g( xj, yj ) * r;
    tmp[2] = clutImage.b( xi, yi ) * (1 - r) + clutImage.b( xj, yj ) * r;

    i = color + level;
    j = color + level + 1;
    pos2xy( i, imageSideSize, xi, yi );
    pos2xy( j, imageSideSize, xj, yj );

    tmp[3] = clutImage.r( xi, yi ) * (1 - r) + clutImage.r( xj, yj ) * r;
    tmp[4] = clutImage.g( xi, yi ) * (1 - r) + clutImage.g( xj, yj ) * r;
    tmp[5] = clutImage.b( xi, yi ) * (1 - r) + clutImage.b( xj, yj ) * r;

    outR = tmp[0] * (1 - g) + tmp[3] * g;
    outG = tmp[1] * (1 - g) + tmp[4] * g;
    outB = tmp[2] * (1 - g) + tmp[5] * g;

    i = color + level * level;
    j = color + level * level + 1;
    pos2xy( i, imageSideSize, xi, yi );
    pos2xy( j, imageSideSize, xj, yj );

    tmp[0] = clutImage.r( xi, yi ) * (1 - r) + clutImage.r( xj, yj ) * r;
    tmp[1] = clutImage.g( xi, yi ) * (1 - r) + clutImage.g( xj, yj ) * r;
    tmp[2] = clutImage.b( xi, yi ) * (1 - r) + clutImage.b( xj, yj ) * r;

    i = color + level + level * level;
    j = color + level + level * level + 1;
    pos2xy( i, imageSideSize, xi, yi );
    pos2xy( j, imageSideSize, xj, yj );

    tmp[3] = clutImage.r( xi, yi ) * (1 - r) + clutImage.r( xj, yj ) * r;
    tmp[4] = clutImage.g( xi, yi ) * (1 - r) + clutImage.g( xj, yj ) * r;
    tmp[5] = clutImage.b( xi, yi ) * (1 - r) + clutImage.b( xj, yj ) * r;

    tmp[0] = tmp[0] * (1 - g) + tmp[3] * g;
    tmp[1] = tmp[1] * (1 - g) + tmp[4] * g;
    tmp[2] = tmp[2] * (1 - g) + tmp[5] * g;

    outR = outR * (1 - b) + tmp[0] * b;
    outG = outG * (1 - b) + tmp[1] * b;
    outB = outB * (1 - b) + tmp[2] * b;
}

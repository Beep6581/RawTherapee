#include "clutstore.h"
#include "rt_math.h"
#include "stdimagesource.h"
#include "../rtgui/options.h"

namespace
{

std::unique_ptr<rtengine::Imagefloat> loadFile(
    const Glib::ustring& filename,
    const Glib::ustring& working_color_space,
    unsigned int& clut_level
)
{
    std::unique_ptr<rtengine::Imagefloat> result;

    rtengine::StdImageSource img_src;

    if (!Glib::file_test(filename, Glib::FILE_TEST_EXISTS) || img_src.load(filename)) {
        return result;
    }

    int fw, fh;
    img_src.getFullSize(fw, fh, TR_NONE);

    bool valid = false;

    if (fw == fh) {
        unsigned int level = 1;
        while (level * level * level < fw) {
            ++level;
        }
        if (level * level * level == fw && level > 1) {
            clut_level = level;
            valid = true;
        }
    }

    if (valid) {
        rtengine::ColorTemp curr_wb = img_src.getWB();
        result = std::unique_ptr<rtengine::Imagefloat>(new rtengine::Imagefloat(fw, fh));
        const PreviewProps pp(0, 0, fw, fh, 1);

        rtengine::procparams::ColorManagementParams icm;
        icm.working = working_color_space;

        img_src.getImage(curr_wb, TR_NONE, result.get(), pp, rtengine::procparams::ToneCurveParams(), icm, rtengine::procparams::RAWParams());

        if (!working_color_space.empty()) {
            img_src.convertColorSpace(result.get(), icm, curr_wb);
        }
    }

    return result;
}

inline void posToXy(unsigned int pos, unsigned int width, unsigned int (&x)[2], unsigned int (&y)[2])
{
    x[0] = pos % width;
    y[0] = pos / width;
    x[1] = (pos + 1) % width;
    y[1] = (pos + 1) / width;
}

}

void rtengine::CLUT::splitClutFilename(
    const Glib::ustring& filename,
    Glib::ustring& name,
    Glib::ustring& extension,
    Glib::ustring& profile_name
)
{
    Glib::ustring basename = Glib::path_get_basename(filename);

    Glib::ustring::size_type last_slash_pos = basename.rfind('/');
    if (last_slash_pos == Glib::ustring::npos) {
        last_slash_pos = basename.rfind('\\');
    }

    const Glib::ustring::size_type last_dot_pos = basename.rfind('.');

    if (last_dot_pos != Glib::ustring::npos) {
        name.assign(basename, 0, last_dot_pos);
        extension.assign(basename, last_dot_pos + 1, Glib::ustring::npos);
    } else {
        name = basename;
    }

    profile_name = "sRGB";

    for (const auto& working_profile : rtengine::getWorkingProfiles()) {
        if ( std::search( name.rbegin(), name.rend(), working_profile.rbegin(), working_profile.rend() ) == name.rbegin() ) {
            profile_name = working_profile;
            name.erase(name.size() - working_profile.size());
            break;
        }
    }
}

rtengine::HaldCLUT::HaldCLUT() :
    clut_level(0),
    clut_profile("sRGB")
{
}

rtengine::HaldCLUT::~HaldCLUT()
{
}

bool rtengine::HaldCLUT::load(const Glib::ustring& filename)
{
    clut_image = loadFile(filename, "", clut_level);
    Glib::ustring name, ext;
    splitClutFilename(filename, name, ext, clut_profile);

    if (clut_image) {
        clut_filename = filename;
        return true;
    }

    return false;
}

rtengine::HaldCLUT::operator bool() const
{
    return static_cast<bool>(clut_image);
}

Glib::ustring rtengine::HaldCLUT::getFilename() const
{
    return clut_filename;
}

Glib::ustring rtengine::HaldCLUT::getProfile() const
{
    return clut_profile;
}

void rtengine::HaldCLUT::getRGB(float r, float g, float b, float& out_r, float& out_g, float& out_b) const
{
    const unsigned int level = clut_level * clut_level;

    const float flevel_minus_one = static_cast<float>(level - 1) / 65535.0f;
    const float flevel_minus_two = static_cast<float>(level - 2);

    const unsigned int red = std::min(flevel_minus_two, r * flevel_minus_one);
    const unsigned int green = std::min(flevel_minus_two, g * flevel_minus_one);
    const unsigned int blue = std::min(flevel_minus_two, b * flevel_minus_one);

    r = r * flevel_minus_one - red;
    g = g * flevel_minus_one - green;
    b = b * flevel_minus_one - blue;

    const unsigned int level_square = level * level;

    const unsigned int color = red + green * level + blue * level_square;

    unsigned int x[2];
    unsigned int y[2];
    posToXy(color, clut_image->getWidth(), x, y);

    float tmp1[4] __attribute__((aligned(16)));
    tmp1[0] = clut_image->r(y[0], x[0]) * (1 - r) + clut_image->r(y[1], x[1]) * r;
    tmp1[1] = clut_image->g(y[0], x[0]) * (1 - r) + clut_image->g(y[1], x[1]) * r;
    tmp1[2] = clut_image->b(y[0], x[0]) * (1 - r) + clut_image->b(y[1], x[1]) * r;

    posToXy(color + level, clut_image->getWidth(), x, y);

    float tmp2[4] __attribute__((aligned(16)));
    tmp2[0] = clut_image->r(y[0], x[0]) * (1 - r) + clut_image->r(y[1], x[1]) * r;
    tmp2[1] = clut_image->g(y[0], x[0]) * (1 - r) + clut_image->g(y[1], x[1]) * r;
    tmp2[2] = clut_image->b(y[0], x[0]) * (1 - r) + clut_image->b(y[1], x[1]) * r;

    float out[4] __attribute__((aligned(16)));
    out[0] = tmp1[0] * (1 - g) + tmp2[0] * g;
    out[1] = tmp1[1] * (1 - g) + tmp2[1] * g;
    out[2] = tmp1[2] * (1 - g) + tmp2[2] * g;

    posToXy(color + level_square, clut_image->getWidth(), x, y);

    tmp1[0] = clut_image->r(y[0], x[0]) * (1 - r) + clut_image->r(y[1], x[1]) * r;
    tmp1[1] = clut_image->g(y[0], x[0]) * (1 - r) + clut_image->g(y[1], x[1]) * r;
    tmp1[2] = clut_image->b(y[0], x[0]) * (1 - r) + clut_image->b(y[1], x[1]) * r;

    posToXy(color + level + level_square, clut_image->getWidth(), x, y);

    tmp2[0] = clut_image->r(y[0], x[0]) * (1 - r) + clut_image->r(y[1], x[1]) * r;
    tmp2[1] = clut_image->g(y[0], x[0]) * (1 - r) + clut_image->g(y[1], x[1]) * r;
    tmp2[2] = clut_image->b(y[0], x[0]) * (1 - r) + clut_image->b(y[1], x[1]) * r;

    tmp1[0] = tmp1[0] * (1 - g) + tmp2[0] * g;
    tmp1[1] = tmp1[1] * (1 - g) + tmp2[1] * g;
    tmp1[2] = tmp1[2] * (1 - g) + tmp2[2] * g;

    out_r = out[0] * (1 - b) + tmp1[0] * b;
    out_g = out[1] * (1 - b) + tmp1[1] * b;
    out_b = out[2] * (1 - b) + tmp1[2] * b;
}

rtengine::CLUTStore& rtengine::CLUTStore::getInstance()
{
    static CLUTStore instance;
    return instance;
}

std::shared_ptr<rtengine::CLUT> rtengine::CLUTStore::getClut(const Glib::ustring& filename)
{
    std::shared_ptr<rtengine::CLUT> result;

    if (!cache.get(filename, result)) {
        std::unique_ptr<rtengine::HaldCLUT> clut(new rtengine::HaldCLUT);
        if (clut->load(filename)) {
            result = std::move(clut);
            cache.insert(filename, result);
        }
    }

    return result;
}

void rtengine::CLUTStore::releaseClut(const std::shared_ptr<rtengine::CLUT>& clut)
{
    cache.remove(clut->getFilename());
}

void rtengine::CLUTStore::clearCache()
{
    cache.clear();
}

rtengine::CLUTStore::CLUTStore() :
    cache(options.clutCacheSize)
{
}

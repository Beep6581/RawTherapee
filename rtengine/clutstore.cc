#include <algorithm>

#ifdef __SSE2__
#include <xmmintrin.h>
#endif

#include "clutstore.h"

#include "imagefloat.h"
#include "stdimagesource.h"
#include "../rtgui/options.h"

namespace
{

bool loadFile(
    const Glib::ustring& filename,
    const Glib::ustring& working_color_space,
    AlignedBuffer<std::uint16_t>& clut_image,
    unsigned int& clut_level
)
{
    rtengine::StdImageSource img_src;

    if (!Glib::file_test(filename, Glib::FILE_TEST_EXISTS) || img_src.load(filename)) {
        return false;
    }

    int fw, fh;
    img_src.getFullSize(fw, fh, TR_NONE);

    bool res = false;

    if (fw == fh) {
        unsigned int level = 1;
        while (level * level * level < fw) {
            ++level;
        }
        if (level * level * level == fw && level > 1) {
            clut_level = level;
            res = true;
        }
    }

    if (res) {
        rtengine::ColorTemp curr_wb = img_src.getWB();
        std::unique_ptr<rtengine::Imagefloat> img_float = std::unique_ptr<rtengine::Imagefloat>(new rtengine::Imagefloat(fw, fh));
        const PreviewProps pp(0, 0, fw, fh, 1);

        rtengine::procparams::ColorManagementParams icm;
        icm.working = working_color_space;

        img_src.getImage(curr_wb, TR_NONE, img_float.get(), pp, rtengine::procparams::ToneCurveParams(), icm, rtengine::procparams::RAWParams());

        if (!working_color_space.empty()) {
            img_src.convertColorSpace(img_float.get(), icm, curr_wb);
        }

        AlignedBuffer<std::uint16_t> image(fw * fh * 4 + 1);

        std::size_t index = 0;
        for (int y = 0; y < fh; ++y) {
            for (int x = 0; x < fw; ++x) {
                image.data[index] = img_float->r(y, x);
                ++index;
                image.data[index] = img_float->g(y, x);
                ++index;
                image.data[index] = img_float->b(y, x);
                index += 2;
            }
        }

        clut_image.swap(image);
    }

    return res;
}

inline void posToIndex(unsigned int pos, size_t (&index)[2])
{
    index[0] = static_cast<size_t>(pos) * 4;
    index[1] = static_cast<size_t>(pos + 1) * 4;
}

#ifdef __SSE2__
inline __m128 getClutValue(const AlignedBuffer<std::uint16_t>& clut_image, size_t index)
{
    return _mm_cvtpu16_ps(*reinterpret_cast<const __m64*>(clut_image.data + index));
}
#endif

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
    flevel_minus_one(0.0f),
    flevel_minus_two(0.0f),
    clut_profile("sRGB")
{
}

rtengine::HaldCLUT::~HaldCLUT()
{
}

bool rtengine::HaldCLUT::load(const Glib::ustring& filename)
{
    if (loadFile(filename, "", clut_image, clut_level)) {
        Glib::ustring name, ext;
        splitClutFilename(filename, name, ext, clut_profile);

        clut_filename = filename;
        clut_level *= clut_level;
        flevel_minus_one = static_cast<float>(clut_level - 1) / 65535.0f;
        flevel_minus_two = static_cast<float>(clut_level - 2);
        return true;
    }

    return false;
}

rtengine::HaldCLUT::operator bool() const
{
    return !clut_image.isEmpty();
}

Glib::ustring rtengine::HaldCLUT::getFilename() const
{
    return clut_filename;
}

Glib::ustring rtengine::HaldCLUT::getProfile() const
{
    return clut_profile;
}

void rtengine::HaldCLUT::getRGB(float r, float g, float b, float out_rgbx[4]) const
{
	const unsigned int level = clut_level; // This is important

    const unsigned int red = std::min(flevel_minus_two, r * flevel_minus_one);
    const unsigned int green = std::min(flevel_minus_two, g * flevel_minus_one);
    const unsigned int blue = std::min(flevel_minus_two, b * flevel_minus_one);

	const unsigned int level_square = level * level;

	const unsigned int color = red + green * level + blue * level_square;

#ifndef __SSE2__
    r = r * flevel_minus_one - red;
    g = g * flevel_minus_one - green;
    b = b * flevel_minus_one - blue;

	size_t index[2];
	posToIndex(color, index);

	float tmp1[4] ALIGNED16;
	tmp1[0] = clut_image.data[index[0]] * (1 - r) + clut_image.data[index[1]] * r;
	tmp1[1] = clut_image.data[index[0] + 1] * (1 - r) + clut_image.data[index[1] + 1] * r;
	tmp1[2] = clut_image.data[index[0] + 2] * (1 - r) + clut_image.data[index[1] + 2] * r;

	posToIndex(color + level, index);

	float tmp2[4] ALIGNED16;
	tmp2[0] = clut_image.data[index[0]] * (1 - r) + clut_image.data[index[1]] * r;
	tmp2[1] = clut_image.data[index[0] + 1] * (1 - r) + clut_image.data[index[1] + 1] * r;
	tmp2[2] = clut_image.data[index[0] + 2] * (1 - r) + clut_image.data[index[1] + 2] * r;

	out_rgbx[0] = tmp1[0] * (1 - g) + tmp2[0] * g;
	out_rgbx[1] = tmp1[1] * (1 - g) + tmp2[1] * g;
	out_rgbx[2] = tmp1[2] * (1 - g) + tmp2[2] * g;

	posToIndex(color + level_square, index);

	tmp1[0] = clut_image.data[index[0]] * (1 - r) + clut_image.data[index[1]] * r;
	tmp1[1] = clut_image.data[index[0] + 1] * (1 - r) + clut_image.data[index[1] + 1] * r;
	tmp1[2] = clut_image.data[index[0] + 2] * (1 - r) + clut_image.data[index[1] + 2] * r;

	posToIndex(color + level + level_square, index);

	tmp2[0] = clut_image.data[index[0]] * (1 - r) + clut_image.data[index[1]] * r;
	tmp2[1] = clut_image.data[index[0] + 1] * (1 - r) + clut_image.data[index[1] + 1] * r;
	tmp2[2] = clut_image.data[index[0] + 2] * (1 - r) + clut_image.data[index[1] + 2] * r;

	tmp1[0] = tmp1[0] * (1 - g) + tmp2[0] * g;
	tmp1[1] = tmp1[1] * (1 - g) + tmp2[1] * g;
	tmp1[2] = tmp1[2] * (1 - g) + tmp2[2] * g;

    out_rgbx[0] = out_rgbx[0] * (1 - b) + tmp1[0] * b;
    out_rgbx[1] = out_rgbx[1] * (1 - b) + tmp1[1] * b;
    out_rgbx[2] = out_rgbx[2] * (1 - b) + tmp1[2] * b;
#else
	const __m128 v_rgb = _mm_set_ps(0.0f, b, g, r) *_mm_load_ps1(&flevel_minus_one) - _mm_set_ps(0.0f, blue, green, red);

	size_t index[2];
	posToIndex(color, index);

	const __m128 v_r = _mm_shuffle_ps(v_rgb, v_rgb, 0x00);

    __m128 v_cv0 = getClutValue(clut_image, index[0]);
    __m128 v_tmp1 = v_r * (getClutValue(clut_image, index[1]) - v_cv0) + v_cv0;

	posToIndex(color + level, index);

    v_cv0 = getClutValue(clut_image, index[0]);
	__m128 v_tmp2 = v_r * (getClutValue(clut_image, index[1]) - v_cv0) + v_cv0;

	const __m128 v_g = _mm_shuffle_ps(v_rgb, v_rgb, 0x55);

	__m128 v_out = v_g * (v_tmp2 - v_tmp1) + v_tmp1;

	posToIndex(color + level_square, index);

    v_cv0 = getClutValue(clut_image, index[0]);
	v_tmp1 = v_r * (getClutValue(clut_image, index[1]) - v_cv0) + v_cv0;

	posToIndex(color + level + level_square, index);

    v_cv0 = getClutValue(clut_image, index[0]);
	v_tmp2 = v_r * (getClutValue(clut_image, index[1]) - v_cv0) + v_cv0;

	v_tmp1 = v_g * (v_tmp2 - v_tmp1) + v_tmp1;

	const __m128 v_b = _mm_shuffle_ps(v_rgb, v_rgb, 0xAA);

	_mm_store_ps(out_rgbx, v_b * (v_tmp1 - v_out) + v_out);
#endif
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

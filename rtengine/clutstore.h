#pragma once

#include <memory>

#include <gtkmm.h>

#include "imagefloat.h"
#include "cache.h"

namespace rtengine
{

class CLUT
{
public:
    CLUT() = default;
    CLUT(const CLUT& other) = delete;
    CLUT& operator =(const CLUT& other) = delete;
    virtual ~CLUT() = default;

    virtual explicit operator bool() const = 0;

    virtual Glib::ustring getFilename() const = 0;
    virtual Glib::ustring getProfile() const = 0;

    virtual void getRGB(float r, float g, float b, float& out_r, float& out_g, float& out_b) const = 0;

    static void splitClutFilename(
        const Glib::ustring& filename,
        Glib::ustring& name,
        Glib::ustring& extension,
        Glib::ustring& profile_name
    );
};

class HaldCLUT
    : public CLUT
{
public:
    HaldCLUT();
    ~HaldCLUT();

    bool load(const Glib::ustring& filename);

    explicit operator bool() const;

    Glib::ustring getFilename() const;
    Glib::ustring getProfile() const;

    void getRGB(float r, float g, float b, float& out_r, float& out_g, float& out_b) const;

private:
    std::unique_ptr<Imagefloat> clut_image;
    unsigned int clut_level;
    Glib::ustring clut_filename;
    Glib::ustring clut_profile;
};

class CLUTStore
{
public:
    static CLUTStore& getInstance();

    CLUTStore(const CLUTStore& other) = delete;
    CLUTStore& operator =(const CLUTStore& other) = delete;

    std::shared_ptr<CLUT> getClut(const Glib::ustring& filename);
    void releaseClut(const std::shared_ptr<CLUT>& clut);

    void clearCache();

private:
    CLUTStore();

    Cache<Glib::ustring, std::shared_ptr<CLUT>> cache;
};

}

/*
 *  This file is part of RawTherapee.
 */
#pragma once

#include <map>
#include <string>
#include <vector>

namespace Glib
{

class ustring;

}

namespace rtengine
{

struct camera_const_levels {
    int levels[4];
};

class CameraConst final
{
private:
    std::string make_model;
    short dcraw_matrix[12];
    int raw_crop[4];
    int raw_mask[2][4];
    int white_max;
    std::map<int, camera_const_levels> mLevels[2];
    std::map<float, float> mApertureScaling;
    std::vector<int> pdafPattern;
    int pdafOffset;
    int globalGreenEquilibration;
    
    CameraConst();
    static bool parseLevels(CameraConst *cc, int bw, const void *ji);
    static bool parseApertureScaling(CameraConst *cc, const void *ji);
    bool get_Levels(camera_const_levels & lvl, int bw, int iso, float fnumber) const;

public:
    static CameraConst *parseEntry(const void *cJSON, const char *make_model);
    bool has_dcrawMatrix(void) const;
    void update_dcrawMatrix(const short *other);
    const short *get_dcrawMatrix(void) const;
    const std::vector<int>& get_pdafPattern() const;
    int get_pdafOffset() const {return pdafOffset;};
    bool has_rawCrop(void) const;
    void get_rawCrop(int& left_margin, int& top_margin, int& width, int& height) const;
    bool has_rawMask(int idx) const;
    void get_rawMask(int idx, int& top, int& left, int& bottom, int& right) const;
    int get_BlackLevel(int idx, int iso_speed) const;
    int get_WhiteLevel(int idx, int iso_speed, float fnumber) const;
    bool has_globalGreenEquilibration() const;
    bool get_globalGreenEquilibration() const;
    void update_Levels(const CameraConst *other);
    void update_Crop(CameraConst *other);
    void update_pdafPattern(const std::vector<int> &other);
    void update_pdafOffset(int other);
    void update_globalGreenEquilibration(bool other);
};

class CameraConstantsStore final
{
private:
    std::map<std::string, CameraConst *> mCameraConstants;

    CameraConstantsStore();
    bool parse_camera_constants_file(const Glib::ustring& filename);

public:
    ~CameraConstantsStore();
    void init(const Glib::ustring& baseDir, const Glib::ustring& userSettingsDir);
    static CameraConstantsStore *getInstance(void);
    const CameraConst *get(const char make[], const char model[]) const;
};

}

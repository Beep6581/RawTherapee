/*
 *  This file is part of RawTherapee.
 */
#ifndef __CAMCONST__
#define __CAMCONST__

#include <glibmm.h>
#include <map>

namespace rtengine
{

struct camera_const_levels {
    int levels[4];
};

class CameraConst
{
private:
    std::string make_model;
    short dcraw_matrix[12];
    int raw_crop[4];
    int raw_mask[8][4];
    int white_max;
    std::map<int, struct camera_const_levels> mLevels[2];
    std::map<float, float> mApertureScaling;
    std::vector<int> pdafPattern;
    int pdafOffset;
    CameraConst();
    static bool parseLevels(CameraConst *cc, int bw, void *ji);
    static bool parseApertureScaling(CameraConst *cc, void *ji);
    bool get_Levels(struct camera_const_levels & lvl, int bw, int iso, float fnumber);

public:
    static CameraConst *parseEntry(void *cJSON, const char *make_model);
    bool has_dcrawMatrix(void);
    bool has_pdafPattern(void);
    void update_dcrawMatrix(const short *other);
    const short *get_dcrawMatrix(void);
    void get_pdafPattern(std::vector<int> &pattern);
    int get_pdafOffset() {return pdafOffset;}
    bool has_rawCrop(void);
    void get_rawCrop(int& left_margin, int& top_margin, int& width, int& height);
    bool has_rawMask(int idx);
    void get_rawMask(int idx, int& top, int& left, int& bottom, int& right);
    int get_BlackLevel(int idx, int iso_speed);
    int get_WhiteLevel(int idx, int iso_speed, float fnumber);
    void update_Levels(const CameraConst *other);
    void update_Crop(CameraConst *other);
};

class CameraConstantsStore
{
private:
    std::map<std::string, CameraConst *> mCameraConstants;

    CameraConstantsStore();
    bool parse_camera_constants_file(Glib::ustring filename);

public:
    ~CameraConstantsStore();
    void init(Glib::ustring baseDir, Glib::ustring userSettingsDir);
    static CameraConstantsStore *getInstance(void);
    CameraConst *get(const char make[], const char model[]);
};

}

#endif

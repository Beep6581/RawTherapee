/*
 *  This file is part of RawTherapee.
 */
#ifndef __CAMCONST__
#define __CAMCONST__

#include <glibmm.h>
#include <map>

namespace rtengine {

struct camera_const_levels {
	int levels[4];
};

class CameraConst {
  private:
    Glib::ustring make_model;
    short dcraw_matrix[12];
    int white_max;
    std::map<int, struct camera_const_levels> mLevels[2];
    std::map<float, float> mApertureScaling;

    CameraConst();
    ~CameraConst();
    static bool parseLevels(CameraConst *cc, int bw, void *ji);
    static bool parseApertureScaling(CameraConst *cc, void *ji);
    bool get_Levels(struct camera_const_levels & lvl, int bw, int iso, float fnumber);

  public:
    static CameraConst *parseEntry(void *cJSON, const char *make_model);
    bool has_dcrawMatrix(void);
    void update_dcrawMatrix(const short *other);
    const short *get_dcrawMatrix(void);
    int get_BlackLevel(int idx, int iso_speed);
    int get_WhiteLevel(int idx, int iso_speed, float fnumber);
    void update_Levels(const CameraConst *other);
};

class CameraConstantsStore {
  private:
    std::map<Glib::ustring, CameraConst *> mCameraConstants;

    CameraConstantsStore();
    bool parse_camera_constants_file(Glib::ustring filename);

  public:
    static void initCameraConstants(Glib::ustring baseDir, Glib::ustring userSettingsDir);
    static CameraConstantsStore *getInstance(void);
    CameraConst *get(const char make[], const char model[]);
};

}

#endif

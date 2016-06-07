/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <www.oliverduis.de>
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

#pragma once

#include <map>
#include <vector>

#include <glibmm.h>

#include "imagefloat.h"
#include "curves.h"
#include "colortemp.h"

#include "../rtgui/threadutils.h"

namespace rtengine
{

class DCPProfile final
{
public:
    struct ApplyState {
        double m2ProPhoto[3][3];
        double m2Work[3][3];
        bool alreadyProPhoto;
        bool useToneCurve;
        bool applyLookTable;
        float blScale;
    };

    struct Illuminants {
        short light_source_1;
        short light_source_2;
        double temperature_1;
        double temperature_2;
        bool will_interpolate;
    };

    DCPProfile(const Glib::ustring& filename);
    ~DCPProfile();

    bool getHasToneCurve() const;
    bool getHasLookTable() const;
    bool getHasHueSatMap() const;
    bool getHasBaselineExposureOffset() const;

    Illuminants getIlluminants() const;

    void apply(
        Imagefloat* img,
        int preferred_illuminant,
        const Glib::ustring& working_space,
        const ColorTemp& white_balance,
        double pre_mul[3],
        double cam_matrix[3][3],
        bool use_tone_curve = false,
        bool apply_hue_sat_map = true,
        bool apply_look_table = false
    ) const;
    void setStep2ApplyState(const Glib::ustring &workingSpace, bool useToneCurve, bool applyLookTable, bool applyBaselineExposure, ApplyState &asOut);
    void step2ApplyTile(float *r, float *g, float *b, int width, int height, int tileWidth, const ApplyState &asIn) const;

private:
    struct HSBModify {
        float hue_shift;
        float sat_scale;
        float val_scale;
    };

    struct HSDTableInfo {
        int hue_divisions;
        int sat_divisions;
        int val_divisions;
        int hue_step;
        int val_step;
        unsigned int array_count;
        bool srgb_gamma;
        struct {
            float h_scale;
            float s_scale;
            float v_scale;
            int max_hue_index0;
            int max_sat_index0;
            int max_val_index0;
            int hue_step;
            int val_step;
        } pc;
    };

    void dngref_XYCoord2Temperature(const double whiteXY[2], double *temp, double *tint) const;
    void dngref_FindXYZtoCamera(const double whiteXY[2], int preferredIlluminant, double (*xyzToCamera)[3]) const;
    void dngref_NeutralToXY(double neutral[3], int preferredIlluminant, double XY[2]) const;
    void makeXyzCam(const ColorTemp &wb, double pre_mul[3], double camWbMatrix[3][3], int preferredIlluminant, double (*mXYZCAM)[3]) const;
    std::vector<HSBModify> makeHueSatMap(const ColorTemp& white_balance, int preferred_illuminant) const;
    void hsdApply(const HSDTableInfo& table_info, const std::vector<HSBModify>& table_base, float& h, float& s, float& v) const;

    double color_matrix_1[3][3];
    double color_matrix_2[3][3];
    bool has_color_matrix_1;
    bool has_color_matrix_2;
    bool has_forward_matrix_1;
    bool has_forward_matrix_2;
    bool has_tone_curve;
    bool has_baseline_exposure_offset;
    bool will_interpolate;
    double forward_matrix_1[3][3];
    double forward_matrix_2[3][3];
    double temperature_1;
    double temperature_2;
    double baseline_exposure_offset;
    std::vector<HSBModify> deltas_1;
    std::vector<HSBModify> deltas_2;
    std::vector<HSBModify> look_table;
    HSDTableInfo delta_info;
    HSDTableInfo look_info;
    short light_source_1;
    short light_source_2;

    AdobeToneCurve tone_curve;
};

class DCPStore final
{
public:
    static DCPStore* getInstance();

    DCPStore(const DCPStore& other) = delete;
    DCPStore& operator =(const DCPStore& other) = delete;

    void init(const Glib::ustring& rt_profile_dir);

    bool isValidDCPFileName(const Glib::ustring& filename) const;

    DCPProfile* getProfile(const Glib::ustring& filename) const;
    DCPProfile* getStdProfile(const Glib::ustring& camShortName) const;

private:
    DCPStore() = default;

    mutable MyMutex mutex;

    // these contain standard profiles from RT. keys are all in uppercase, file path is value
    std::map<Glib::ustring, Glib::ustring> file_std_profiles;

    // Maps file name to profile as cache
    mutable std::map<Glib::ustring, DCPProfile*> profile_cache;
};

}

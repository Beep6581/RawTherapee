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
#include <array>
#include <memory>

#include <glibmm.h>

#include "../rtgui/threadutils.h"

#include "imagefloat.h"
#include "curves.h"
#include "colortemp.h"
#include "noncopyable.h"

namespace rtengine
{

class DCPProfile final
{
public:
    class ApplyState final
    {
    public:
        ApplyState();
        ~ApplyState();

    private:
        struct Data;

        std::unique_ptr<Data> data;

        friend class DCPProfile;
    };

    struct Illuminants {
        short light_source_1;
        short light_source_2;
        double temperature_1;
        double temperature_2;
        bool will_interpolate;
    };

    using Triple = std::array<double, 3>;
    using Matrix = std::array<Triple, 3>;

    explicit DCPProfile(const Glib::ustring& filename);
    ~DCPProfile();

    explicit operator bool() const;

    bool getHasToneCurve() const;
    bool getHasLookTable() const;
    bool getHasHueSatMap() const;
    bool getHasBaselineExposureOffset() const;

    Illuminants getIlluminants() const;
    bool isValid();

    void apply(
        Imagefloat* img,
        int preferred_illuminant,
        const Glib::ustring& working_space,
        const ColorTemp& white_balance,
        const Triple& pre_mul,
        const Matrix& cam_wb_matrix,
        bool apply_hue_sat_map = true
    ) const;
    void setStep2ApplyState(const Glib::ustring& working_space, bool use_tone_curve, bool apply_look_table, bool apply_baseline_exposure, ApplyState& as_out);
    void step2ApplyTile(float* r, float* g, float* b, int width, int height, int tile_width, const ApplyState& as_in) const;

private:
    struct HsbModify {
        float hue_shift;
        float sat_scale;
        float val_scale;
    };

    struct HsdTableInfo {
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

    Matrix findXyztoCamera(const std::array<double, 2>& white_xy, int preferred_illuminant) const;
    std::array<double, 2> neutralToXy(const Triple& neutral, int preferred_illuminant) const;
    Matrix makeXyzCam(const ColorTemp& white_balance, const Triple& pre_mul, const Matrix& cam_wb_matrix, int preferred_illuminant) const;
    std::vector<HsbModify> makeHueSatMap(const ColorTemp& white_balance, int preferred_illuminant) const;
    void hsdApply(const HsdTableInfo& table_info, const std::vector<HsbModify>& table_base, float& h, float& s, float& v) const;

    Matrix color_matrix_1;
    Matrix color_matrix_2;
    bool has_color_matrix_1;
    bool has_color_matrix_2;
    bool has_forward_matrix_1;
    bool has_forward_matrix_2;
    bool has_tone_curve;
    bool has_baseline_exposure_offset;
    bool will_interpolate;
    bool valid;
    Matrix forward_matrix_1;
    Matrix forward_matrix_2;
    double temperature_1;
    double temperature_2;
    double baseline_exposure_offset;
    std::vector<HsbModify> deltas_1;
    std::vector<HsbModify> deltas_2;
    std::vector<HsbModify> look_table;
    HsdTableInfo delta_info;
    HsdTableInfo look_info;
    short light_source_1;
    short light_source_2;

    AdobeToneCurve tone_curve;
};

class DCPStore final :
    public NonCopyable
{
public:
    ~DCPStore();
    static DCPStore* getInstance();

    void init(const Glib::ustring& rt_profile_dir, bool loadAll = true);

    bool isValidDCPFileName(const Glib::ustring& filename) const;

    DCPProfile* getProfile(const Glib::ustring& filename) const;
    DCPProfile* getStdProfile(const Glib::ustring& camShortName) const;

private:
    DCPStore() = default;

    mutable MyMutex mutex;
    std::vector<Glib::ustring> profileDir;

    // these contain standard profiles from RT. keys are all in uppercase, file path is value
    std::map<Glib::ustring, Glib::ustring> file_std_profiles;

    // Maps file name to profile as cache
    mutable std::map<Glib::ustring, DCPProfile*> profile_cache;
};

}

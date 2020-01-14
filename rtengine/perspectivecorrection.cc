/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

// taken from darktable (src/iop/ashift.c)
/*
  This file is part of darktable,
  copyright (c) 2016 Ulrich Pegelow.

  darktable is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  darktable is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/
// Inspiration to this module comes from the program ShiftN (http://www.shiftn.de) by
// Marcus Hebel.

// Thanks to Marcus for his support when implementing part of the ShiftN functionality
// to darktable.


#include "perspectivecorrection.h"
#include "improcfun.h"
#include "rt_math.h"
#include <string.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../rtgui/threadutils.h"
#include "settings.h"

namespace rtengine { extern const Settings *settings; }

#define _(msg) (msg)
#define dt_control_log(msg) \
    if (settings->verbose) { \
        printf("%s\n", msg);       \
        fflush(stdout);            \
    }


namespace rtengine {

namespace {

inline int mat3inv(float *const dst, const float *const src)
{
    std::array<std::array<float, 3>, 3> tmpsrc;
    std::array<std::array<float, 3>, 3> tmpdst;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            tmpsrc[i][j] = src[3 * i + j];
        }
    }
    if (invertMatrix(tmpsrc, tmpdst)) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dst[3 * i + j] = tmpdst[i][j];
            }
        }
        return 0;
    } else {
        return 1;
    }
}


// the darktable ashift iop (adapted to RT), which does most of the work
#include "ashift_dt.c"


procparams::PerspectiveParams import_meta(const procparams::PerspectiveParams &pp, const FramesMetaData *metadata)
{
    procparams::PerspectiveParams ret(pp);
    if (metadata && ret.flength == 0) {
        double f = metadata->getFocalLen();
        double f35 = metadata->getFocalLen35mm();
        if (f > 0 && f35 > 0) {
            ret.flength = f;
            ret.cropfactor = f35 / f;
        } else if (f > 0) {
            ret.flength = f;
        }
    }
    return ret;
}

} // namespace


PerspectiveCorrection::PerspectiveCorrection():
    ok_(false),
    scale_(1.0),
    offx_(0.0),
    offy_(0.0)
{
}


void PerspectiveCorrection::init(int width, int height, const procparams::PerspectiveParams &params, bool fill, const FramesMetaData *metadata)
{
    if (params.enabled) {
        auto pp = import_meta(params, metadata);
        homography((float *)ihomograph_, params.angle, params.vertical / 100.0, -params.horizontal / 100.0, params.shear / 100.0, params.flength * params.cropfactor, 100.f, params.aspect, width, height, ASHIFT_HOMOGRAPH_INVERTED);

        ok_ = true;
        calc_scale(width, height, pp, fill);
    } else {
        ok_ = false;
    }
}


inline void PerspectiveCorrection::correct(double &x, double &y, double scale, double offx, double offy)
{
    if (ok_) {       
        float pin[3], pout[3];
        pout[0] = x;
        pout[1] = y;
        pout[0] *= scale;
        pout[1] *= scale;
        pout[0] += offx;
        pout[1] += offy;
        pout[2] = 1.f;
        mat3mulv(pin, (float *)ihomograph_, pout);
        pin[0] /= pin[2];
        pin[1] /= pin[2];
        x = pin[0];
        y = pin[1];
    }
}


void PerspectiveCorrection::operator()(double &x, double &y)
{
    correct(x, y, scale_, offx_, offy_);
}


namespace {

std::vector<Coord2D> get_corners(int w, int h)
{
    int x1 = 0, y1 = 0;
    int x2 = w, y2 = h;

    std::vector<Coord2D> corners = {
        Coord2D(x1, y1),
        Coord2D(x1, y2),
        Coord2D(x2, y2),
        Coord2D(x2, y1)
    };
    return corners;
}

void init_dt_structures(dt_iop_ashift_params_t *p, dt_iop_ashift_gui_data_t *g,
                        const procparams::PerspectiveParams *params)
{
    dt_iop_ashift_params_t dp = {
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        DEFAULT_F_LENGTH,
        1.f,
        0.0f,
        1.0f,
        ASHIFT_MODE_SPECIFIC,
        0,
        ASHIFT_CROP_OFF,
        0.0f,
        1.0f,
        0.0f,
        1.0f
    };
    *p = dp;

    g->buf = NULL;
    g->buf_width = 0;
    g->buf_height = 0;
    g->buf_x_off = 0;
    g->buf_y_off = 0;
    g->buf_scale = 1.0f;
    g->buf_hash = 0;
    g->isflipped = 0;
    g->lastfit = ASHIFT_FIT_NONE;
    g->fitting = 0;
    g->lines = NULL;
    g->lines_count =0;
    g->horizontal_count = 0;
    g->vertical_count = 0;
    g->grid_hash = 0;
    g->lines_hash = 0;
    g->rotation_range = ROTATION_RANGE_SOFT;
    g->lensshift_v_range = LENSSHIFT_RANGE_SOFT;
    g->lensshift_h_range = LENSSHIFT_RANGE_SOFT;
    g->shear_range = SHEAR_RANGE_SOFT;
    g->lines_suppressed = 0;
    g->lines_version = 0;
    g->show_guides = 0;
    g->isselecting = 0;
    g->isdeselecting = 0;
    g->isbounding = ASHIFT_BOUNDING_OFF;
    g->near_delta = 0;
    g->selecting_lines_version = 0;
    g->points = NULL;
    g->points_idx = NULL;
    g->points_lines_count = 0;
    g->points_version = 0;
    g->jobcode = ASHIFT_JOBCODE_NONE;
    g->jobparams = 0;
    g->adjust_crop = FALSE;
    g->lastx = g->lasty = -1.0f;
    g->crop_cx = g->crop_cy = 1.0f;

    if (params) {
        p->rotation = params->angle;
        p->lensshift_v = params->vertical / 100.0;
        p->lensshift_h = -params->horizontal / 100.0;
        p->shear = params->shear / 100.0;
        p->f_length = params->flength;
        p->crop_factor = params->cropfactor;
        p->aspect = params->aspect;
    }
}


void get_view_size(int w, int h, const procparams::PerspectiveParams &params, double &cw, double &ch)
{
    double min_x = RT_INFINITY, max_x = -RT_INFINITY;
    double min_y = RT_INFINITY, max_y = -RT_INFINITY;

    auto corners = get_corners(w, h);

    float homo[3][3];
    homography((float *)homo, params.angle, params.vertical / 100.0, -params.horizontal / 100.0, params.shear / 100.0, params.flength * params.cropfactor, 100.f, params.aspect, w, h, ASHIFT_HOMOGRAPH_FORWARD);
    
    for (auto &c : corners) {
        float pin[3] = { float(c.x), float(c.y), 1.f };
        float pout[3];
        mat3mulv(pout, (float *)homo, pin);
        double x = pout[0] / pout[2];
        double y = pout[1] / pout[2];
        min_x = min(min_x, x);
        max_x = max(max_x, x);
        min_y = min(min_y, y);
        max_y = max(max_y, y);
    }

    cw = max_x - min_x;
    ch = max_y - min_y;
}    

} // namespace


void PerspectiveCorrection::calc_scale(int w, int h, const procparams::PerspectiveParams &params, bool fill)
{
    double cw, ch;
    get_view_size(w, h, params, cw, ch);

    if (!fill) {
        scale_ = max(cw / double(w), ch / double(h));
        offx_ = (cw - w * scale_) * 0.5;
        offy_ = (ch - h * scale_) * 0.5;
    } else {
        dt_iop_ashift_params_t p;
        dt_iop_ashift_gui_data_t g;
        init_dt_structures(&p, &g, &params);
        dt_iop_module_t module = { &g, false };
        g.buf_width = w;
        g.buf_height = h;
        p.cropmode = ASHIFT_CROP_ASPECT;
        do_crop(&module, &p);
        offx_ = p.cl * cw;
        offy_ = p.ct * ch;
        scale_ = (p.cr - p.cl) * cw/double(w);
    }
}


procparams::PerspectiveParams PerspectiveCorrection::autocompute(ImageSource *src, Direction dir, const procparams::ProcParams *pparams, const FramesMetaData *metadata)
{
    auto pcp = import_meta(pparams->perspective, metadata);
    procparams::PerspectiveParams dflt;
    pcp.horizontal = dflt.horizontal;
    pcp.vertical = dflt.vertical;
    pcp.angle = dflt.angle;
    pcp.shear = dflt.shear;
    
    dt_iop_ashift_params_t p;
    dt_iop_ashift_gui_data_t g;
    init_dt_structures(&p, &g, &pcp);
    dt_iop_module_t module;
    module.gui_data = &g;
    module.is_raw = src->isRAW();

    int tr = getCoarseBitMask(pparams->coarse);
    int fw, fh;
    src->getFullSize(fw, fh, tr);
    int skip = max(float(max(fw, fh)) / 900.f + 0.5f, 1.f);
    PreviewProps pp(0, 0, fw, fh, skip);
    int w, h;
    src->getSize(pp, w, h);
    std::unique_ptr<Imagefloat> img(new Imagefloat(w, h));

    ProcParams neutral;
    neutral.raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
    neutral.raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
    neutral.icm.outputProfile = ColorManagementParams::NoICMString;    
    src->getImage(src->getWB(), tr, img.get(), pp, neutral.exposure, neutral.raw);
    src->convertColorSpace(img.get(), pparams->icm, src->getWB());

    neutral.rotate = pparams->rotate;
    neutral.distortion = pparams->distortion;
    neutral.lensProf = pparams->lensProf;
    ImProcFunctions ipf(&neutral, true);
    if (ipf.needsTransform()) {
        Imagefloat *tmp = new Imagefloat(w, h);
        ipf.transform(img.get(), tmp, 0, 0, 0, 0, w, h, w, h,
                      src->getMetaData(), src->getRotateDegree(), false);
        img.reset(tmp);
    }

    // allocate the gui buffer
    g.buf = static_cast<float *>(malloc(sizeof(float) * w * h * 4));
    g.buf_width = w;
    g.buf_height = h;

    img->normalizeFloatTo1();
    
#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int i = (y * w + x) * 4;
            g.buf[i] = img->r(y, x);
            g.buf[i+1] = img->g(y, x);
            g.buf[i+2] = img->b(y, x);
            g.buf[i+3] = 1.f;
        }
    }

    dt_iop_ashift_fitaxis_t fitaxis = ASHIFT_FIT_NONE;
    switch (dir) {
    case HORIZONTAL:
        fitaxis = ASHIFT_FIT_HORIZONTALLY;
        break;  
    case VERTICAL:
        fitaxis = ASHIFT_FIT_VERTICALLY;
        break;
    default:
        fitaxis = ASHIFT_FIT_BOTH_SHEAR;
        break;
    }

    // reset the pseudo-random seed for repeatability -- ashift_dt uses rand()
    // internally!
    srand(1);
    
    auto res = do_get_structure(&module, &p, ASHIFT_ENHANCE_EDGES) && do_fit(&module, &p, fitaxis);
    procparams::PerspectiveParams retval = pparams->perspective;

    // cleanup the gui
    if (g.lines) free(g.lines);
    if (g.points) free(g.points);
    if (g.points_idx) free(g.points_idx);
    free(g.buf);

    if (res) {
        retval.horizontal = -p.lensshift_h * 100;
        retval.vertical = p.lensshift_v * 100;
        retval.angle = p.rotation;
        retval.shear = p.shear * 100;
    }
    return retval;
}


void PerspectiveCorrection::autocrop(int width, int height, bool fixratio, const procparams::PerspectiveParams &params, const FramesMetaData *metadata, int &x, int &y, int &w, int &h)
{
    auto pp = import_meta(params, metadata);
    double cw, ch;
    get_view_size(width, height, params, cw, ch);
    double s = min(double(width)/cw, double(height)/ch);
    dt_iop_ashift_params_t p;
    dt_iop_ashift_gui_data_t g;
    init_dt_structures(&p, &g, &pp);
    dt_iop_module_t module = { &g, false };
    g.buf_width = width;
    g.buf_height = height;
    p.cropmode = fixratio ? ASHIFT_CROP_ASPECT : ASHIFT_CROP_LARGEST;
    do_crop(&module, &p);
    cw *= s;
    ch *= s;
    double ox = p.cl * cw;
    double oy = p.ct * ch;
    x = ox - (cw - width)/2.0 + 0.5;
    y = oy - (ch - height)/2.0 + 0.5;
    w = (p.cr - p.cl) * cw;
    h = (p.cb - p.ct) * ch;
}

} // namespace rtengine

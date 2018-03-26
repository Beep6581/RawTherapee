/*
 *  This file is part of RawTherapee.
 */
#include "camconst.h"
#include "settings.h"
#include "rt_math.h"
#include <cstdio>
#include <cstring>

// cJSON is a very minimal JSON parser lib in C, not for threaded stuff etc, so if we're going to use JSON more than just
// here we should probably replace cJSON with something beefier.
#include "cJSON.h"
#include <errno.h>
#include <assert.h>
#include <inttypes.h>

namespace rtengine
{

extern const Settings* settings;

CameraConst::CameraConst() : pdafOffset(0)
{
    memset(dcraw_matrix, 0, sizeof(dcraw_matrix));
    memset(raw_crop, 0, sizeof(raw_crop));
    memset(raw_mask, 0, sizeof(raw_mask));
    white_max = 0;
}


bool
CameraConst::parseApertureScaling(CameraConst *cc, void *ji_)
{
    cJSON *ji = (cJSON *)ji_;

    if (ji->type != cJSON_Array) {
        fprintf(stderr, "\"ranges\":\"aperture_scaling\" must be an array\n");
        return false;
    }

    for (ji = ji->child; ji != nullptr; ji = ji->next) {
        cJSON *js = cJSON_GetObjectItem(ji, "aperture");

        if (!js) {
            fprintf(stderr, "missing \"ranges\":\"aperture_scaling\":\"aperture\" object item.\n");
            return false;
        } else if (js->type != cJSON_Number) {
            fprintf(stderr, "\"ranges\":\"aperture_scaling\":\"aperture\" must be a number.\n");
            return false;
        }

        float aperture = (float)js->valuedouble;
        js = cJSON_GetObjectItem(ji, "scale_factor");

        if (!js) {
            fprintf(stderr, "missing \"ranges\":\"aperture_scaling\":\"scale_factor\" object item.\n");
            return false;
        } else if (js->type != cJSON_Number) {
            fprintf(stderr, "\"ranges\":\"aperture_scaling\":\"scale_factor\" must be a number.\n");
            return false;
        }

        float scale_factor = (float)js->valuedouble;
        cc->mApertureScaling.insert(std::pair<float, float>(aperture, scale_factor));
    }

    return true;
}

bool
CameraConst::parseLevels(CameraConst *cc, int bw, void *ji_)
{
    cJSON *ji = (cJSON *)ji_;

    if (ji->type == cJSON_Number) {
        struct camera_const_levels lvl;
        lvl.levels[0] = lvl.levels[1] = lvl.levels[2] = lvl.levels[3] = ji->valueint;
        cc->mLevels[bw].insert(std::pair<int, struct camera_const_levels>(0, lvl));
        return true;
    } else if (ji->type != cJSON_Array) {
        fprintf(stderr, "\"ranges\":\"%s\" must be a number or an array\n", bw ? "white" : "black");
        return false;
    }

    if (ji->child->type == cJSON_Number) {
        struct camera_const_levels lvl;
        int i;
        cJSON *js;

        for (js = ji->child, i = 0; js != nullptr && i < 4; js = js->next, i++) {
            lvl.levels[i] = js->valueint;
        }

        if (i == 3) {
            lvl.levels[3] = lvl.levels[1]; // G2 = G1
        } else if (i == 1) {
            lvl.levels[3] = lvl.levels[2] = lvl.levels[1] = lvl.levels[0];
        } else if (i != 4 || js != nullptr) {
            fprintf(stderr, "\"ranges\":\"%s\" array must have 1, 3 or 4 numbers.\n", bw ? "white" : "black");
            return false;
        }

        cc->mLevels[bw].insert(std::pair<int, struct camera_const_levels>(0, lvl));
        return true;
    }

    for (ji = ji->child; ji != nullptr; ji = ji->next) {
        int iso[1000] = { 0 };
        int iso_count = 0;
        cJSON *js = cJSON_GetObjectItem(ji, "iso");

        if (!js) {
            fprintf(stderr, "missing \"ranges\":\"%s\":\"iso\" object item.\n", bw ? "white" : "black");
            return false;
        } else if (js->type == cJSON_Number) {
            iso[0] = js->valueint;
            iso_count = 1;
        } else if (js->type == cJSON_Array) {
            int i;

            for (js = js->child, i = 0; js != nullptr && i < 1000; js = js->next, i++) {
                if (js->type != cJSON_Number) {
                    fprintf(stderr, "\"ranges\":\"%s\":\"iso\" must be a number or an array of numbers.\n", bw ? "white" : "black");
                    return false;
                }

                iso[i] = js->valueint;
            }

            iso_count = i;
        } else {
            fprintf(stderr, "\"ranges\":\"%s\":\"iso\" must be an array or a number.\n", bw ? "white" : "black");
            return false;
        }

        js = cJSON_GetObjectItem(ji, "levels");

        if (!js) {
            fprintf(stderr, "missing \"ranges\":\"%s\":\"levels\".\n", bw ? "white" : "black");
            return false;
        }

        struct camera_const_levels lvl;

        if (js->type == cJSON_Number) {
            lvl.levels[0] = lvl.levels[1] = lvl.levels[2] = lvl.levels[3] = js->valueint;
        } else if (js->type == cJSON_Array) {
            int i;

            for (js = js->child, i = 0; js != nullptr && i < 4; js = js->next, i++) {
                if (js->type != cJSON_Number) {
                    fprintf(stderr, "\"ranges\":\"%s\":\"levels\" must be a number or an array of numbers.\n", bw ? "white" : "black");
                    return false;
                }

                lvl.levels[i] = js->valueint;
            }

            if (i == 3) {
                lvl.levels[3] = lvl.levels[1]; // G2 = G1
            } else if (i == 1) {
                lvl.levels[3] = lvl.levels[2] = lvl.levels[1] = lvl.levels[0];
            } else if (i != 4 || js != nullptr) {
                fprintf(stderr, "\"ranges\":\"%s\":\"levels\" array must have 1, 3 or 4 numbers.\n", bw ? "white" : "black");
                return false;
            }
        } else {
            fprintf(stderr, "\"ranges\":\"%s\":\"levels\" must be a number or an array of numbers.\n", bw ? "white" : "black");
            return false;
        }

        for (int i = 0; i < iso_count; i++) {
            cc->mLevels[bw].insert(std::pair<int, struct camera_const_levels>(iso[i], lvl));
        }
    }

    return true;
}

CameraConst *
CameraConst::parseEntry(void *cJSON_, const char *make_model)
{
    cJSON *js, *ji, *jranges;
    js = (cJSON *)cJSON_;

    CameraConst *cc = new CameraConst;
    cc->make_model = make_model;

    ji = cJSON_GetObjectItem(js, "dcraw_matrix");

    if (ji) {
        if (ji->type != cJSON_Array) {
            fprintf(stderr, "\"dcraw_matrix\" must be an array\n");
            goto parse_error;
        }

        int i;

        for (i = 0, ji = ji->child; i < 12 && ji != nullptr; i++, ji = ji->next) {
            if (ji->type != cJSON_Number) {
                fprintf(stderr, "\"dcraw_matrix\" array must contain numbers\n");
                goto parse_error;
            }

            cc->dcraw_matrix[i] = (short)ji->valueint;
        }
    }

    ji = cJSON_GetObjectItem(js, "raw_crop");

    if (ji) {
        if (ji->type != cJSON_Array) {
            fprintf(stderr, "\"raw_crop\" must be an array\n");
            goto parse_error;
        }

        int i;

        for (i = 0, ji = ji->child; i < 4 && ji != nullptr; i++, ji = ji->next) {
            if (ji->type != cJSON_Number) {
                fprintf(stderr, "\"raw_crop\" array must contain numbers\n");
                goto parse_error;
            }

            cc->raw_crop[i] = ji->valueint;
        }

        if (i != 4 || ji != nullptr) {
            fprintf(stderr, "\"raw_crop\" must contain 4 numbers\n");
            goto parse_error;
        }
    }

    ji = cJSON_GetObjectItem(js, "masked_areas");

    if (ji) {
        if (ji->type != cJSON_Array) {
            fprintf(stderr, "\"masked_areas\" must be an array\n");
            goto parse_error;
        }

        int i;

        for (i = 0, ji = ji->child; i < 8 * 4 && ji != nullptr; i++, ji = ji->next) {
            if (ji->type != cJSON_Number) {
                fprintf(stderr, "\"masked_areas\" array must contain numbers\n");
                goto parse_error;
            }

            cc->raw_mask[i / 4][i % 4] = ji->valueint;
        }

        if (i % 4 != 0) {
            fprintf(stderr, "\"masked_areas\" array length must be divisable by 4\n");
            goto parse_error;
        }
    }

    jranges = cJSON_GetObjectItem(js, "ranges");

    if (jranges) {
        ji = cJSON_GetObjectItem(jranges, "black");

        if (ji) {
            if (!parseLevels(cc, 0, ji)) {
                goto parse_error;
            }
        }

        ji = cJSON_GetObjectItem(jranges, "white");

        if (ji) {
            if (!parseLevels(cc, 1, ji)) {
                goto parse_error;
            }
        }

        ji = cJSON_GetObjectItem(jranges, "white_max");

        if (ji) {
            if (ji->type != cJSON_Number) {
                fprintf(stderr, "\"ranges\":\"white_max\" must be a number\n");
                goto parse_error;
            }

            cc->white_max = (int)ji->valueint;
        }

        ji = cJSON_GetObjectItem(jranges, "aperture_scaling");

        if (ji) {
            if (!parseApertureScaling(cc, ji)) {
                goto parse_error;
            }
        }
    }

    for (int bw = 0; bw < 2; bw++) {
        struct camera_const_levels lvl;

        if (!cc->get_Levels(lvl, bw, 0, 0)) {
            std::map<int, struct camera_const_levels>::iterator it;
            it = cc->mLevels[bw].begin();

            if (it != cc->mLevels[bw].end()) {
                // insert levels with lowest iso as the default (iso 0)
                struct camera_const_levels lvl = it->second;
                cc->mLevels[bw].insert(std::pair<int, struct camera_const_levels>(0, lvl));
            }
        }
    }

    ji = cJSON_GetObjectItem(js, "pdafPattern");

    if (ji) {
        if (ji->type != cJSON_Array) {
            fprintf(stderr, "\"pdafPattern\" must be an array\n");
            goto parse_error;
        }

        for (ji = ji->child; ji != nullptr; ji = ji->next) {
            if (ji->type != cJSON_Number) {
                fprintf(stderr, "\"pdafPattern\" array must contain numbers\n");
                goto parse_error;
            }

            cc->pdafPattern.push_back(ji->valueint);
        }
    }

    ji = cJSON_GetObjectItem(js, "pdafOffset");

    if (ji) {
        if (ji->type != cJSON_Number) {
            fprintf(stderr, "\"pdafOffset\" must contain a number\n");
            goto parse_error;
        }

        cc->pdafOffset = ji->valueint;
    }

    return cc;

parse_error:
    delete cc;
    return nullptr;
}

bool
CameraConst::has_dcrawMatrix()
{
    return dcraw_matrix[0] != 0;
}

void
CameraConst::update_dcrawMatrix(const short *other)
{
    if (!other) {
        return;
    }

    for (int i = 0; i < 12; ++i) {
        dcraw_matrix[i] = other[i];
    }
}

const short *
CameraConst::get_dcrawMatrix()
{
    if (!has_dcrawMatrix()) {
        return nullptr;
    }

    return dcraw_matrix;
}

bool
CameraConst::has_pdafPattern()
{
    return pdafPattern.size() > 0;
}

std::vector<int>
CameraConst::get_pdafPattern()
{
    return pdafPattern;
}

void
CameraConst::update_pdafPattern(const std::vector<int> &other)
{
    if (other.empty()) {
        return;
    }
    pdafPattern = other;
}

void
CameraConst::update_pdafOffset(int other)
{
    if (other == 0) {
        return;
    }
    pdafOffset = other;
}

bool
CameraConst::has_rawCrop()
{
    return raw_crop[0] != 0 || raw_crop[1] != 0 || raw_crop[2] != 0 || raw_crop[3] != 0;
}

void
CameraConst::get_rawCrop(int& left_margin, int& top_margin, int& width, int& height)
{
    left_margin = raw_crop[0];
    top_margin = raw_crop[1];
    width = raw_crop[2];
    height = raw_crop[3];
}

bool
CameraConst::has_rawMask(int idx)
{
    if (idx < 0 || idx > 7) {
        return false;
    }

    return (raw_mask[idx][0] | raw_mask[idx][1] | raw_mask[idx][2] | raw_mask[idx][3]) != 0;
}

void
CameraConst::get_rawMask(int idx, int& top, int& left, int& bottom, int& right)
{
    top = left = bottom = right = 0;

    if (idx < 0 || idx > 7) {
        return;
    }

    top =    raw_mask[idx][0];
    left =   raw_mask[idx][1];
    bottom = raw_mask[idx][2];
    right =  raw_mask[idx][3];
}

void
CameraConst::update_Levels(const CameraConst *other)
{
    if (!other) {
        return;
    }

    if (other->mLevels[0].size()) {
        mLevels[0].clear();
        mLevels[0] = other->mLevels[0];
    }

    if (other->mLevels[1].size()) {
        mLevels[1].clear();
        mLevels[1] = other->mLevels[1];
    }

    if (other->mApertureScaling.size()) {
        mApertureScaling.clear();
        mApertureScaling = other->mApertureScaling;
    }

    if (other->white_max) {
        white_max = other->white_max;
    }

//  for (std::map<int, struct camera_const_levels>::iterator i=other->mLevels[0].begin(); i!=other->mLevels[0].end(); i++) {
//  }
}

void
CameraConst::update_Crop(CameraConst *other)
{
    if (!other) {
        return;
    }

    if (other->has_rawCrop()) {
        other->get_rawCrop(raw_crop[0], raw_crop[1], raw_crop[2], raw_crop[3]);
    }
}

bool
CameraConst::get_Levels(struct camera_const_levels & lvl, int bw, int iso, float fnumber)
{
    std::map<int, struct camera_const_levels>::iterator it;
    it = mLevels[bw].find(iso);

    if (it == mLevels[bw].end()) {
        std::map<int, struct camera_const_levels>::iterator best_it = mLevels[bw].begin();

        if (iso > 0) {
            for (it = mLevels[bw].begin(); it != mLevels[bw].end(); ++it) {
                if (abs(it->first - iso) <= abs(best_it->first - iso)) {
                    best_it = it;
                } else {
                    break;
                }
            }
        }

        it = best_it;

        if (it == mLevels[bw].end()) {
            return false;
        }
    }

    lvl = it->second;

    if (bw == 1 && fnumber > 0 && mApertureScaling.size() > 0) {
        std::map<float, float>::iterator it;
        it = mApertureScaling.find(fnumber);

        if (it == mApertureScaling.end()) {
            // fnumber may be an exact aperture, eg 1.414, or a rounded eg 1.4. In our map we
            // should have rounded numbers so we translate and retry the lookup

            // table with traditional 1/3 stop f-number rounding used by most cameras, we only
            // have in the range 0.7 - 10.0, but aperture scaling rarely happen past f/4.0
            const float fn_tab[8][3] = {
                { 0.7, 0.8, 0.9 },
                { 1.0, 1.1, 1.2 },
                { 1.4, 1.6, 1.8 },
                { 2.0, 2.2, 2.5 },
                { 2.8, 3.2, 3.5 },
                { 4.0, 4.5, 5.0 },
                { 5.6, 6.3, 7.1 },
                { 8.0, 9.0, 10.0 }
            };

            for (int avh = 0; avh < 8; avh++) {
                for (int k = 0; k < 3; k++) {
                    float av = (avh - 1) + (float)k / 3;
                    float aperture = sqrtf(powf(2, av));

                    if (fnumber > aperture * 0.97 && fnumber < aperture / 0.97) {
                        fnumber = fn_tab[avh][k];
                        it = mApertureScaling.find(fnumber);
                        avh = 7;
                        break;
                    }
                }
            }
        }

        float scaling = 1.0;

        if (it == mApertureScaling.end()) {
            std::map<float, float>::reverse_iterator it;

            for (it = mApertureScaling.rbegin(); it != mApertureScaling.rend(); ++it) {
                if (it->first > fnumber) {
                    scaling = it->second;
                } else {
                    break;
                }
            }
        } else {
            scaling = it->second;
        }

        if (scaling > 1.0) {
            for (int i = 0; i < 4; i++) {
                lvl.levels[i] *= scaling;

                if (white_max > 0 && lvl.levels[i] > white_max) {
                    lvl.levels[i] = white_max;
                }
            }
        }
    }

    return true;
}

int
CameraConst::get_BlackLevel(const int idx, const int iso_speed)
{
    assert(idx >= 0 && idx <= 3);
    struct camera_const_levels lvl;

    if (!get_Levels(lvl, 0, iso_speed, 0.0)) {
        return -1;
    }

    return lvl.levels[idx];
}

int
CameraConst::get_WhiteLevel(const int idx, const int iso_speed, const float fnumber)
{
    assert(idx >= 0 && idx <= 3);
    struct camera_const_levels lvl;

    if (!get_Levels(lvl, 1, iso_speed, fnumber)) {
        return -1;
    }

    return lvl.levels[idx];
}

bool
CameraConstantsStore::parse_camera_constants_file(Glib::ustring filename_)
{
    // read the file into a single long string
    const char *filename = filename_.c_str();
    FILE *stream = fopen(filename, "rt");

    if (stream == nullptr) {
        fprintf(stderr, "Could not open camera constants file \"%s\": %s\n", filename, strerror(errno));
        return false;
    }

    size_t bufsize = 16384;
    size_t increment = 2 * bufsize;
    size_t datasize = 0, ret;
    char *buf = (char *)malloc(bufsize);

    while ((ret = fread(&buf[datasize], 1, bufsize - datasize, stream)) != 0) {
        datasize += ret;

        if (datasize == bufsize) { // we need more memory
            bufsize += increment;
            void *temp = realloc(buf, bufsize); // try to realloc buffer with new size
            if(!temp) { // realloc failed
                temp = malloc(bufsize); // alloc now buffer
                if (temp) { // alloc worked
                    memcpy(temp, buf, bufsize - increment); // copy old buffer content to new buffer
                    free(buf); // free old buffer
                } else { // alloc didn't work, break
                    break;
                }
            }
            buf = (char *)temp; // assign new buffer
            increment *= 2; // double increment
        }
    }

    if (!feof(stream)) {
        fclose(stream);
        free(buf);
        fprintf(stderr, "Failed to read camera constants file \"%s\"\n", filename);
        return false;
    }

    fclose(stream);

    if(datasize == bufsize) {
        buf = (char *)realloc(buf, datasize + 1);
    }

    buf[datasize] = '\0';

    // remove comments
    cJSON_Minify(buf);

    // parse
    cJSON *jsroot = cJSON_Parse(buf);

    if (!jsroot) {
        char str[128];
        const char *ep = cJSON_GetErrorPtr() - 10;

        if ((uintptr_t)ep < (uintptr_t)buf) {
            ep = buf;
        }

        strncpy(str, ep, sizeof(str));
        str[sizeof(str) - 1] = '\0';
        fprintf(stderr, "JSON parse error in file \"%s\" near '%s'\n", filename, str);
        free(buf);
        return false;
    }

    free(buf);
    /*{
        char *js_str = cJSON_Print(jsroot);
        printf("%s\n", js_str);
        free(js_str);
    }*/
    cJSON *js = cJSON_GetObjectItem(jsroot, "camera_constants");

    if (!js) {
        fprintf(stderr, "missing \"camera_constants\" object item\n");
        goto parse_error;
    }

    for (js = js->child; js != nullptr; js = js->next) {
        cJSON *ji = cJSON_GetObjectItem(js, "make_model");

        if (!ji) {
            fprintf(stderr, "missing \"make_model\" object item\n");
            goto parse_error;
        }

        bool is_array = false;

        if (ji->type == cJSON_Array) {
            ji = ji->child;
            is_array = true;
        }

        while (ji != nullptr) {
            if (ji->type != cJSON_String) {
                fprintf(stderr, "\"make_model\" must be a string or an array of strings\n");
                goto parse_error;
            }

            CameraConst *cc = CameraConst::parseEntry((void *)js, ji->valuestring);

            if (!cc) {
                goto parse_error;
            }

            Glib::ustring make_model(ji->valuestring);
            make_model = make_model.uppercase();

            const auto ret = mCameraConstants.emplace(make_model, cc);

            if(ret.second) { // entry inserted into map
                if (settings->verbose) {
                    printf("Add camera constants for \"%s\"\n", make_model.c_str());
                }
            } else {
                // The CameraConst already exist for this camera make/model -> we merge the values
                CameraConst *existingcc = ret.first->second;

                // updating the dcraw matrix
                existingcc->update_dcrawMatrix(cc->get_dcrawMatrix());
                // deleting all the existing levels, replaced by the new ones
                existingcc->update_Levels(cc);
                existingcc->update_Crop(cc);
                existingcc->update_pdafPattern(cc->get_pdafPattern());
                existingcc->update_pdafOffset(cc->get_pdafOffset());

                if (settings->verbose) {
                    printf("Merging camera constants for \"%s\"\n", make_model.c_str());
                }
            }

            if (is_array) {
                ji = ji->next;
            } else {
                ji = nullptr;
            }
        }
    }

    cJSON_Delete(jsroot);
    return true;

parse_error:
    fprintf(stderr, "failed to parse camera constants file \"%s\"\n", filename);
    mCameraConstants.clear();
    cJSON_Delete(jsroot);
    return false;
}

CameraConstantsStore::CameraConstantsStore()
{
}


CameraConstantsStore::~CameraConstantsStore()
{
    for (auto &p : mCameraConstants) {
        delete p.second;
    }
}

void CameraConstantsStore::init(Glib::ustring baseDir, Glib::ustring userSettingsDir)
{
    parse_camera_constants_file(Glib::build_filename(baseDir, "camconst.json"));

    Glib::ustring userFile(Glib::build_filename(userSettingsDir, "camconst.json"));

    if (Glib::file_test(userFile, Glib::FILE_TEST_EXISTS)) {
        parse_camera_constants_file(userFile);
    }
}

CameraConstantsStore *
CameraConstantsStore::getInstance()
{
    static CameraConstantsStore instance_;
    return &instance_;
}

CameraConst *
CameraConstantsStore::get(const char make[], const char model[])
{
    Glib::ustring key(make);
    key += " ";
    key += model;
    key = key.uppercase();
    std::map<std::string, CameraConst *>::iterator it;
    it = mCameraConstants.find(key);

    if (it == mCameraConstants.end()) {
        return nullptr;
    }

    return it->second;
}

} // namespace rtengine

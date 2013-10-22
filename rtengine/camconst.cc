/*
 *  This file is part of RawTherapee.
 */
#include "camconst.h"
#include "safegtk.h"
#include "rt_math.h"
#include <cstdio>
#include <cstring>

// cJSON is a very minimal JSON parser lib in C, not for threaded stuff etc, so if we're going to use JSON more than just
// here we should probably replace cJSON with something beefier.
#include "cJSON.h"
#include <errno.h>
#include <assert.h>
#include <inttypes.h>

namespace rtengine {

CameraConst::CameraConst()
{
	memset(dcraw_matrix, 0, sizeof(dcraw_matrix));
	white_max = 0;
}

CameraConst::~CameraConst()
{
}

bool
CameraConst::parseApertureScaling(CameraConst *cc, void *ji_)
{
	cJSON *ji = (cJSON *)ji_;
	if (ji->type != cJSON_Array) {
		fprintf(stderr, "\"ranges\":\"aperture_scaling\" must be an array\n");
		return false;
	}
	for (ji = ji->child; ji != NULL; ji = ji->next) {
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
		cc->mApertureScaling.insert(std::pair<float,float>(aperture, scale_factor));
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
		cc->mLevels[bw].insert(std::pair<int,struct camera_const_levels>(0, lvl));
		return true;
	} else if (ji->type != cJSON_Array) {
		fprintf(stderr, "\"ranges\":\"%s\" must be a number or an array\n", bw ? "white" : "black");
		return false;
	}

	if (ji->child->type == cJSON_Number) {
		struct camera_const_levels lvl;
		int i;
		cJSON *js;
		for (js = ji->child, i = 0; js != NULL && i < 4; js = js->next, i++) {
			lvl.levels[i] = js->valueint;
		}
		if (i == 3) {
			lvl.levels[3] = lvl.levels[1]; // G2 = G1
		} else if (i == 1) {
			lvl.levels[3] = lvl.levels[2] = lvl.levels[1] = lvl.levels[0];
		} else if (i != 4 || js != NULL) {
			fprintf(stderr, "\"ranges\":\"%s\" array must have 1, 3 or 4 numbers.\n", bw ? "white" : "black");
			return false;
		}
		cc->mLevels[bw].insert(std::pair<int,struct camera_const_levels>(0, lvl));
		return true;
	}

	for (ji = ji->child; ji != NULL; ji = ji->next) {
		int iso = 0;
		cJSON *js = cJSON_GetObjectItem(ji, "iso");
		if (!js) {
			fprintf(stderr, "missing \"ranges\":\"%s\":\"iso\" object item.\n", bw ? "white" : "black");
			return false;
		} else if (js->type != cJSON_Number) {
			fprintf(stderr, "\"ranges\":\"%s\":\"iso\" must be a a number.\n", bw ? "white" : "black");
			return false;
		}
		iso = js->valueint;
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
			for (js = js->child, i = 0; js != NULL && i < 4; js = js->next, i++) {
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
			} else if (i != 4 || js != NULL) {
				fprintf(stderr, "\"ranges\":\"%s\":\"levels\" array must have 1, 3 or 4 numbers.\n", bw ? "white" : "black");
				return false;
			}
		} else {
			fprintf(stderr, "\"ranges\":\"%s\":\"levels\" must be a number or an array of numbers.\n", bw ? "white" : "black");
			return false;
		}
		cc->mLevels[bw].insert(std::pair<int,struct camera_const_levels>(iso, lvl));
	}
	return true;
}

CameraConst *
CameraConst::parseEntry(void *cJSON_, const char *make_model)
{
	CameraConst *cc = 0;
	cJSON *js, *ji, *jranges;
	js = (cJSON *)cJSON_;

	cc = new CameraConst;
	cc->make_model = Glib::ustring(make_model);

	ji = cJSON_GetObjectItem(js, "dcraw_matrix");
	if (ji) {
		if (ji->type != cJSON_Array) {
			fprintf(stderr, "\"dcraw_matrix\" must be an array\n");
			goto parse_error;
		}
		int i;
		for (i = 0, ji = ji->child; i < 12 && ji != NULL; i++, ji = ji->next) {
			if (ji->type != cJSON_Number) {
				fprintf(stderr, "\"dcraw_matrix\" array must contain numbers\n");
				goto parse_error;
			}
			cc->dcraw_matrix[i] = (short)ji->valueint;
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
				cc->mLevels[bw].insert(std::pair<int,struct camera_const_levels>(0, lvl));
			}
		}
	}
	return cc;

parse_error:
	return 0;
}

bool
CameraConst::has_dcrawMatrix(void)
{
	return dcraw_matrix[0] != 0;
}

void
CameraConst::update_dcrawMatrix(const short *other) {
	if (!other)
		return;

	for (int i=0; i<12; ++i)
		dcraw_matrix[i] = other[i];
}

const short *
CameraConst::get_dcrawMatrix(void)
{
	if (!has_dcrawMatrix()) {
		return 0;
	}
	return dcraw_matrix;
}

void
CameraConst::update_Levels(const CameraConst *other) {
	if (!other)
		return;

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
	if (other->white_max)
		white_max = other->white_max;

//	for (std::map<int, struct camera_const_levels>::iterator i=other->mLevels[0].begin(); i!=other->mLevels[0].end(); i++) {
//	}
}

bool
CameraConst::get_Levels(struct camera_const_levels & lvl, int bw, int iso, float fnumber)
{
	std::map<int, struct camera_const_levels>::iterator it;
	it = mLevels[bw].find(iso);
	if (it == mLevels[bw].end()) {
		std::map<int, struct camera_const_levels>::iterator best_it = mLevels[bw].begin();
		if (iso > 0) {
			for (it = mLevels[bw].begin(); it != mLevels[bw].end(); it++) {
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
					float av = (avh-1) + (float)k / 3;
					float aperture = sqrtf(powf(2, av));
					if (fnumber > aperture*0.97 && fnumber < aperture/0.97) {
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
			for (it = mApertureScaling.rbegin(); it != mApertureScaling.rend(); it++) {
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
	if (stream == NULL) {
		fprintf(stderr, "Could not open camera constants file \"%s\": %s\n", filename, strerror(errno));
		return false;
	}
	size_t bufsize = 4096;
	size_t datasize = 0, ret;
	char *buf = (char *)malloc(bufsize);
	while ((ret = fread(&buf[datasize], 1, bufsize - datasize, stream)) != 0) {
		datasize += ret;
		if (datasize == bufsize) {
			bufsize += 4096;
			buf = (char *)realloc(buf, bufsize);
		}
	}
	if (!feof(stream)) {
		fclose(stream);
		free(buf);
		fprintf(stderr, "Failed to read camera constants file \"%s\"\n", filename);
		return false;
	}
	fclose(stream);
	buf = (char *)realloc(buf, datasize + 1);
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
		str[sizeof(str)-1] = '\0';
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
	for (js = js->child; js != NULL; js = js->next) {
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
		while (ji != NULL) {
			if (ji->type != cJSON_String) {
				fprintf(stderr, "\"make_model\" must be a string or an array of strings\n");
				goto parse_error;
			}
			CameraConst *cc = CameraConst::parseEntry((void *)js, ji->valuestring);
			if (!cc) {
				goto parse_error;
			}
			Glib::ustring make_model(ji->valuestring);
			std::map<Glib::ustring, CameraConst *>::iterator existingccIter = mCameraConstants.find(make_model);

			if (existingccIter == mCameraConstants.end())
				// add the new CamConst to the map
				mCameraConstants.insert(std::pair<Glib::ustring,CameraConst *>(make_model, cc));
			else {
				// The CameraConst already exist for this camera make/model -> we merge the values
				CameraConst *existingcc = existingccIter->second;

				// updating the dcraw matrix
				existingcc->update_dcrawMatrix(cc->get_dcrawMatrix());
				// deleting all the existing levels, replaced by the new ones
				existingcc->update_Levels(cc);
			}
			if (is_array) {
				ji = ji->next;
			} else {
				ji = NULL;
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

static CameraConstantsStore *global_instance;

void CameraConstantsStore::initCameraConstants(Glib::ustring baseDir, Glib::ustring userSettingsDir)
{
	if (global_instance) {
		// should only be called once during init.
		abort();
	}
	global_instance = new CameraConstantsStore();
	global_instance->parse_camera_constants_file(Glib::build_filename(baseDir, "camconst.json"));

	Glib::ustring userFile(Glib::build_filename(userSettingsDir, "camconst.json"));
	if (safe_file_test(userFile, Glib::FILE_TEST_EXISTS))
		global_instance->parse_camera_constants_file(userFile);
}

CameraConstantsStore *
CameraConstantsStore::getInstance(void)
{
	return global_instance;
}

CameraConst *
CameraConstantsStore::get(const char make[], const char model[])
{
	Glib::ustring key(make);
	key += " ";
	key += model;
	std::map<Glib::ustring, CameraConst *>::iterator it;
	it = mCameraConstants.find(key);
	if (it == mCameraConstants.end()) {
		return 0;
	}
	return it->second;
}

} // namespace rtengine

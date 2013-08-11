/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2011 Jean-Christophe FRISCH <natureh@free.fr>
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

#include "rtimage.h"
#include "../rtengine/safegtk.h"
#include "../rtengine/safekeyfile.h"

extern Glib::ustring argv0;
extern Options options;

std::vector<Glib::ustring> imagesPaths;
std::vector<RTImage*>       imagesList;			// List of images in order to live update them on theme switch

/*
 * RTImage is a derived class of Gtk::Image, in order to handle theme related iconsets
 */
RTImage::RTImage(Glib::ustring fileName, Glib::ustring rtlFileName) : Gtk::Image() {
	Glib::ustring path;
	if (rtlFileName.length()) {
		const Gtk::TextDirection dir = get_direction();
		if (dir == Gtk::TEXT_DIR_RTL)
			path = findIconAbsolutePath(rtlFileName);
		else
			path = findIconAbsolutePath(fileName);
	}
	else
		path = findIconAbsolutePath(fileName);

	set(path);
	imagesList.push_back(this);
}

RTImage::~RTImage() {
	// Remove the image from the global images list
	std::vector<RTImage*>::iterator i = std::find (imagesList.begin(), imagesList.end(), this);
	if (i!=imagesList.end())
		imagesList.erase(i);
}

void RTImage::updateImages() {
	for (unsigned int i=0; i<imagesList.size(); i++) {
		Glib::ustring oldPath = imagesList[i]->property_file();
		Glib::ustring fileName = Glib::path_get_basename(oldPath);
		imagesList[i]->clear();
		Glib::ustring fullPath = findIconAbsolutePath(fileName);
		imagesList[i]->set(fullPath);
	}
}

// TODO: Maybe this could be optimized: in order to avoid looking up for an icon file in the filesystem on each popupmenu selection, maybe we could find a way to copy the image data from another RTImage
void RTImage::changeImage(Glib::ustring &newImage) {
	clear();
	Glib::ustring fullPath = findIconAbsolutePath(newImage);
	set(fullPath);
}

Glib::ustring RTImage::findIconAbsolutePath(const Glib::ustring &iconFName) {
	Glib::ustring path;
	for (unsigned int i=0; i<imagesPaths.size(); i++) {
		path = Glib::build_filename(imagesPaths[i], iconFName);
		//printf("\"%s\" \n", path.c_str());
		if (safe_file_test(path, Glib::FILE_TEST_EXISTS)) {
			//printf("Found!\n");
			return path;
		}
	}
	printf("\"%s\" not found!\n", iconFName.c_str());
	return "";
}

void RTImage::setPaths(Options &opt) {
	Glib::ustring configFilename;
	rtengine::SafeKeyFile keyFile;
	bool hasKeyFile = true;

	imagesPaths.clear();

	// system theme will use the theme set in system.iconset
	if (opt.useSystemTheme) {
		configFilename = Glib::build_filename(argv0, Glib::build_filename("themes","system.iconset"));
	}
	// Gtk theme will use the theme set in it's *.iconset fiel, if it exists
	else {
		configFilename = Glib::build_filename(argv0, Glib::build_filename("themes", Glib::ustring::format(opt.theme, ".iconset")));
	}
	try {
		if (!safe_file_test(configFilename, Glib::FILE_TEST_EXISTS) || !keyFile.load_from_file (configFilename)) {
			// ...otherwise fallback to the iconset set in default.iconset
			configFilename = Glib::build_filename(argv0, Glib::build_filename("themes", "Default.iconset"));
			if (!keyFile.load_from_file (configFilename)) {
				hasKeyFile = false;
			}
		}
	}
	catch (Glib::Error &err) {
		if (options.rtSettings.verbose)
			printf("RTImage::setPaths / Error code %d while reading values from \"%s\":\n%s\n", err.code(), configFilename.c_str(), err.what().c_str());
	}
	catch (...) {
		if (options.rtSettings.verbose)
			printf("RTImage::setPaths / Unknown exception while trying to load \"%s\"!\n", configFilename.c_str());
	}

	if (hasKeyFile && keyFile.has_group ("General")) {
		Glib::ustring iSet;

		if (keyFile.has_key ("General", "Iconset"))
			iSet = keyFile.get_string ("General", "Iconset");
		if (iSet.length()) {
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "actions"))));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", iSet)));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "devices"))));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "places"))));
		}

		iSet.clear();
		if (keyFile.has_key ("General", "FallbackIconset"))
			iSet = keyFile.get_string ("General", "FallbackIconset");
		if (iSet.length()) {
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "actions"))));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", iSet)));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "devices"))));
			imagesPaths.push_back (Glib::build_filename(argv0, Glib::build_filename("images", Glib::build_filename(iSet, "places"))));
		}
	}
	// The images/ folder is the second fallback solution
	imagesPaths.push_back (Glib::build_filename(argv0, "images"));
}


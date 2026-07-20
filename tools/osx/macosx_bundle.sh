#!/usr/bin/env bash

# Required variables, must be set in rtdata/CMakeLists.txt
# - PROJECT_NAME
# - PROJECT_VERSION (if without git)
# - PROJECT_SOURCE_DIR
# - CMAKE_BUILD_TYPE
# - PROC_BIT_DEPTH
# - GTK_PREFIX

# Formatting
fNormal="$(printf "\e[0m")"
fBold="$(printf "\e[1m")"
# Colors depend upon the user's terminal emulator color scheme - what is readable for you may be not readable for someone else.
fMagenta="$(printf "\e[1;35m")"
fRed="$(printf "\e[1;31m")"

function msg {
    printf "\\n${fBold}-- %s${fNormal}\\n" "${@}"
}

function msgError {
    printf "\\n${fBold}Error:${fNormal}\\n%s\\n" "${@}"
}

function GetDependencies {
    otool -L "$1" | awk 'NR >= 2 && $1 !~ /^(\/usr\/lib|\/System|@executable_path|@loader_path|@rpath)\// { print $1 }'  2>&1
}

function CheckLink {
    GetDependencies "$1" | while read -r; do
        local dest="${LIB}/$(basename "${REPLY}")"
        test -f "${dest}" || { ditto --arch "${arch}" "${REPLY}" "${dest}"; CheckLink "${dest}"; }
    done
}

function ModifyInstallNames {
    find -E "${CONTENTS}" -type f -regex '.*/(rawtherapee-cli|rawtherapee|.*\.(dylib|so))' | while read -r x; do
        msg "Modifying install names: ${x}"
        {
            # id
            if [[ ${x:(-6)} == ".dylib" ]] || [[ ${x:(-3)} == ".so" ]]; then
                install_name_tool -id /Applications/"${LIB}"/"$(basename "${x}")" "${x}" 2>/dev/null
            fi
            GetDependencies "${x}" | while read -r y
            do
                install_name_tool -change ${y} /Applications/"${LIB}"/$(basename ${y}) ${x} 2>/dev/null
            done
        } | bash -v
    done
}

die() {
    msgError "$*"
    exit 1
}

copy_tree() {
    local src=$1
    local dst=$2

    [[ -d "$src" ]] || die "Missing directory: $src"
    install -d "$dst"
    ditto "$src" "$dst" ||
        die "Failed to copy directory: $src -> $dst"
}

# Copy a resource hierarchy as real files. Package-manager share directories
# commonly contain relative symlinks into a Cellar or ports tree; preserving
# those links makes the application bundle non-relocatable and invalidates it
# for codesigning. The /. suffix gives cp stable merge semantics when dst
# already exists.
copy_tree_dereference() {
    local src=$1
    local dst=$2

    [[ -d "$src" ]] || die "Missing directory: $src"
    install -d "$dst"
    cp -RL "${src}/." "${dst}/" ||
        die "Failed to dereference directory: $src -> $dst"
}

copy_macho() {
    local src=$1
    local dst=$2

    # -e follows a command-line symlink and verifies its target exists.
    [[ -e "$src" ]] || die "Missing library: $src"
    install -d "$(dirname "$dst")"
    ditto --arch "$arch" "$src" "$dst" ||
        die "Failed to copy architecture $arch: $src -> $dst"
    [[ -f "$dst" ]] || die "Copy did not create: $dst"
}

# Store the first existing candidate in RESOLVED_FILE. This avoids assuming a
# package-manager prefix (/opt/local, /usr/local, or /opt/homebrew).
resolve_file() {
    local description=$1
    shift
    local candidate

    RESOLVED_FILE=""
    for candidate in "$@"; do
        if [[ -n "$candidate" && -e "$candidate" ]]; then
            RESOLVED_FILE="$candidate"
            return
        fi
    done

    die "Could not locate ${description}"
}

resolve_directory() {
    local description=$1
    shift
    local candidate

    RESOLVED_DIRECTORY=""
    for candidate in "$@"; do
        if [[ -n "$candidate" && -d "$candidate" ]]; then
            RESOLVED_DIRECTORY="$candidate"
            return
        fi
    done

    die "Could not locate ${description}"
}

# Source check
if [[ ! -d $CMAKE_BUILD_TYPE ]]; then
    msgError "${PWD}/${CMAKE_BUILD_TYPE} folder does not exist. Please execute 'make install' first."
    exit 1
fi

# Update project version
if [[ -x $(which git) && -d $PROJECT_SOURCE_DIR/.git ]]; then
    ### This section is copied from tools/generateReleaseInfo
    # Get version description.
    # Depending on whether you checked out a branch (dev) or a tag (release),
    # "git describe" will return "5.0-gtk2-2-g12345678" or "5.0-gtk2", respectively.
    gitDescribe="$(git describe --tags --always)"

    # Apple requires a numeric version of the form n.n.n
    # https://goo.gl/eWDQv6

    # Get number of commits since tagging. This is what gitDescribe uses.
    # Works when checking out branch, tag or commit.
    gitCommitsSinceTag="$(git rev-list --count HEAD --not $(git tag --merged HEAD))"

    # Create numeric version.
    # This version is nonsense, either don't use it at all or use it only where you have no other choice, e.g. Inno Setup's VersionInfoVersion.
    # Strip everything after hyphen, e.g. "5.0-gtk2" -> "5.0", "5.1-rc1" -> "5.1" (ergo BS).
    if [[ -z $gitCommitsSinceTag ]]; then
        gitVersionNumericBS="0.0.0"
    else
        gitVersionNumericBS="${gitDescribe%%-*}" # Remove everything after first hyphen.
        gitVersionNumericBS="${gitVersionNumericBS}.${gitCommitsSinceTag}" # Remove everything until after first hyphen: 5.0
    fi
    ### Copy end.

    PROJECT_FULL_VERSION="$gitDescribe"
    PROJECT_VERSION="$gitVersionNumericBS"
fi

#In: CMAKE_OSX_DEPLOYMENT_TARGET=11.3
#Out: 11.3
CMAKE_OSX_DEPLOYMENT_TARGET="$(cmake .. -L -N | grep CMAKE_OSX_DEPLOYMENT_TARGET)"; CMAKE_OSX_DEPLOYMENT_TARGET="${CMAKE_OSX_DEPLOYMENT_TARGET#*=}"

#In: CMAKE_OSX_ARCHITECTURES=x86_64
#Out: x86_64
CMAKE_OSX_ARCHITECTURES="$(cmake .. -L -N | grep CMAKE_OSX_ARCHITECTURES)"; CMAKE_OSX_ARCHITECTURES="${CMAKE_OSX_ARCHITECTURES#*=}"

MINIMUM_SYSTEM_VERSION="$(otool -l "${CMAKE_BUILD_TYPE}"/MacOS/rawtherapee | grep -A2 'LC_VERSION_MIN_MACOSX' | awk '$1 ~ /version/ { printf $2 }')"
if [[ -z $MINIMUM_SYSTEM_VERSION ]]; then
    MINIMUM_SYSTEM_VERSION=${CMAKE_OSX_DEPLOYMENT_TARGET}
fi

arch=${CMAKE_OSX_ARCHITECTURES}

cat <<__EOS__
PROJECT_NAME:           ${PROJECT_NAME}
PROJECT_VERSION:        ${PROJECT_VERSION}
PROJECT_SOURCE_DIR:     ${PROJECT_SOURCE_DIR}
CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}
PROC_BIT_DEPTH:         ${PROC_BIT_DEPTH}
MINIMUM_SYSTEM_VERSION: ${MINIMUM_SYSTEM_VERSION}
GTK_PREFIX:             ${GTK_PREFIX}
PWD:                    ${PWD}
__EOS__

minimum_macos_version=${MINIMUM_SYSTEM_VERSION}

# Retrieve cached values from cmake

#In: LOCAL_PREFIX:STRING=/opt
#Out: /opt
LOCAL_PREFIX="$(cmake .. -L -N | grep LOCAL_PREFIX)"; LOCAL_PREFIX="${LOCAL_PREFIX#*=}"

#In: OSX_UNIVERSAL_URL=file:/// etc.
#Out: file:/// etc.
UNIVERSAL_URL="$(cmake .. -L -N | grep OSX_UNIVERSAL_URL)"; UNIVERSAL_URL="${UNIVERSAL_URL#*=}"
if [[ -n $UNIVERSAL_URL ]]; then
    echo "Universal app is ON. The URL is ${UNIVERSAL_URL}"
fi

#In: pkgcfg_lib_EXPAT_expat:FILEPATH=/opt/local/lib/libexpat.dylib
#Out: /opt/local/lib/libexpat.dylib
EXPATLIB="$(cmake .. -LA -N | grep pkgcfg_lib_EXPAT_expat)"; EXPATLIB="${EXPATLIB#*=}"

#In: CODESIGNID:STRING=Developer ID Application: Doctor Who (1234567890)
#Out: Developer ID Application: Doctor Who (1234567890)
CODESIGNID="$(cmake .. -L -N | grep CODESIGNID)"; CODESIGNID="${CODESIGNID#*=}"

#In: NOTARY:STRING="--apple-id drwho@bbc.com --password abcd-efgh-hijk-lmno --team-id ABCDE12345"
#Out: --apple-id drwho@bbc.com --password abcd-efgh-hijk-lmno --team-id ABCDE12345
NOTARY="$(cmake .. -L -N | grep NOTARY)"; NOTARY="${NOTARY#*=}"

# In: FANCY_DMG:BOOL=ON
# Out: ON
FANCY_DMG="$(cmake .. -L -N | grep FANCY_DMG)"; FANCY_DMG="${FANCY_DMG#*=}"
if [[ -n $FANCY_DMG ]]; then
    echo "Fancy .dmg build is ON."
fi

# In: OSX_NIGHTLY:BOOL=ON
# Out: ON
OSX_NIGHTLY="$(cmake .. -L -N | grep OSX_NIGHTLY)"; NIGHTLY="${OSX_NIGHTLY#*=}"
if [[ -n $NIGHTLY ]]; then
    echo "Nightly/generically-named zip is ON."
fi

# In: OSX_CONTINUOUS:BOOL=ON
# Out: ON
OSX_CONTINUOUS="$(cmake .. -L -N | grep OSX_CONTINUOUS)"; CONTINUOUS="${OSX_CONTINUOUS#*=}"
if [[ -n $CONTINUOUS ]]; then
    NIGHTLY="${OSX_CONTINUOUS#*=}"
    echo "Continuous/generically-named zip is ON."
fi

APP="${PROJECT_NAME}.app"
CONTENTS="${APP}/Contents"
RESOURCES="${CONTENTS}/Resources"
MACOS="${CONTENTS}/MacOS"
LIB="${CONTENTS}/Frameworks"
ETC="${RESOURCES}/etc"
EXECUTABLE="${MACOS}/rawtherapee"

msg "Removing old files:"
rm -rf "${APP}" *.dmg *.zip *.app

msg "Creating bundle container:"
install -d "${RESOURCES}"
install -d "${MACOS}"
install -d "${LIB}"
install -d "${ETC}"

msg "Copying binary executable files."
ditto "${CMAKE_BUILD_TYPE}/MacOS" "${MACOS}"

msg "Copying Resources directory."
#cp AboutThisBuild.txt "${RESOURCES}"
ditto "${CMAKE_BUILD_TYPE}/Resources" "${RESOURCES}"

echo "\n--------\n" >> "${RESOURCES}/AboutThisBuild.txt"
echo "Bundle system: $(sysctl -n machdep.cpu.brand_string)" >> "${RESOURCES}/AboutThisBuild.txt"
echo "Bundle OS:     $(sw_vers -productName) $(sw_vers -productVersion) $(sw_vers -buildVersion) $(uname -mrs)" >> "${RESOURCES}/AboutThisBuild.txt"
echo "Bundle date:   $(date -Ru) UTC" >> "${RESOURCES}/AboutThisBuild.txt"
echo "Bundle epoch:  $(date +%s)" >> "${RESOURCES}/AboutThisBuild.txt"
echo "Bundle UUID:   $(uuidgen|tr 'A-Z' 'a-z')" >> "${RESOURCES}/AboutThisBuild.txt"

# Discover package locations through pkg-config where available. LOCAL_PREFIX
# remains the fallback for MacPorts and linked Homebrew formulae.
lensfun_prefix="$(pkg-config --variable=prefix lensfun 2>/dev/null)"
lensfun_libdir="$(pkg-config --variable=libdir lensfun 2>/dev/null)"
libpng_libdir="$(pkg-config --variable=libdir libpng 2>/dev/null)"
libtiff_libdir="$(pkg-config --variable=libdir libtiff-4 2>/dev/null)"
libsharpyuv_libdir="$(pkg-config --variable=libdir libsharpyuv 2>/dev/null)"
libjxl_libdir="$(pkg-config --variable=libdir libjxl 2>/dev/null)"

# Copy the Lensfun database and library into the app bundle.
lensfunversion="$(pkg-config --modversion lensfun | cut -f3 -d'.')"
if [[ "$lensfunversion" == 95 ]]; then
    lensfun_database_version="version_2"
    lensfun_library_name="liblensfun.2.dylib"
else
    lensfun_database_version="version_1"
    lensfun_library_name="liblensfun.1.dylib"
fi
resolve_directory "Lensfun ${lensfun_database_version} database" \
    "${LOCAL_PREFIX}/share/lensfun/${lensfun_database_version}" \
    "${lensfun_prefix}/share/lensfun/${lensfun_database_version}"
lensfun_data_dir="$RESOLVED_DIRECTORY"
resolve_file "$lensfun_library_name" \
    "${LOCAL_PREFIX}/lib/${lensfun_library_name}" \
    "${lensfun_libdir}/${lensfun_library_name}"
lensfun_library="$RESOLVED_FILE"
copy_tree "$lensfun_data_dir" "${RESOURCES}/share/lensfun"
copy_macho "$lensfun_library" "${LIB}/${lensfun_library_name}"

# libomp can be directly linked, nested under lib/ by MacPorts, or keg-only
# under opt/ in Homebrew. LOCAL_PREFIX is authoritative for every layout.
resolve_file "libomp.dylib" \
    "${LOCAL_PREFIX}/lib/libomp.dylib" \
    "${LOCAL_PREFIX}/lib/libomp/libomp.dylib" \
    "${LOCAL_PREFIX}/opt/libomp/lib/libomp.dylib" \
    "${LOCAL_PREFIX}"/libexec/llvm-*/lib/libomp.dylib
copy_macho "$RESOLVED_FILE" "${LIB}/libomp.dylib"

msg "Copying dependencies from ${GTK_PREFIX}."
CheckLink "${EXECUTABLE}" 2>&1

# Copy libraries needed by optional modules and codecs.
resolve_file "libpng16.16.dylib" \
    "${LOCAL_PREFIX}/lib/libpng16.16.dylib" \
    "${libpng_libdir}/libpng16.16.dylib"
copy_macho "$RESOLVED_FILE" "${LIB}/libpng16.16.dylib"

resolve_file "libgraphite2.3.dylib" \
    "${LOCAL_PREFIX}/lib/libgraphite2.3.dylib"
copy_macho "$RESOLVED_FILE" "${LIB}/libgraphite2.3.dylib"

resolve_file "libtiff.6.dylib" \
    "${LOCAL_PREFIX}/lib/libtiff.6.dylib" \
    "${libtiff_libdir}/libtiff.6.dylib"
copy_macho "$RESOLVED_FILE" "${LIB}/libtiff.6.dylib"

resolve_file "libsharpyuv.0.dylib" \
    "${LOCAL_PREFIX}/lib/libsharpyuv.0.dylib" \
    "${libsharpyuv_libdir}/libsharpyuv.0.dylib"
copy_macho "$RESOLVED_FILE" "${LIB}/libsharpyuv.0.dylib"

resolve_directory "libjxl library directory" \
    "$libjxl_libdir" \
    "${LOCAL_PREFIX}/lib"
libjxl_libdir="$RESOLVED_DIRECTORY"
libjxl_cms_count=0
for src in "${libjxl_libdir}"/libjxl_cms.*.dylib; do
    [[ -e "$src" ]] || continue
    copy_macho "$src" "${LIB}/$(basename "$src")"
    ((libjxl_cms_count += 1))
done
((libjxl_cms_count > 0)) ||
    die "No libjxl_cms libraries found in ${libjxl_libdir}"

# Copy GTK modules directly from their installed locations. Tree-copying these
# directories first can preserve package-manager symlinks and break them after
# the temporary hierarchy is removed.
GDK_LOADERS_DIR="$(pkg-config --variable=gdk_pixbuf_moduledir gdk-pixbuf-2.0)"
GTK_LIBDIR="$(pkg-config --variable=libdir gtk+-3.0)"

[[ -d "$GDK_LOADERS_DIR" ]] ||
    die "Could not locate the gdk-pixbuf loaders directory: ${GDK_LOADERS_DIR}"
[[ -d "${GTK_LIBDIR}/gtk-3.0" ]] ||
    die "Could not locate the GTK module root: ${GTK_LIBDIR}/gtk-3.0"

GTK_IMMODULES_DIR="$(
    find "${GTK_LIBDIR}/gtk-3.0" \
        -type d -name immodules -print -quit
)"

[[ -n "$GTK_IMMODULES_DIR" ]] ||
    die "Could not locate the GTK immodules directory"

gdk_loader_count=0
for src in "${GDK_LOADERS_DIR}"/*.so; do
    [[ -e "$src" ]] || continue

    copy_macho \
        "$src" \
        "${LIB}/$(basename "$src")"

    ((gdk_loader_count += 1))
done

((gdk_loader_count > 0)) ||
    die "No usable gdk-pixbuf loaders found in ${GDK_LOADERS_DIR}"

gtk_immodule_count=0
for src in \
    "${GTK_IMMODULES_DIR}"/*.so \
    "${GTK_IMMODULES_DIR}"/*.dylib
do
    [[ -e "$src" ]] || continue

    copy_macho \
        "$src" \
        "${LIB}/$(basename "$src")"

    ((gtk_immodule_count += 1))
done

((gtk_immodule_count > 0)) ||
    die "No usable GTK immodules found in ${GTK_IMMODULES_DIR}"

# GTK+3 themes
msg "Copy GTK+3 theme and icon resources:"
ditto {"${LOCAL_PREFIX}","${RESOURCES}"}/share/themes/Mac/gtk-3.0/gtk-keys.css
ditto {"${LOCAL_PREFIX}","${RESOURCES}"}/share/themes/Default/gtk-3.0/gtk-keys.css

# Adwaita icons
msg "Copy Adwaita icons"
copy_tree_dereference \
    "${LOCAL_PREFIX}/share/icons/Adwaita" \
    "${RESOURCES}/share/icons/Adwaita"
"${LOCAL_PREFIX}/bin/gtk-update-icon-cache" "${RESOURCES}/share/icons/Adwaita" || "${LOCAL_PREFIX}/bin/gtk-update-icon-cache-3.0" "${RESOURCES}/share/icons/Adwaita"
copy_tree_dereference \
    "${LOCAL_PREFIX}/share/icons/hicolor" \
    "${RESOURCES}/share/icons/hicolor"

# fix libfreetype install name
for lib in "${LIB}"/*; do
    install_name_tool -change libfreetype.6.dylib "${LIB}"/libfreetype.6.dylib "${lib}" 2>/dev/null
done

# Prepare GTK3 pixbuf loaders & immodules
msg "Build GTK3 databases:"
mkdir -p "${RESOURCES}"/share/gtk-3.0
mkdir -p "${ETC}"/gtk-3.0
# Change a relative path for the SVG pixbufloader
install_name_tool -delete_rpath @loader_path/../lib "${LIB}"/libpixbufloader_svg.so
install_name_tool -change @rpath/librsvg-2.2.dylib "${LOCAL_PREFIX}"/lib/librsvg-2.2.dylib  "${LIB}"/libpixbufloader_svg.so
# codesign Frameworks
sudo codesign --sign "${CODESIGNID}" --force "${LIB}"/*
# Build databases
"${LOCAL_PREFIX}"/bin/gdk-pixbuf-query-loaders "${LIB}"/libpixbufloader*[^dylib] > "${ETC}"/gtk-3.0/gdk-pixbuf.loaders
"${LOCAL_PREFIX}"/bin/gtk-query-immodules-3.0 "${LIB}"/im-* > "${ETC}"/gtk-3.0/gtk.immodules || "${LOCAL_PREFIX}"/bin/gtk-query-immodules "${LIB}"/im-* > "${ETC}"/gtk-3.0/gtk.immodules
sed -i.bak -e "s|${PWD}/RawTherapee.app/Contents/|/Applications/RawTherapee.app/Contents/|" "${ETC}"/gtk-3.0/gdk-pixbuf.loaders "${ETC}/gtk-3.0/gtk.immodules"
sed -i.bak -e "s|${LOCAL_PREFIX}/share/|/Applications/RawTherapee.app/Contents/Resources/share/|" "${ETC}"/gtk-3.0/gtk.immodules
sed -i.bak -e "s|${LOCAL_PREFIX}/|/Applications/RawTherapee.app/Contents/Frameworks/|" "${ETC}"/gtk-3.0/gtk.immodules
rm "${ETC}"/*/*.bak
# Remove a relative path for the SVG pixbufloader
install_name_tool -change @rpath/librsvg-2.2.dylib /Applications/RawTherapee.app/Contents/Frameworks/librsvg-2.2.dylib "${LIB}"/libpixbufloader_svg.so
# Modify the libpixbufloader_svg librsvg install_name
install_name_tool -change "${PWD}"/"${LIB}"/librsvg-2.2.dylib /Applications/"${LIB}"/librsvg-2.2.dylib "${LIB}"/libpixbufloader_svg.so

# Install names
ModifyInstallNames 2>/dev/null

# Mime directory
msg "Copying shared files from ${GTK_PREFIX}:"
copy_tree_dereference \
    "${LOCAL_PREFIX}/share/mime" \
    "${RESOURCES}/share/mime"

msg "Installing required application bundle files:"
if [[ -d "${PROJECT_SOURCE_DIR}/rtdata/fonts" ]]; then
    copy_tree "${PROJECT_SOURCE_DIR}/rtdata/fonts" "${ETC}/fonts"
else
    msg "No rtdata/fonts directory in this source tree; skipping it."
fi

# App bundle resources
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}/tools/osx"
ditto "${PROJECT_SOURCE_DATA_DIR}/"{rawtherapee,profile}.icns "${RESOURCES}"

update-mime-database -V  "${RESOURCES}/share/mime"
copy_tree_dereference \
    "${LOCAL_PREFIX}/share/locale" \
    "${RESOURCES}/share/locale"

msg "Build glib database:"
mkdir -p "${RESOURCES}/share/glib-2.0"
copy_tree_dereference \
    "${LOCAL_PREFIX}/share/glib-2.0/schemas" \
    "${RESOURCES}/share/glib-2.0/schemas"
"${LOCAL_PREFIX}/bin/glib-compile-schemas" "${RESOURCES}/share/glib-2.0/schemas"

# Append an LC_RPATH
msg "Registering @rpath into the main executable."
install_name_tool -add_rpath /Applications/"${LIB}" "${EXECUTABLE}" 2>/dev/null

ModifyInstallNames 2>/dev/null

# fix @rpath in Frameworks
msg "Registering @rpath in Frameworks folder."
for frameworklibs in "${LIB}"/*{dylib,so,cli}; do
    install_name_tool -delete_rpath ${LOCAL_PREFIX}/lib "${frameworklibs}" 2>/dev/null
    install_name_tool -add_rpath /Applications/"${LIB}" "${frameworklibs}" 2>/dev/null
done
install_name_tool -delete_rpath RawTherapee.app/Contents/Frameworks "${EXECUTABLE}"-cli 2>/dev/null
install_name_tool -add_rpath /Applications/"${LIB}" "${EXECUTABLE}"-cli 2>/dev/null

# Link to libomp instead of libgomp
sudo install_name_tool -change /Applications/RawTherapee.app/Contents/Frameworks/libgomp.1.dylib /Applications/RawTherapee.app/Contents/Frameworks/libomp.dylib RawTherapee.app/Contents/Frameworks/libfftw3f_omp.3.dylib
if [[ -e "${LIB}/libgomp.1.dylib" ]]; then
    rm "${LIB}/libgomp.1.dylib"
fi

# Merge the app with the other architecture to create the Universal app.
if [[ -n $UNIVERSAL_URL ]]; then
    msg "Getting Universal countercomponent."
    curl -L ${UNIVERSAL_URL} -o univ.zip
    msg "Extracting app."
    unzip univ.zip -d univapp
    hdiutil attach -mountpoint ./RawTherapeeuniv univapp/*folder/*dmg
    if [[ $arch = "arm64" ]]; then
        cp -R RawTherapee.app RawTherapee-arm64.app
        minimum_arm64_version=$(f=$(cat RawTherapee-arm64.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
        cp -R RawTherapeeuniv/RawTherapee.app RawTherapee-x86_64.app
        minimum_x86_64_version=$(f=$(cat RawTherapee-x86_64.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
        echo "\n\n=====================================\n\n" >> RawTherapee.app/Contents/Resources/AboutThisBuild.txt
        cat RawTherapee-x86_64.app/Contents/Resources/AboutThisBuild.txt >> RawTherapee.app/Contents/Resources/AboutThisBuild.txt
    else
        cp -R RawTherapee.app RawTherapee-x86_64.app
        minimum_x86_64_version=$(f=$(cat RawTherapee-x86_64.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
        cp -R RawTherapeeuniv/RawTherapee.app RawTherapee-arm64.app
        minimum_arm64_version=$(f=$(cat RawTherapee-arm64.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
        echo "\n\n=====================================\n\n" >> RawTherapee.app/Contents/Resources/AboutThisBuild.txt
        cat RawTherapee-arm64.app/Contents/Resources/AboutThisBuild.txt >> RawTherapee.app/Contents/Resources/AboutThisBuild.txt
    fi
    cmake -DPROJECT_SOURCE_DATA_DIR=${PROJECT_SOURCE_DATA_DIR} -DCONTENTS=${CONTENTS} -Dversion=${PROJECT_FULL_VERSION} -DshortVersion=${PROJECT_VERSION} -Dminimum_arm64_version=${minimum_arm64_version} -Dminimum_x86_64_version=${minimum_x86_64_version} -Darch=${arch} -P ${PROJECT_SOURCE_DATA_DIR}/info-plist.cmake
    plutil -convert xml1 ${APP}/Contents/Info.plist
    hdiutil unmount ./RawTherapeeuniv
    rm -r univapp
    # Create the fat main RawTherapee binary and move it into the new bundle
    lipo -create -output RawTherapee RawTherapee-arm64.app/Contents/MacOS/rawtherapee RawTherapee-x86_64.app/Contents/MacOS/rawtherapee
    lipo -create -output rawtherapee-cli RawTherapee-arm64.app/Contents/MacOS/rawtherapee-cli RawTherapee-x86_64.app/Contents/MacOS/rawtherapee-cli
    mv rawtherapee RawTherapee.app/Contents/MacOS
    # Create all the fat dependencies and move them into the bundle
    for lib in RawTherapee-arm64.app/Contents/Frameworks/* ; do
        lipo -create -output $(basename $lib) RawTherapee-arm64.app/Contents/Frameworks/$(basename $lib) RawTherapee-x86_64.app/Contents/Frameworks/$(basename $lib)
    done
    sudo mv *so *dylib RawTherapee.app/Contents/Frameworks
    sudo mv *-cli RawTherapee.app/Contents/MacOS
    rm -r RawTherapee-arm64.app
    rm -r RawTherapee-x86_64.app
else
    minimum_arm64_version=$(f=$(cat RawTherapee.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
    minimum_x86_64_version=${minimum_arm64_version}
        cmake -DPROJECT_SOURCE_DATA_DIR=${PROJECT_SOURCE_DATA_DIR} -DCONTENTS=${CONTENTS} -Dversion=${PROJECT_FULL_VERSION} -DshortVersion=${PROJECT_VERSION} -Dminimum_arm64_version=${minimum_arm64_version} -Dminimum_x86_64_version=${minimum_x86_64_version} -Darch=${arch} -P ${PROJECT_SOURCE_DATA_DIR}/info-plist.cmake
fi

# Package-manager links must never escape from Resources into a build-machine
# prefix. This catches that condition before codesign reports only the generic
# "invalid destination for symbolic link in bundle" error.
resource_symlink="$(find "${RESOURCES}" -type l -print -quit)"
[[ -z "$resource_symlink" ]] ||
    die "Application Resources still contains a symbolic link: ${resource_symlink}"

# Codesign the app
if [[ -n $CODESIGNID ]]; then
    msg "Codesigning Application."
    iconv -f UTF-8 -t ASCII "${PROJECT_SOURCE_DATA_DIR}"/rt.entitlements > "${CMAKE_BUILD_TYPE}"/rt.entitlements
    plutil -convert xml1 "${CMAKE_BUILD_TYPE}"/rt.entitlements
    for frame in ${APP}/Contents/Frameworks/* ; do
        echo $frame
        codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "$frame"
    done
    for resource in ${APP}/Contents/Resources/* ; do
        echo $resource
        if [ ! -d "$resource" ]; then
            codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "$resource"
        else
            for subresource in ${APP}/Contents/Resources/$(basename "$resource")/* ; do
                if [ ! -d "$subresource" ]; then
                    codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "$subresource"
                fi
            done
        fi
    done
    codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "${APP}"/Contents/MacOS/rawtherapee-cli
    codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "${APP}"/Contents/MacOS/rawtherapee
    codesign --preserve-metadata=identifier --digest-algorithm=sha256 --force --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "${APP}"
    spctl -a -vvvv "${APP}"
fi

# Notarize the app
if [[ -n $NOTARY ]]; then
    msg "Notarizing the application:"
    ditto -c -k --sequesterRsrc --keepParent "${APP}" "${APP}.zip"
    echo "Uploading..."
    sudo xcrun notarytool submit "${APP}.zip" ${NOTARY} --wait
    sudo xcrun stapler staple "${APP}"
fi

function CreateDmg {
    local srcDir="$(mktemp -dt $$.XXXXXXXXXXXX)"

    msg "Preparing disk image sources at ${srcDir}:"
    cp -R "${APP}" "${srcDir}"
    ln -s /Applications "${srcDir}"

    # Web bookmarks
    function CreateWebloc {
        defaults write "${srcDir}/$1" URL "$2"
        mv "${srcDir}/$1".{plist,webloc}
    }
    CreateWebloc       'Website' 'https://www.rawtherapee.com/'
    CreateWebloc 'Documentation' 'https://rawpedia.rawtherapee.com/'
    CreateWebloc         'Forum' 'https://discuss.pixls.us/c/software/rawtherapee'
    CreateWebloc    'Report Bug' 'https://github.com/RawTherapee/RawTherapee/issues/new?template=bug_report.md'

    # Disk image name
    if [[ -n $UNIVERSAL_URL ]]; then
        arch="Universal"
    fi
    dmg_name="${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}"
    lower_build_type="$(tr '[:upper:]' '[:lower:]' <<< "$CMAKE_BUILD_TYPE")"
    if [[ $lower_build_type != release ]]; then
        dmg_name="${dmg_name}_${lower_build_type}"
    fi

    msg "Creating disk image:"
    if [[ $FANCY_DMG == "ON" ]]; then
        echo "Building Fancy .dmg"
        touch message
        MESSAGE="$(cat message)"
        magick "${PROJECT_SOURCE_DATA_DIR}"/rtdmg-bkgd.png -pointsize 80 -font "/System/Library/Fonts/Supplemental/Arial.ttf" -fill Black -draw "text 14,1307 '${PROJECT_FULL_VERSION}'" -fill Salmon -font "/System/Library/Fonts/Supplemental/Arial.ttf" -draw "text 10,1300 '${PROJECT_FULL_VERSION}'" 1rtdmg-bkgd.png
        magick 1rtdmg-bkgd.png -pointsize 90 -fill Black -gravity center -font "/System/Library/Fonts/Supplemental/Arial.ttf" -draw "text 5,120 \"$MESSAGE\"" -fill Red -gravity center -font "/System/Library/Fonts/Supplemental/Arial.ttf" -draw "text 1,124 \"$MESSAGE\"" rtdmg-bkgd.png
        create-dmg \
        --background rtdmg-bkgd.png \
        --volname ${PROJECT_NAME}_${PROJECT_FULL_VERSION} \
        --volicon ${PROJECT_SOURCE_DATA_DIR}/rtdmg.icns \
        --window-pos 72 72 \
        --window-size 1000 697 \
        --text-size 16 \
        --icon-size 80 \
        --icon RawTherapee.app 250 238 \
        --icon Applications 700 238 \
        --icon Website.webloc 300 487 \
        --icon Forum.webloc 420 487 \
        --icon Report\ Bug.webloc 540 487 \
        --icon Documentation.webloc 680 487 \
        --no-internet-enable \
        --eula ${PROJECT_SOURCE_DATA_DIR}/../../LICENSE \
        --hdiutil-verbose \
        --hide-extension Website.webloc \
        --hide-extension Report\ Bug.webloc \
        --hide-extension Forum.webloc \
        --hide-extension Documentation.webloc \
        --filesystem APFS \
        ${dmg_name}.dmg ${srcDir}
    else
        hdiutil create -format UDBZ -fs HFS+ -srcdir "${srcDir}" -volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" "${dmg_name}.dmg"
    fi

    # Sign disk image
    if [[ -n $CODESIGNID ]]; then
        msg "Signing disk image"
        codesign  --digest-algorithm=sha256 --force -v -s "${CODESIGNID}" --timestamp "${dmg_name}.dmg"
    fi

    # Notarize the dmg
    if ! test -z "$NOTARY"; then
        msg "Notarizing the dmg:"
        zip "${dmg_name}.dmg.zip" "${dmg_name}.dmg"
        echo "Uploading..."
        sudo xcrun notarytool submit "${dmg_name}.dmg.zip" ${NOTARY} --wait
    fi

    # Zip disk image for redistribution
    msg "Zipping disk image for redistribution:"
    mkdir "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder"
    cp {"${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.dmg","${PROJECT_NAME}.app/Contents/MacOS/rawtherapee-cli","${PROJECT_SOURCE_DATA_DIR}/INSTALL.readme.rtf"} "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder"
    mv "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder/INSTALL.readme.rtf" "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder/install-readme.txt"
    codesign -s "${CODESIGNID}" -i com.rawtherapee.rawtherapee-cli -f "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder/rawtherapee-cli"
    zip -r "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.zip" "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder/"
    if [[ -n $NIGHTLY ]]; then
        cp "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.zip" "${PROJECT_NAME}_macOS_${arch}_latest.zip"
    fi
    if [[ -n $CONTINUOUS ]]; then
        BRANCH=$(git branch --show-current)
        if test -z "${BRANCH}"; then
            BRANCH=$(git rev-parse --short HEAD)
        fi
        mv "${PROJECT_NAME}_macOS_${arch}_latest.zip" "${PROJECT_NAME}_${BRANCH}_macOS_${arch}_${CMAKE_BUILD_TYPE}.zip"
    fi
}

CreateDmg
msg "Finishing build:"
# Clean up items
rm *app.zip
rm *dmg.zip
rm univ.zip
echo "Script complete."

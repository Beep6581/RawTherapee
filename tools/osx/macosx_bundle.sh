#!/usr/bin/env bash

# Required variables, must be set in rtdata/CMakeLists.txt
# - PROJECT_NAME
# - PROJECT_VERSION (if without git)
# - PROJECT_SOURCE_DIR
# - CMAKE_BUILD_TYPE
# - PROC_BIT_DEPTH
# - GTK_PREFIX

# Formatting
fNormal="$(tput sgr0)"
fBold="$(tput bold)"
# Colors depend upon the user's terminal emulator color scheme - what is readable for you may be not readable for someone else.
fMagenta="$(tput setaf 5)"
fRed="$(tput setaf 1)"

function msg {
    printf "\n${fBold}-- %s${fNormal}\n" "${@}"
}

function msgError {
    printf "\n${fBold}Error:${fNormal}\n%s\n" "${@}"
}

function GetDependencies {
    otool -L "$1" | awk 'NR >= 2 && $1 !~ /^(\/usr\/lib|\/System|@executable_path|@rpath)\// { print $1 }'
}

function CheckLink {
    GetDependencies "$1" | while read; do
        local dest="${LIB}/$(basename "${REPLY}")"
        test -f "${dest}" || { ditto --arch ${arch} "${REPLY}" "${dest}"; CheckLink "${dest}"; }
    done
}

# Source check
if [[ ! -d "${CMAKE_BUILD_TYPE}" ]]; then
    msgError "${PWD}/${CMAKE_BUILD_TYPE} folder does not exist. Please execute 'make install' first."
    exit 1
fi

# Update project version
if [[ -x "$(which git)" && -d "${PROJECT_SOURCE_DIR}/.git" ]]; then
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

MINIMUM_SYSTEM_VERSION="$(otool -l "${CMAKE_BUILD_TYPE}"/MacOS/rawtherapee | grep -A2 'LC_VERSION_MIN_MACOSX' | awk '$1 ~ /version/ { printf $2 }')"
if [[ -z "${MINIMUM_SYSTEM_VERSION}" ]]; then
    MINIMUM_SYSTEM_VERSION="$(sw_vers -productVersion | cut -d. -f-2)"
fi

case ${PROC_BIT_DEPTH} in
    64) arch=x86_64;;
    32) arch=i386;;
esac

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

APP="${PROJECT_NAME}.app"
CONTENTS="${APP}/Contents"
RESOURCES="${CONTENTS}/Resources"
MACOS="${CONTENTS}/MacOS"
LIB="${CONTENTS}/Frameworks"
ETC="${RESOURCES}/etc"
EXECUTABLE="${MACOS}/rawtherapee"

msg "Removing old files:"
rm -rf "${APP}" "${PROJECT_NAME}_*.dmg"

msg "Creating bundle container:"
install -d  "${RESOURCES}" \
    "${MACOS}" \
    "${LIB}" \
    "${ETC}"

msg "Copying release files:"
ditto "${CMAKE_BUILD_TYPE}/MacOS" "${MACOS}"
ditto "${CMAKE_BUILD_TYPE}/Resources" "${RESOURCES}"

msg "Copying dependencies from ${GTK_PREFIX}:"
CheckLink "${EXECUTABLE}"

msg "Copying library modules from ${GTK_PREFIX}:"
ditto --arch "${arch}" {"${GTK_PREFIX}/lib","${LIB}"}/gdk-pixbuf-2.0
ditto --arch "${arch}" {"${GTK_PREFIX}/lib","${LIB}"}/gtk-3.0

msg "Removing static libraries and cache files:"
find -E "${LIB}" -type f -regex '.*\.(a|la|cache)$' | while read; do rm "${REPLY}"; done

msg "Copying configuration files from ${GTK_PREFIX}:"
install -d "${ETC}/gtk-3.0"
cp "${GTK_PREFIX}/etc/gtk-3.0/im-multipress.conf" "${ETC}/gtk-3.0"
"${GTK_PREFIX}/bin/gdk-pixbuf-query-loaders" "${LIB}"/gdk-pixbuf-2.0/*/loaders/*.so > "${ETC}/gtk-3.0/gdk-pixbuf.loaders"
"${GTK_PREFIX}/bin/gtk-query-immodules-3.0"  "${LIB}"/gtk-3.0/*/immodules/*.so      > "${ETC}/gtk-3.0/gtk.immodules"
sed -i "" -e "s|${PWD}|/tmp|" "${ETC}/gtk-3.0/gdk-pixbuf.loaders" "${ETC}/gtk-3.0/gtk.immodules"

ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/glib-2.0/schemas
"${GTK_PREFIX}/bin/glib-compile-schemas" "${RESOURCES}/share/glib-2.0/schemas"

msg "Copying shared files from ${GTK_PREFIX}:"
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/mime
# GTK3 themes
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/themes/Mac/gtk-3.0/gtk-keys.css
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/themes/Default/gtk-3.0/gtk-keys.css
# Adwaita icons
iconfolders=("16x16/actions" "16x16/devices" "16x16/mimetypes" "16x16/places" "16x16/status" "48x48/devices")
for f in "${iconfolders[@]}"; do
    ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/icons/Adwaita/"$f"
done
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/icons/Adwaita/index.theme
"${GTK_PREFIX}/bin/gtk-update-icon-cache-3.0" "${RESOURCES}/share/icons/Adwaita"

### Pending deletion:
# fontconfig files (X11 backend only)
# if otool -L "${EXECUTABLE}" | grep -sq 'libgtk-x11-2.0'; then
#     msg "Installing fontconfig files (Using X11 backend. FONTCONFIG_PATH will be set by executable loader.)"
#     cp -RL "${GTK_PREFIX}/etc/fonts" "${ETC}"
# fi

# Install names
find -E "${CONTENTS}" -type f -regex '.*/(rawtherapee|.*\.(dylib|so))' | while read x; do
    msg "Modifying install names: ${x}"
    {
        # id
        case ${x} in *.dylib) echo "   install_name_tool -id '@rpath/$(basename "${x}")' '${x}'";; esac
        # names
        GetDependencies "${x}" | while read y; do
            echo "   install_name_tool -change '${y}' '@rpath/$(basename "${y}")' '${x}'"
        done
    } | bash -v
done

msg "Registering @loader_path into the executable:"
echo "   install_name_tool -add_rpath @loader_path/../Frameworks '${EXECUTABLE}'" | bash -v

msg "Installing required application bundle files:"
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}/tools/osx"

# Executable loader
# Note: executable is renamed to 'rawtherapee-bin'.
mv "${MACOS}/rawtherapee" "${MACOS}/rawtherapee-bin"
install -m 0755 "${PROJECT_SOURCE_DATA_DIR}/executable_loader.in" "${MACOS}/rawtherapee"
# App bundle resources
cp "${PROJECT_SOURCE_DATA_DIR}/"{rawtherapee,profile}.icns "${RESOURCES}"
cp "${PROJECT_SOURCE_DATA_DIR}/PkgInfo" "${CONTENTS}"
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}/Info.plist.in" "${CONTENTS}/Info.plist"
sed -i "" -e "s|@version@|${PROJECT_FULL_VERSION}|
s|@shortVersion@|${PROJECT_VERSION}|
s|@arch@|${arch}|" \
    "${CONTENTS}/Info.plist"
plutil -convert binary1 "${CONTENTS}/Info.plist"



function CreateDmg {
    local srcDir="$(mktemp -dt $$)"

    msg "Preparing disk image sources at ${srcDir}:"
    cp -R "${APP}" "${srcDir}"
    cp AboutThisBuild.txt "${srcDir}"
    ln -s /Applications "${srcDir}"

    # Web bookmarks
    function CreateWebloc {
        defaults write "${srcDir}/$1" URL "$2"
        mv "${srcDir}/$1".{plist,webloc}
    }
    CreateWebloc 'Website' 'http://www.rawtherapee.com/'
    CreateWebloc 'Manual'  'http://rawpedia.rawtherapee.com/'

    # Disk image name
    dmg_name="${PROJECT_NAME// /_}_OSX_${MINIMUM_SYSTEM_VERSION}_${PROC_BIT_DEPTH}_${PROJECT_FULL_VERSION}"
    lower_build_type="$(tr '[:upper:]' '[:lower:]' <<< "$CMAKE_BUILD_TYPE")"
    if [[ ${lower_build_type} != release ]]; then
        dmg_name="${dmg_name}_${lower_build_type}"
    fi

    msg "Creating disk image:"
    hdiutil create -format UDBZ -srcdir "${srcDir}" -volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" "${dmg_name}.dmg"

    # Zip disk image for redistribution
    zip "${dmg_name}.zip" "${dmg_name}.dmg" AboutThisBuild.txt
    rm "${dmg_name}.dmg"

    msg "Removing disk image caches:"
    rm -rf "${srcDir}"
}
CreateDmg

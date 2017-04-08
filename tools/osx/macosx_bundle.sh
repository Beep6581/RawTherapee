#!/bin/bash

# Required variables
# ------------------
# these are very important variables. Must be set into rtdata/CMakeLists.txt!
# - PROJECT_NAME
# - PROJECT_SOURCE_DIR
# - PROJECT_VERSION (if without git)
# - CMAKE_BUILD_TYPE
# - PROC_BIT_DEPTH
# - GTK_PREFIX

function message {
    printf '\e[35m-- %s\e[m\n' "$*"
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

# source check
if test ! -d "${CMAKE_BUILD_TYPE}"; then
    printf "\e[31m${PWD}/${CMAKE_BUILD_TYPE} directory is not found. Please execute 'make install' first.\e[m\n"
    exit 1
fi

# update project version
if test -x "$(which git)" -a -d "${PROJECT_SOURCE_DIR}/.git"; then
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

# if not specify CMAKE_OSX_DEPLOYMENT_TARGET when compiling,
# 'MINIMUM_VERSION' will be used host OS X version.
MINIMUM_SYSTEM_VERSION="$(otool -l "${CMAKE_BUILD_TYPE}"/MacOS/rawtherapee | grep -A2 'LC_VERSION_MIN_MACOSX' | awk '$1 ~ /version/ { printf $2 }')"
if test ! -n "${MINIMUM_SYSTEM_VERSION}"; then
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

APP="${PROJECT_NAME}".app
CONTENTS="${APP}"/Contents
RESOURCES="${CONTENTS}"/Resources
MACOS="${CONTENTS}"/MacOS
LIB="${CONTENTS}"/Frameworks
ETC="${RESOURCES}"/etc
EXECUTABLE="${MACOS}"/rawtherapee

message "Removing old files"
rm -rf "${APP}" "${PROJECT_NAME}_*.dmg"

message "Creating bundle container"
install -d  "${RESOURCES}" \
"${MACOS}" \
"${LIB}" \
"${ETC}"

message "Copying release files"
ditto "${CMAKE_BUILD_TYPE}"/MacOS "${MACOS}"
ditto "${CMAKE_BUILD_TYPE}"/Resources "${RESOURCES}"

message "Copying dependencies from ${GTK_PREFIX}"
CheckLink "${EXECUTABLE}"

message "Copying library modules from ${GTK_PREFIX}"
ditto --arch "${arch}" {"${GTK_PREFIX}"/lib,"${LIB}"}/gdk-pixbuf-2.0
ditto --arch "${arch}" {"${GTK_PREFIX}"/lib,"${LIB}"}/gtk-3.0

message "Removing static libraries and cache files"
find -E "${LIB}" -type f -regex '.*\.(a|la|cache)$' | while read; do rm "${REPLY}"; done

message "Copying configuration files from ${GTK_PREFIX}"
install -d "${ETC}"/gtk-3.0
cp "${GTK_PREFIX}"/etc/gtk-3.0/im-multipress.conf "${ETC}"/gtk-3.0
"${GTK_PREFIX}"/bin/gdk-pixbuf-query-loaders "${LIB}"/gdk-pixbuf-2.0/*/loaders/*.so > "${ETC}"/gtk-3.0/gdk-pixbuf.loaders
"${GTK_PREFIX}"/bin/gtk-query-immodules-3.0  "${LIB}"/gtk-3.0/*/immodules/*.so      > "${ETC}"/gtk-3.0/gtk.immodules
sed -i "" -e "s|${PWD}|/tmp|" "${ETC}"/gtk-3.0/gdk-pixbuf.loaders \
"${ETC}"/gtk-3.0/gtk.immodules

ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/glib-2.0/schemas
"${GTK_PREFIX}"/bin/glib-compile-schemas "${RESOURCES}"/share/glib-2.0/schemas

message "Copying shared files from ${GTK_PREFIX}"
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/mime
# gtk3 themes
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/themes/Mac/gtk-3.0/gtk-keys.css
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/themes/Default/gtk-3.0/gtk-keys.css
# Adwaita icons
iconfolders=("16x16/actions" "16x16/devices" "16x16/mimetypes" "16x16/places" "16x16/status" "48x48/devices")
for f in "${iconfolders[@]}"; do ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/icons/Adwaita/"$f"; done
ditto {"${GTK_PREFIX}","${RESOURCES}"}/share/icons/Adwaita/index.theme
"${GTK_PREFIX}"/bin/gtk-update-icon-cache-3.0 "${RESOURCES}"/share/icons/Adwaita
# fontconfig files (X11 backend only)
# if otool -L "${EXECUTABLE}" | grep -sq 'libgtk-x11-2.0'; then
#     message "Installing fontconfig files (Your library is X11 backend. 'FONTCONFIG_PATH' will be set by executable loader.)"
#     cp -RL "${GTK_PREFIX}"/etc/fonts "${ETC}"
# fi



# install names
find -E "${MACOS}" -type f -regex '.*/(rawtherapee|.*\.(dylib|so))' | while read x; do
    message "Modifying install names: ${x}"
    {
        # id
    case ${x} in *.dylib) echo "   install_name_tool -id '@rpath/$(basename "${x}")' '${x}'";; esac
    # names
    GetDependencies "${x}" | while read y; do
        echo "   install_name_tool -change '${y}' '@rpath/$(basename "${y}")' '${x}'"
    done
} | bash -v
done

message "Registering @loader_path into the executable"
echo "   install_name_tool -add_rpath @loader_path/lib '${EXECUTABLE}'" | bash -v



message "Installing required application bundle files"
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}"/tools/osx
# executable loader
# note: executable is renamed to 'rawtherapee-bin'.
mv "${MACOS}"/rawtherapee{,-bin}
install -m 0755 "${PROJECT_SOURCE_DATA_DIR}"/executable_loader.in "${MACOS}"/rawtherapee
# app bundle resources
cp "${PROJECT_SOURCE_DATA_DIR}"/{rawtherapee,profile}.icns "${RESOURCES}"
cp "${PROJECT_SOURCE_DATA_DIR}"/PkgInfo "${CONTENTS}"
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}"/Info.plist.in "${CONTENTS}"/Info.plist
sed -i "" -e "s|@version@|${PROJECT_FULL_VERSION}|
s|@shortVersion@|${PROJECT_VERSION}|
s|@arch@|${arch}|" \
"${CONTENTS}"/Info.plist
plutil -convert binary1 "${CONTENTS}"/Info.plist



function CreateDmg {
    local srcdir="$(mktemp -dt $$)"

    message "Preparing disk image sources at ${srcdir}"
    cp -R "${APP}" "${srcdir}"
    cp AboutThisBuild.txt "${srcdir}"
    ln -s /Applications "${srcdir}"

    # web bookmarks
    function CreateWebloc {
        defaults write "${srcdir}/$1" URL "$2"
        mv "${srcdir}/$1".{plist,webloc}
    }
    CreateWebloc 'RawTherapee Blog' 'http://www.rawtherapee.com'
    CreateWebloc 'Online Manual'    'http://rawpedia.rawtherapee.com/'

    # disk image name
    dmg_name="${PROJECT_NAME// /_}_OSX_${MINIMUM_SYSTEM_VERSION}_${PROC_BIT_DEPTH}_${PROJECT_FULL_VERSION}"
    if ! echo "${CMAKE_BUILD_TYPE}" | grep -sqi "release"; then
        dmg_name="${dmg_name}_$(echo ${CMAKE_BUILD_TYPE} | tr '[:upper:]' '[:lower:]')"
    fi

    message "Creating disk image"
    hdiutil create -format UDBZ -srcdir "${srcdir}" -volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" "${dmg_name}".dmg

    # zip .dmg for re-distribution
    zip "${dmg_name}.zip" "${dmg_name}.dmg" AboutThisBuild.txt
    rm "${dmg_name}.dmg"

    message "Removing disk image caches"
    rm -rf "${srcdir}"
}
CreateDmg

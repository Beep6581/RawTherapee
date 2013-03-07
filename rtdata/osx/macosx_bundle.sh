#!/bin/bash

# Required variables
# ------------------
# these are very important variables. Must be set into rtdata/CMakeLists.txt!
# - PROJECT_NAME
# - PROJECT_SOURCE_DIR
# - PROJECT_VERSION (if without mercurial)
# - CMAKE_BUILD_TYPE
# - PROC_BIT_DEPTH
# - GTK_PREFIX


# --------------------------------------
# Functions
# --------------------------------------
function message {
  printf '\e[34m-- %s\e[m\n' "$*"
}
function get_dependencies {
  otool -L "$1" | awk 'NR >= 2 && $1 !~ /^(\/usr\/lib|\/System|@executable_path|@rpath)\// { print $1 }'
}
function checkLink {
  get_dependencies "$1" |
  while read x
  do
    if test ! -f "${LIB}/${x##*/}"; then
      ditto --arch ${arch} "${x}" "$_"
      checkLink "$_"
    fi
  done
}


# source check
if test ! -d "${CMAKE_BUILD_TYPE}"; then
  printf "\e[31m${PWD}/${CMAKE_BUILD_TYPE} directory is not found. Please execute 'make install' first.\e[m\n"
  exit 1
fi

# update project version
if test -x `which hg` -a -d "${PROJECT_SOURCE_DIR}/.hg"; then
  PROJECT_VERSION=`hg -R "${PROJECT_SOURCE_DIR}" parents --template "{latesttag}.{latesttagdistance}"`
fi

# if not specify CMAKE_OSX_DEPLOYMENT_TARGET when compiling,
# 'MINIMUM_VERSION' will be used host OS X version.
MINIMUM_SYSTEM_VERSION=`otool -l "${CMAKE_BUILD_TYPE}"/rawtherapee | grep -A2 'LC_VERSION_MIN_MACOSX' | awk '$1 ~ /version/ { printf $2 }'`
if test ! -n "${MINIMUM_SYSTEM_VERSION}"; then
  MINIMUM_SYSTEM_VERSION=`sw_vers -productVersion | cut -d. -f-2`
fi

case ${PROC_BIT_DEPTH} in
  64) arch=x86_64;;
  32) arch=i386;;
esac

echo "PROJECT_NAME:           ${PROJECT_NAME}"
echo "PROJECT_VERSION:        ${PROJECT_VERSION}"
echo "PROJECT_SOURCE_DIR:     ${PROJECT_SOURCE_DIR}"
echo "CMAKE_BUILD_TYPE:       ${CMAKE_BUILD_TYPE}"
echo "PROC_BIT_DEPTH:         ${PROC_BIT_DEPTH}"
echo "MINIMUM_SYSTEM_VERSION: ${MINIMUM_SYSTEM_VERSION}"
echo "GTK_PREFIX:             ${GTK_PREFIX}"
echo "PWD:                    ${PWD}"

APP="${PROJECT_NAME}".app
CONTENTS="${APP}"/Contents
RESOURCES="${CONTENTS}"/Resources
MACOS="${CONTENTS}"/MacOS
LIB="${MACOS}"/lib
ETC="${MACOS}"/etc
EXECUTABLE="${MACOS}"/rawtherapee

message "Removing old files"
rm -rf "${APP}" ${PROJECT_NAME}_*.dmg

message "Creating bundle container"
install -d  "${RESOURCES}" \
            "${MACOS}" \
            "${LIB}" \
            "${ETC}"

# set bundle attribute
setFile -a B "${APP}"

message "Copying release files"
ditto --arch ${arch} "${CMAKE_BUILD_TYPE}" "${MACOS}"

message "Copying dependencies from ${GTK_PREFIX}"
checkLink "${EXECUTABLE}"

message "Copying library modules from ${GTK_PREFIX}"
ditto --arch ${arch} {"${GTK_PREFIX}"/lib,"${LIB}"}/gdk-pixbuf-2.0
ditto --arch ${arch} {"${GTK_PREFIX}"/lib,"${LIB}"}/gtk-2.0
ditto --arch ${arch} {"${GTK_PREFIX}"/lib,"${LIB}"}/pango

message "Removing static libraries and cache files"
find -E "${LIB}" -type f -regex '.*\.(a|la|cache)$' | while read; do rm "${REPLY}"; done

message "Copying configuration files from ${GTK_PREFIX}"
install -d "${ETC}"/{gtk-2.0,pango}
cp "${GTK_PREFIX}"/etc/gtk-2.0/im-multipress.conf "${ETC}"/gtk-2.0
"${GTK_PREFIX}"/bin/gdk-pixbuf-query-loaders "${LIB}"/gdk-pixbuf-2.0/*/loaders/*.so > "${ETC}"/gtk-2.0/gdk-pixbuf.loaders
"${GTK_PREFIX}"/bin/gtk-query-immodules-2.0  "${LIB}"/gtk-2.0/*/immodules/*.so      > "${ETC}"/gtk-2.0/gtk.immodules
"${GTK_PREFIX}"/bin/pango-querymodules       "${LIB}"/pango/*/modules/*.so          > "${ETC}"/pango/pango.modules
sed -i "" -e "s|${PWD}|/tmp|" "${ETC}"/gtk-2.0/gdk-pixbuf.loaders \
                              "${ETC}"/gtk-2.0/gtk.immodules \
                              "${ETC}"/pango/pango.modules
printf "[Pango]\nModuleFiles = /tmp/${ETC}/pango/pango.modules" > "${ETC}"/pango/pangorc

message "Copying shared files from ${GTK_PREFIX}"
cp -R "${GTK_PREFIX}"/share/mime "${MACOS}"/share
# gtk themes
ditto {/opt/local,"${MACOS}"}/share/themes/Mac/gtk-2.0-key/gtkrc
ditto {/opt/local,"${MACOS}"}/share/themes/Clearlooks/gtk-2.0/gtkrc
install -d "${MACOS}"/share/themes/Raleigh/gtk-2.0 && (
  cd "$_" &&
  ln -s ../../Clearlooks/gtk-2.0/gtkrc
)

# fontconfig files (X11 backend only)
if otool -L "${EXECUTABLE}" | grep -sq 'libgtk-x11-2.0'; then
    message "Installing fontconfig files (Your library is X11 backend. 'FONTCONFIG_PATH' will be set by executable loader.)"
    cp -RL "${GTK_PREFIX}"/etc/fonts "${ETC}"
fi

# --------------------------------------
# Install names
# --------------------------------------
## install_name_tool
# this command replaces install names. you will not need to set 'DYLD_LIBRARY_PATH'
# '@rpath' is replaced '@loader_path' when execute.

# @rpath
find -E "${MACOS}" -type f -regex '.*(rawtherapee|\.dylib|\.so)$' | while read x; do
  message "Modifying install names: ${x}"
  {
    # id (dylib file only)
    case "${x}" in (*.dylib)
      printf '   install_name_tool -id "%s" "%s"\n' "@rpath/${x##*/}" "${x}";;
    esac
    
    # dependencies
    get_dependencies "${x}" | while read y; do
      printf '   install_name_tool -change "%s" "%s" "%s"\n' "${y}" "@rpath/${y##*/}" "${x}"
    done
  } | bash -v
done

# @loader_path
message "Registering @loader_path into the executable"
printf '   install_name_tool -add_rpath @loader_path/lib "%s"' "${EXECUTABLE}" | bash -v


# --------------------------------------
# required files
# --------------------------------------
message "Installing required application bundle files"
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}"/rtdata/osx

## executable and executable loader
# note: executable is renamed to 'rawtherapee-bin'.
mv "${MACOS}"/rawtherapee{,-bin}
install -m 0755 "${PROJECT_SOURCE_DATA_DIR}"/executable_loader.in "${MACOS}"/rawtherapee
# icons
cp "${PROJECT_SOURCE_DATA_DIR}"/{rawtherapee,profile}.icns "${RESOURCES}"
# PkgInfo
cp "${PROJECT_SOURCE_DATA_DIR}"/PkgInfo "${CONTENTS}"
# Info.plist
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}"/Info.plist.in "${CONTENTS}"/Info.plist
# configure property values
sed -i "" -e "s|@version@|${PROJECT_VERSION}|" \
          -e "s|@shortVersion@|`echo ${PROJECT_VERSION} | cut -d. -f-3`|" \
          -e "s|@arch@|${arch}|" \
          "${CONTENTS}"/Info.plist
# specify min version
#defaults write "${PWD}/${CONTENTS}"/Info LSMinimumSystemVersion "${MINIMUM_SYSTEM_VERSION}"
plutil -convert binary1 "${CONTENTS}"/Info.plist



# --------------------------------------
# disk image
# --------------------------------------
## cache directory
DMG_SOURCE_DIR=/tmp/`uuidgen`

## disk image name
# space and comma will be replaced.
# example: RawTherapee_OSX108_4.0.10.1_Release.dmg
DMG="${PROJECT_NAME// /_}_OSX${MINIMUM_SYSTEM_VERSION//.}_${PROJECT_VERSION}_${CMAKE_BUILD_TYPE}".dmg
# --------------------------------------
message "Preparing disk image sources at ${DMG_SOURCE_DIR}"
install -d ${DMG_SOURCE_DIR}
mv "${APP}" ${DMG_SOURCE_DIR}
cp AboutThisBuild.txt ${DMG_SOURCE_DIR}
# symlink to /Applications dir
ln -s /Applications ${DMG_SOURCE_DIR}
# .webloc files
defaults write ${DMG_SOURCE_DIR}/"RawTherapee Blog" URL 'http://www.rawtherapee.com'
mv ${DMG_SOURCE_DIR}/"RawTherapee Blog".{plist,webloc}
defaults write ${DMG_SOURCE_DIR}/"Online Manual" URL 'https://docs.google.com/document/d/1DHLb_6xNQsEInxiuU8pz1-sWNinnj09bpBUA4_Vl8w8/edit'
mv ${DMG_SOURCE_DIR}/"Online Manual".{plist,webloc}

message "Creating disk image"
## hdiutil
# -format: defaut is 'UDZO'(zip compress). 'UDBZ' is bzip2 compress.
# -volname: example, RawTherapee_4.0.10.1
hdiutil create  -format UDBZ \
                -srcdir ${DMG_SOURCE_DIR} \
                -volname "${PROJECT_NAME}_${PROJECT_VERSION}" \
                "${DMG}"

message "Cleaning disk image caches"
rm -rf ${DMG_SOURCE_DIR}

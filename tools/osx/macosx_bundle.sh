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
    printf "\\n${fBold}-- %s${fNormal}\\n" "${@}"
}

function msgError {
    printf "\\n${fBold}Error:${fNormal}\\n%s\\n" "${@}"
}

function GetDependencies {
    otool -L "$1" | awk 'NR >= 2 && $1 !~ /^(\/usr\/lib|\/System|@executable_path|@rpath)\// { print $1 }'
}

function CheckLink {
    GetDependencies "$1" | while read -r; do
        local dest="${LIB}/$(basename "${REPLY}")"
        test -f "${dest}" || { ditto --arch "${arch}" "${REPLY}" "${dest}"; CheckLink "${dest}"; }
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

arch=x86_64

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

LOCAL_PREFIX="$(cmake .. -LA -N | grep "LOCAL_PREFIX" | cut -d "=" -f2)"
EXPATLIB="$(cmake .. -LA -N | grep "pkgcfg_lib_EXPAT_expat" | cut -d "=" -f2)"
 
APP="${PROJECT_NAME}.app"
CONTENTS="${APP}/Contents"
RESOURCES="${CONTENTS}/Resources"
MACOS="${CONTENTS}/MacOS"
LIB="${CONTENTS}/Frameworks"
ETC="${RESOURCES}/etc"
EXECUTABLE="${MACOS}/rawtherapee"
GDK_PREFIX="${LOCAL_PREFIX}/local/"

msg "Removing old files:"
rm -rf "${APP}" "${PROJECT_NAME}_*.dmg" "*zip"

msg "Creating bundle container:"
install -d  "${RESOURCES}" \
    "${MACOS}" \
    "${LIB}" \
    "${ETC}"

msg "Copying release files:"
ditto "${CMAKE_BUILD_TYPE}/MacOS" "${MACOS}"
ditto "Resources" "${RESOURCES}"

msg "Copying dependencies from ${GTK_PREFIX}:"
CheckLink "${EXECUTABLE}"

msg "Copying library modules from ${GTK_PREFIX}:"
ditto --arch "${arch}" {"${GTK_PREFIX}/lib","${LIB}"}/gdk-pixbuf-2.0
ditto --arch "${arch}" {"${GTK_PREFIX}/lib","${LIB}"}/gtk-3.0

msg "Removing static libraries and cache files:"
find -E "${LIB}" -type f -regex '.*\.(a|la|cache)$' | while read -r; do rm "${REPLY}"; done

msg "Copying configuration files from ${GTK_PREFIX}:"
install -d "${ETC}/gtk-3.0"

# Make Frameworks folder flat
ditto "${LIB}"/gdk-pixbuf-2.0/2*/loaders/*.so "${LIB}"
ditto "${LIB}"/gtk-3.0/3*/immodules/*.{dylib,so} "${LIB}"
rm -r "${LIB}"/gtk-3.0
rm -r "${LIB}"/gdk-pixbuf-2.0

"${LOCAL_PREFIX}"/local/bin/gdk-pixbuf-query-loaders "${LIB}"/libpix*.so > "${ETC}"/gtk-3.0/gdk-pixbuf.loaders
"${LOCAL_PREFIX}"/local/bin/gtk-query-immodules-3.0 "${LIB}"/im-* > "${ETC}"/gtk-3.0/gtk.immodules
sed -i "" -e "s|${PWD}/RawTherapee.app/Contents/|/Applications/RawTherapee.app/Contents/|" "${ETC}/gtk-3.0/gdk-pixbuf.loaders" "${ETC}/gtk-3.0/gtk.immodules"
sed -i "" -e "s|/opt/local/|/usr/|" "${ETC}/gtk-3.0/gtk.immodules"

mkdir -p ${RESOURCES}/share/glib-2.0
ditto {"${LOCAL_PREFIX}/local","${RESOURCES}"}/share/glib-2.0/schemas
"${LOCAL_PREFIX}/local/bin/glib-compile-schemas" "${RESOURCES}/share/glib-2.0/schemas"

msg "Copying shared files from ${GTK_PREFIX}:"
ditto {"${LOCAL_PREFIX}/local","${RESOURCES}"}/share/mime

# GTK3 themes
ditto {"${LOCAL_PREFIX}/local","${RESOURCES}"}/share/themes/Mac/gtk-3.0/gtk-keys.css
ditto {"${LOCAL_PREFIX}/local","${RESOURCES}"}/share/themes/Default/gtk-3.0/gtk-keys.css
# Adwaita icons
iconfolders=("16x16/actions" "16x16/devices" "16x16/mimetypes" "16x16/places" "16x16/status" "48x48/devices")
for f in "${iconfolders[@]}"; do
    mkdir -p ${RESOURCES}/share/icons/Adwaita/${f}
    ditto ${LOCAL_PREFIX}/local/share/icons/Adwaita/${f}/* "${RESOURCES}"/share/icons/Adwaita/${f}
done
ditto {"${LOCAL_PREFIX}/local","${RESOURCES}"}/share/icons/Adwaita/index.theme
"${LOCAL_PREFIX}/local/bin/gtk-update-icon-cache" "${RESOURCES}/share/icons/Adwaita"
ditto "${LOCAL_PREFIX}/local/share/icons/hicolor" "${RESOURCES}/share/icons/hicolor"

# Copy libjpeg-turbo ("62") into the app bundle
ditto ${LOCAL_PREFIX}/local/lib/libjpeg.62.dylib "${CONTENTS}/Frameworks/libjpeg.62.dylib"

# Copy libexpat into the app bundle (which is keg-only)

if [[ -d /usr/local/Cellar/expat ]]; then ditto /usr/local/Cellar/expat/*/lib/libexpat.1.dylib "${CONTENTS}/Frameworks"; else ditto "${EXPATLIB}" "${CONTENTS}/Frameworks/libexpat.1.dylib"; fi

# Copy libz into the app bundle
ditto ${LOCAL_PREFIX}/local/lib/libz.1.dylib "${CONTENTS}/Frameworks"

# Copy libtiff 5 into the app bundle
ditto ${LOCAL_PREFIX}/local/lib/libtiff.5.dylib "${CONTENTS}/Frameworks/libtiff.5.dylib"

# Copy the Lensfun database into the app bundle
mkdir -p "${RESOURCES}/share/lensfun"
ditto ${LOCAL_PREFIX}/local/share/lensfun/version_2/* "${RESOURCES}/share/lensfun"

# Copy liblensfun to Frameworks
ditto ${LOCAL_PREFIX}/local/lib/liblensfun.2.dylib "${CONTENTS}/Frameworks/liblensfun.2.dylib"

# Copy libomp to Frameworks
ditto ${LOCAL_PREFIX}/local/lib/libomp.dylib "${CONTENTS}/Frameworks"

# Install names
find -E "${CONTENTS}" -type f -regex '.*/(rawtherapee-cli|rawtherapee|.*\.(dylib|so))' | while read -r x; do
    msg "Modifying install names: ${x}"
    {
        # id
        case ${x} in *.dylib) echo "   install_name_tool -id '@rpath/$(basename "${x}")' '${x}'";; esac
        # names
        GetDependencies "${x}" | while read -r y; do
            echo "   install_name_tool -change '${y}' '@rpath/$(basename "${y}")' '${x}'"
        done
    } | bash -v
done



msg "Installing required application bundle files:"
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}/tools/osx"
ditto "${PROJECT_SOURCE_DIR}/build/Resources" "${RESOURCES}"
# Executable loader
# Note: executable is renamed to 'rawtherapee-bin'.
mkdir "${MACOS}/bin"
ditto "${MACOS}/rawtherapee" "${MACOS}/bin/rawtherapee-bin"
rm "${MACOS}/rawtherappe"
install -m 0755 "${PROJECT_SOURCE_DATA_DIR}/executable_loader.in" "${MACOS}/rawtherapee"
# App bundle resources
ditto "${PROJECT_SOURCE_DATA_DIR}/"{rawtherapee,profile}.icns "${RESOURCES}"
ditto "${PROJECT_SOURCE_DATA_DIR}/PkgInfo" "${CONTENTS}"
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}/Info.plist.in" "${CONTENTS}/Info.plist"
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}/Info.plist-bin.in" "${CONTENTS}/MacOS/bin/Info.plist"
sed -i "" -e "s|@version@|${PROJECT_FULL_VERSION}|
s|@shortVersion@|${PROJECT_VERSION}|
s|@arch@|${arch}|" \
    "${CONTENTS}/Info.plist"
plutil -convert xml1 "${CONTENTS}/Info.plist"
plutil -convert xml1 "${CONTENTS}/MacOS/bin/Info.plist"
update-mime-database -V  "${CONTENTS}/Resources/share/mime"

msg "Registering @rpath into the executable:"
echo "   install_name_tool -add_rpath /Applications/RawTherapee.app/Contents/Frameworks '${MACOS}/bin/rawtherapee-bin'" | bash -v
echo "   install_name_tool -add_rpath /Applications/RawTherapee.app/Contents/Frameworks '${EXECUTABLE}-cli'" | bash -v
for frameworklibs in ${CONTENTS}/Frameworks/* ; do
    echo "   install_name_tool -delete_rpath /opt/local/lib '${frameworklibs}'" | bash -v
    echo "   install_name_tool -add_rpath /Applications/RawTherapee.app/Contents/Frameworks '${frameworklibs}'" | bash -v
done

# Sign the app
CODESIGNID="$(cmake .. -LA -N | grep "CODESIGNID" | cut -d "=" -f2)"
if ! test -z "$CODESIGNID" ; then
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}/rt.entitlements" "${CONTENTS}/Entitlements.plist"
plutil -convert xml1 "${CONTENTS}/Entitlements.plist"
install -m 0644 "${PROJECT_SOURCE_DATA_DIR}/rt-bin.entitlements" "${CONTENTS}/MacOS/bin/Entitlements.plist"
plutil -convert xml1 "${CONTENTS}/MacOS/bin/Entitlements.plist"
codesign -v -s "${CODESIGNID}" -i "com.rawtherapee.rawtherapee-bin" -o runtime --timestamp --entitlements "${APP}/Contents/MacOS/bin/Entitlements.plist" "${APP}/Contents/MacOS/bin/rawtherapee-bin"
for frameworklibs in ${CONTENTS}/Frameworks/* ; do
    codesign -v -s "${CODESIGNID}" -i "com.rawtherapee.rawtherapee-bin" -o runtime --timestamp "${frameworklibs}"
done
codesign --deep --preserve-metadata=identifier,entitlements,runtime --timestamp --strict -v -s "${CODESIGNID}" -i "com.rawtherapee.rawtherapee" -o runtime --entitlements "${CONTENTS}/Entitlements.plist" "${APP}"
spctl -a -vvvv "${APP}"
fi

# Notarize the app

NOTARY="$(cmake .. -LA -N | grep "NOTARY" | cut -d "=" -f2)"
if ! test -z "$NOTARY" ; then
    msg "Notarizing the application:"
    ditto -c -k --sequesterRsrc --keepParent "${APP}" "${APP}.zip"
    uuid=`xcrun altool --notarize-app --primary-bundle-id "com.rawtherapee.rawtherapee" ${NOTARY} --file "${APP}.zip" 2>&1 | grep 'RequestUUID' | awk '{ print $3 }'`
    echo "Result= $uuid" # Display identifier string
    sleep 15
    while :
        do
        fullstatus=`xcrun altool --notarization-info "$uuid" ${NOTARY}  2>&1`  # get the status
        status1=`echo "$fullstatus" | grep 'Status\:' | awk '{ print $2 }'`
        if [ "$status1" = "success" ]; then
            xcrun stapler staple *app   #  staple the ticket
            xcrun stapler validate -v *app
            echo "Notarization success"
            break
        elif [ "$status1" = "in" ]; then
            echo "Notarization still in progress, sleeping for 15 seconds and trying again"
            sleep 15
        else
            echo "Notarization failed fullstatus below"
            echo "$fullstatus"
            exit 1
        fi
    done
fi

function CreateDmg {
    local srcDir="$(mktemp -dt $$)"

    msg "Preparing disk image sources at ${srcDir}:"
    cp -R "${APP}" "${srcDir}"
    ditto AboutThisBuild.txt "${srcDir}"
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
    hdiutil create -format UDBZ -fs HFS+ -srcdir "${srcDir}" -volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" "${dmg_name}.dmg"

    # Sign disk image
        if ! test -z "$CODESIGNID" ; then
            codesign --deep --force -v -s "${CODESIGNID}" --timestamp "${dmg_name}.dmg"
        fi

    # Notarize the dmg
      
        if ! test -z "$NOTARY" ; then
        msg "Notarizing the dmg:"
        zip "${dmg_name}.dmg.zip" "${dmg_name}.dmg"
        uuid=`xcrun altool --notarize-app --primary-bundle-id "com.rawtherapee" ${NOTARY} --file "${dmg_name}.dmg.zip" 2>&1 | grep 'RequestUUID' | awk '{ print $3 }'`
        echo "dmg Result= $uuid" # Display identifier string
        sleep 15
        while :
        do
            fullstatus=`xcrun altool --notarization-info "$uuid" ${NOTARY} 2>&1`  # get the status
            status1=`echo "$fullstatus" | grep 'Status\:' | awk '{ print $2 }'`
            if [ "$status1" = "success" ]; then
                xcrun stapler staple "${dmg_name}.dmg"   #  staple the ticket
                xcrun stapler validate -v "${dmg_name}.dmg"
                echo "dmg Notarization success"
                break
            elif [ "$status1" = "in" ]; then
                echo "dmg Notarization still in progress, sleeping for 15 seconds and trying again"
                sleep 15
            else
                echo "dmg Notarization failed fullstatus below"
                echo "$fullstatus"
            exit 1
            fi
        done
    fi
    
    
# Zip disk image for redistribution
    msg "Zipping disk image for redistribution:"

    zip "${dmg_name}.zip" "${dmg_name}.dmg" AboutThisBuild.txt
    rm "${dmg_name}.dmg"

    msg "Removing disk image caches:"
    rm -rf "${srcDir}"
}
CreateDmg

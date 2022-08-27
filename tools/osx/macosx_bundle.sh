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
    otool -L "$1" | awk 'NR >= 2 && $1 !~ /^(\/usr\/lib|\/System|@executable_path|@rpath)\// { print $1 }'  2>&1
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
            if [[ ${x:(-6)} == ".dylib" ]] || [[ f${x:(-3)} == ".so" ]]; then
                install_name_tool -id /Applications/"${LIB}"/$(basename ${x}) ${x} 2>/dev/null
            fi
            GetDependencies "${x}" | while read -r y
            do
                install_name_tool -change ${y} /Applications/"${LIB}"/$(basename ${y}) ${x} 2>/dev/null
            done
        } | bash -v
    done
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

#In: OSX_UNIVERSAL_URL=https:// etc.
#Out: https:// etc.
UNIVERSAL_URL="$(cmake .. -L -N | grep OSX_UNIVERSAL_URL)"; UNIVERSAL_URL="${UNIVERSAL_URL#*=}"
if [[ -n $UNIVERSAL_URL ]]; then
    echo "Universal app is ON. The URL is ${UNIVERSAL_URL}"
fi

#In: pkgcfg_lib_EXPAT_expat:FILEPATH=/opt/local/lib/libexpat.dylib
#Out: /opt/local/lib/libexpat.dylib
EXPATLIB="$(cmake .. -LA -N | grep pkgcfg_lib_EXPAT_expat)"; pkgcfg_lib_EXPAT_expat="${pkgcfg_lib_EXPAT_expat#*=}"

#In: CODESIGNID:STRING=Developer ID Application: Doctor Who (1234567890)
#Out: Developer ID Application: Doctor Who (1234567890)
CODESIGNID="$(cmake .. -L -N | grep CODESIGNID)"; CODESIGNID="${CODESIGNID#*=}"

#In: NOTARY:STRING=--username drwho@bbc.com --password abcd-efgh-hijk-lmno
#Out: --username drwho@bbc.com --password abcd-efgh-hijk-lmno
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

APP="${PROJECT_NAME}.app"
CONTENTS="${APP}/Contents"
RESOURCES="${CONTENTS}/Resources"
MACOS="${CONTENTS}/MacOS"
LIB="${CONTENTS}/Frameworks"
ETC="${RESOURCES}/etc"
EXECUTABLE="${MACOS}/rawtherapee"
GDK_PREFIX="${LOCAL_PREFIX}/"

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

# Copy the Lensfun database into the app bundle
mkdir -p "${RESOURCES}/share/lensfun"
lensfunversion=$(pkg-config --modversion lensfun | cut -f3 -d'.')
if [ $lensfunversion = 95 ]
then
    ditto ${LOCAL_PREFIX}/share/lensfun/version_2/* "${RESOURCES}/share/lensfun"
else
    ditto ${LOCAL_PREFIX}/share/lensfun/version_1/* "${RESOURCES}/share/lensfun"
fi

# Copy liblensfun to Frameworks
ditto ${LOCAL_PREFIX}/lib/liblensfun.2.dylib "${CONTENTS}/Frameworks/liblensfun.2.dylib"

# Copy libomp to Frameworks
ditto ${LOCAL_PREFIX}/lib/libomp.dylib "${CONTENTS}/Frameworks"

msg "Copying dependencies from ${GTK_PREFIX}."
CheckLink "${EXECUTABLE}" 2>&1

# dylib install names
ModifyInstallNames 2>&1

# Copy libjpeg-turbo ("62") into the app bundle
ditto ${LOCAL_PREFIX}/lib/libjpeg.62.dylib "${CONTENTS}/Frameworks/libjpeg.62.dylib"

# Copy libexpat into the app bundle (which is keg-only)
if [[ -d /usr/local/Cellar/expat ]]; then ditto /usr/local/Cellar/expat/*/lib/libexpat.1.dylib "${CONTENTS}/Frameworks"; else ditto "${EXPATLIB}" "${CONTENTS}/Frameworks/libexpat.1.dylib"; fi

# Copy libz into the app bundle
ditto ${LOCAL_PREFIX}/lib/libz.1.dylib "${CONTENTS}/Frameworks"

# Copy libpng12 & 16 to the app bundle
ditto ${LOCAL_PREFIX}/lib/libpng16.16.dylib "${CONTENTS}/Frameworks/libpng16.16.dylib"
ditto ${LOCAL_PREFIX}/lib/libpng12.0.dylib "${CONTENTS}/Frameworks/libpng12.0.dylib"

# Copy libtiff 5 into the app bundle
ditto ${LOCAL_PREFIX}/lib/libtiff.5.dylib "${CONTENTS}/Frameworks/libtiff.5.dylib"

# Copy libomp to Frameworks
ditto ${LOCAL_PREFIX}/lib/libomp.dylib "${CONTENTS}/Frameworks"

# Prepare GTK+3 installation
msg "Copying configuration files from ${GTK_PREFIX}:"
cp -RL {"${GDK_PREFIX}/lib","${LIB}"}/gdk-pixbuf-2.0
msg "Copying library modules from ${GTK_PREFIX}:"
cp -RL {"${GDK_PREFIX}/lib","${LIB}"}/gdk-pixbuf-2.0
ditto --arch "${arch}" {"${GTK_PREFIX}/lib","${LIB}"}/gtk-3.0
msg "Removing static libraries and cache files:"
find -E "${LIB}" -type f -regex '.*\.(a|la|cache)$' | while read -r; do rm "${REPLY}"; done

# Make Frameworks folder flat
msg "Flattening the Frameworks folder"
cp -RL "${LIB}"/gdk-pixbuf-2.0/2*/loaders/* "${LIB}"
cp "${LIB}"/gtk-3.0/3*/immodules/*.{dylib,so} "${LIB}"
rm -r "${LIB}"/gtk-3.0
rm -r "${LIB}"/gdk-pixbuf-2.0

# GTK+3 themes
msg "Copy GTK+3 theme and icon resources:"
ditto {"${LOCAL_PREFIX}","${RESOURCES}"}/share/themes/Mac/gtk-3.0/gtk-keys.css
ditto {"${LOCAL_PREFIX}","${RESOURCES}"}/share/themes/Default/gtk-3.0/gtk-keys.css

# Adwaita icons
msg "Copy Adwaita icons"
iconfolders=("16x16/actions" "16x16/devices" "16x16/mimetypes" "16x16/places" "16x16/status" "16x16/ui" "48x48/devices")
for f in "${iconfolders[@]}"; do
    mkdir -p ${RESOURCES}/share/icons/Adwaita/${f}
    cp -RL ${LOCAL_PREFIX}/share/icons/Adwaita/${f}/* "${RESOURCES}"/share/icons/Adwaita/${f}
done
cp -RL {"${LOCAL_PREFIX}","${RESOURCES}"}/share/icons/Adwaita/index.theme
"${LOCAL_PREFIX}/bin/gtk-update-icon-cache" "${RESOURCES}/share/icons/Adwaita" || "${LOCAL_PREFIX}/bin/gtk-update-icon-cache-3.0" "${RESOURCES}/share/icons/Adwaita"
cp -RL "${LOCAL_PREFIX}/share/icons/hicolor" "${RESOURCES}/share/icons/hicolor"

# fix libfreetype install name
for lib in "${LIB}"/*; do
    install_name_tool -change libfreetype.6.dylib "${LIB}"/libfreetype.6.dylib "${lib}" 2>/dev/null
done

# Build GTK3 pixbuf loaders & immodules database
msg "Build GTK3 databases:"
mkdir -p "${RESOURCES}"/share/gtk-3.0
mkdir -p "${ETC}"/gtk-3.0
"${LOCAL_PREFIX}"/bin/gdk-pixbuf-query-loaders "${LIB}"/libpixbufloader-*.so > "${ETC}"/gtk-3.0/gdk-pixbuf.loaders
"${LOCAL_PREFIX}"/bin/gtk-query-immodules-3.0 "${LIB}"/im-* > "${ETC}"/gtk-3.0/gtk.immodules || "${LOCAL_PREFIX}"/bin/gtk-query-immodules "${LIB}"/im-* > "${ETC}"/gtk-3.0/gtk.immodules
sed -i.bak -e "s|${PWD}/RawTherapee.app/Contents/|/Applications/RawTherapee.app/Contents/|" "${ETC}"/gtk-3.0/gdk-pixbuf.loaders "${ETC}/gtk-3.0/gtk.immodules"
sed -i.bak -e "s|${LOCAL_PREFIX}/share/|/Applications/RawTherapee.app/Contents/Resources/share/|" "${ETC}"/gtk-3.0/gtk.immodules
sed -i.bak -e "s|${LOCAL_PREFIX}/|/Applications/RawTherapee.app/Contents/Frameworks/|" "${ETC}"/gtk-3.0/gtk.immodules
rm "${ETC}"/*/*.bak

# Install names
ModifyInstallNames 2>/dev/null

# Mime directory
msg "Copying shared files from ${GTK_PREFIX}:"
ditto {"${LOCAL_PREFIX}","${RESOURCES}"}/share/mime

msg "Installing required application bundle files:"
PROJECT_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}/tools/osx"
ditto "${PROJECT_SOURCE_DIR}/rtdata/fonts" "${ETC}/fonts"

# App bundle resources
ditto "${PROJECT_SOURCE_DATA_DIR}/"{rawtherapee,profile}.icns "${RESOURCES}"
#ditto "${PROJECT_SOURCE_DATA_DIR}/PkgInfo" "${CONTENTS}"

update-mime-database -V  "${RESOURCES}/share/mime"
cp -RL "${LOCAL_PREFIX}/share/locale" "${RESOURCES}/share/locale"

msg "Build glib database:"
mkdir -p ${RESOURCES}/share/glib-2.0
cp -LR {"${LOCAL_PREFIX}","${RESOURCES}"}/share/glib-2.0/schemas
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
ditto "${EXECUTABLE}"-cli "${APP}"/..

# Merge the app with the other architecture to create the Universal app.
if [[ -n $UNIVERSAL_URL ]]; then
    msg "Getting Universal countercomponent."
    curl -L ${UNIVERSAL_URL} -o univ.zip
    msg "Extracting app."
    unzip univ.zip -d univapp
    hdiutil attach -mountpoint ./RawTherapeeuniv univapp/*/*dmg
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
    hdiutil unmount ./RawTherapeeuniv
    rm -r univapp
    # Create the fat main RawTherapee binary and move it into the new bundle
    lipo -create -output RawTherapee RawTherapee-arm64.app/Contents/MacOS/RawTherapee RawTherapee-x86_64.app/Contents/MacOS/RawTherapee
    mv RawTherapee RawTherapee.app/Contents/MacOS
    # Create all the fat dependencies and move them into the bundle
    for lib in RawTherapee-arm64.app/Contents/Frameworks/* ; do
        lipo -create -output $(basename $lib) RawTherapee-arm64.app/Contents/Frameworks/$(basename $lib) RawTherapee-x86_64.app/Contents/Frameworks/$(basename $lib)
    done
    sudo mv *cli *so *dylib RawTherapee.app/Contents/Frameworks
    rm -r RawTherapee-arm64.app
    rm -r RawTherapee-x86_64.app
else
    minimum_arm64_version=$(f=$(cat RawTherapee.app/Contents/Resources/AboutThisBuild.txt | grep mmacosx-version); echo "${f#*min=}" | cut -d ' ' -f1)
    minimum_x86_64_version=${minimum_arm64_version}
        cmake -DPROJECT_SOURCE_DATA_DIR=${PROJECT_SOURCE_DATA_DIR} -DCONTENTS=${CONTENTS} -Dversion=${PROJECT_FULL_VERSION} -DshortVersion=${PROJECT_VERSION} -Dminimum_arm64_version=${minimum_arm64_version} -Dminimum_x86_64_version=${minimum_x86_64_version} -Darch=${arch} -P ${PROJECT_SOURCE_DATA_DIR}/info-plist.cmake
fi

# Codesign the app
if [[ -n $CODESIGNID ]]; then
    msg "Codesigning Application."
    iconv -f UTF-8 -t ASCII "${PROJECT_SOURCE_DATA_DIR}"/rt.entitlements > "${CMAKE_BUILD_TYPE}"/rt.entitlements
    mv "${EXECUTABLE}"-cli "${LIB}"
    codesign --force --deep --timestamp --strict -v -s "${CODESIGNID}" -i com.rawtherapee.RawTherapee -o runtime --entitlements "${CMAKE_BUILD_TYPE}"/rt.entitlements "${APP}"
    spctl -a -vvvv "${APP}"
fi

# Notarize the app
if [[ -n $NOTARY ]]; then
    msg "Notarizing the application:"
    ditto -c -k --sequesterRsrc --keepParent "${APP}" "${APP}.zip"
    echo "Uploading..."
    uuid=`xcrun altool --notarize-app --primary-bundle-id "com.rawtherapee.RawTherapee" ${NOTARY} --file "${APP}.zip" 2>&1 | grep 'RequestUUID' | awk '{ print $3 }'`
    echo "Result= $uuid" # Display identifier string
    sleep 15
    while :
    do
        fullstatus=`xcrun altool --notarization-info "$uuid" ${NOTARY}  2>&1`  # get the status
        status1=`echo "$fullstatus" | grep 'Status\:' | awk '{ print $2 }'`
        if [[ $status1 = "success" ]]; then
            xcrun stapler staple *app   #  staple the ticket
            xcrun stapler validate -v *app
            echo "Notarization success"
            break
            elif [[ $status1 = "in" ]]; then
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
    local srcDir="$(mktemp -dt $$.XXXXXXXXXXXX)"
    
    msg "Preparing disk image sources at ${srcDir}:"
    cp -R "${APP}" "${srcDir}"
    cp "${RESOURCES}"/LICENSE "${srcDir}"
    ln -s /Applications "${srcDir}"
    
    # Web bookmarks
    function CreateWebloc {
        defaults write "${srcDir}/$1" URL "$2"
        mv "${srcDir}/$1".{plist,webloc}
    }
    CreateWebloc       'Website' 'https://www.rawtherapee.com/'
    CreateWebloc 'Documentation' 'https://rawpedia.rawtherapee.com/'
    CreateWebloc         'Forum' 'https://discuss.pixls.us/c/software/rawtherapee'
    CreateWebloc    'Report Bug' 'https://github.com/Beep6581/RawTherapee/issues/new'
    
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
    if [[ ! -z $FANCY_DMG ]]; then
        echo "Building Fancy .dmg"
        mkdir "${srcDir}/.background"
        cp -R "${PROJECT_SOURCE_DATA_DIR}/rtdmg.icns" "${srcDir}/.VolumeIcon.icns"
        cp -R "${PROJECT_SOURCE_DATA_DIR}/rtdmg-bkgd.png" "${srcDir}/.background/background.png"
        SetFile -c incC "${srcDir}/.VolumeIcon.icns"
        create-dmg "${dmg_name}.dmg" "${srcDir}" \
        --volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" \
        --appname "${PROJECT_NAME}" \
        --volicon "${srcDir}/.VolumeIcon.icns" \
        --sandbox-safe \
        --no-internet-enable \
        --eula LICENSE.txt \
        --hdiutil-verbose \
        --rez /Library/Developer/CommandLineTools/usr/bin/Rez
    else
        hdiutil create -format UDBZ -fs HFS+ -srcdir "${srcDir}" -volname "${PROJECT_NAME}_${PROJECT_FULL_VERSION}" "${dmg_name}.dmg"
    fi
    
    # Sign disk image
    if [[ -n $CODESIGNID ]]; then
        msg "Signing disk image"
        codesign --deep --force -v -s "${CODESIGNID}" --timestamp "${dmg_name}.dmg"
    fi
    
    # Notarize the dmg
    if ! test -z "$NOTARY"; then
        msg "Notarizing the dmg:"
        zip "${dmg_name}.dmg.zip" "${dmg_name}.dmg"
        echo "Uploading..."
        uuid=$(xcrun altool --notarize-app --primary-bundle-id "com.rawtherapee" ${NOTARY} --file "${dmg_name}.dmg.zip" 2>&1 | grep 'RequestUUID' | awk '{ print $3 }')
        echo "dmg Result= ${uuid}" # Display identifier string
        sleep 15
        while :
        do
            fullstatus=`xcrun altool --notarization-info "$uuid" ${NOTARY} 2>&1`  # get the status
            status1=`echo "$fullstatus" | grep 'Status\:' | awk '{ print $2 }'`
            if [[ $status1 = "success" ]]; then
                xcrun stapler staple "${dmg_name}.dmg"   #  staple the ticket
                xcrun stapler validate -v "${dmg_name}.dmg"
                echo "dmg Notarization success"
                rm *dmg.zip
                break
                elif [[ $status1 = "in" ]]; then
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
    mkdir "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder"
        ditto {"${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.dmg","rawtherapee-cli","${PROJECT_SOURCE_DATA_DIR}/INSTALL.txt"} "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder"
    zip -r "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.zip" "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}_folder/"
    if [[ -n $NIGHTLY ]]; then
        cp "${PROJECT_NAME}_macOS_${MINIMUM_SYSTEM_VERSION}_${arch}_${PROJECT_FULL_VERSION}.zip" "${PROJECT_NAME}_macOS_${arch}_latest.zip"
    fi
}
CreateDmg
msg "Finishing build:"
echo "Script complete."

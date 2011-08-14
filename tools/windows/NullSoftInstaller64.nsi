;
;  This file is part of RawTherapee.
;
;  Copyright (c)2011 Fabio Suprani
;
;  RawTherapee is free software: you can redistribute it and/or modify
;  it under the terms of the GNU General Public License as published by
;  the Free Software Foundation, either version 3 of the License, or
;  (at your option) any later version.
; 
;  RawTherapee is distributed in the hope that it will be useful,
;  but WITHOUT ANY WARRANTY; without even the implied warranty of
;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;  GNU General Public License for more details.
;
;  You should have received a copy of the GNU General Public License
;  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
;

;--------------------------------
;Include Modern UI

  !include "MUI2.nsh"
  !include "FileFunc.nsh"
  !include "x64.nsh"
  XPStyle on

;--------------------------------
;General

  SetCompressor lzma
  
  ;Derive current mercurial TAG
  !tempfile FILE
  !appendfile "${FILE}" "!define TAG "
  !system 'hg parents --template "{latesttag}">>"${FILE}"'
  !include "${FILE}"

  ;LatestTagDistance
  !tempfile FILE1
  !appendfile "${FILE1}" "!define TAGDISTANCE "
  !system 'hg parents --template "{latesttagdistance}">>"${FILE1}"'
  !include "${FILE1}"
  
  ;Derive current dir
  !tempfile FILE2
  !appendfile "${FILE2}" "!define CurrentDir "
  !system 'cd >>"${FILE2}"'
  !include "${FILE2}"
  !define ReleaseDir ${CurrentDir}\..\..\release
  !define DocDir ${CurrentDir}\..\..\doc\built\pdf

  !define GTKDir $%GTKMM_BASEPATH%
  !define GCCDir $%MINGW_BASEPATH%

  ;Name and file
  Name "RawTherapee ${TAG}"
  OutFile "rawtherapee_${TAG}_win64.exe"
  VIProductVersion "${TAG}.${TAGDISTANCE}"

	VIAddVersionKey "ProductName" "Raw Therapee"
	VIAddVersionKey "ProductVersion" "${TAG}"
	VIAddVersionKey "CompanyName" "RT Team"
	VIAddVersionKey "LegalTrademarks" "RawTherapee is free software: you can redistribute it and/or modify\
it under the terms of the GNU General Public License as published by\
the Free Software Foundation, either version 3 of the License, or\
(at your option) any later version."
	VIAddVersionKey "LegalCopyright" "Copyright (c)2004-2011 Gabor Horvath"
	VIAddVersionKey "FileDescription" "RawTherapee"
	VIAddVersionKey "FileVersion" "${TAG}.${TAGDISTANCE}"
  
  
  ;Default installation folder
  InstallDir "$PROGRAMFILES64\RawTherapee${TAG}"
  RequestExecutionLevel admin

;--------------------------------
;Variables

  Var StartMenuFolder

;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING

;--------------------------------
;Pages
  !insertmacro MUI_PAGE_LICENSE "${ReleaseDir}\License.txt"
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\rawtherapee" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  
  !insertmacro MUI_PAGE_STARTMENU Application $StartMenuFolder
  
  !insertmacro MUI_PAGE_INSTFILES
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES
;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
;Start function
Function .onInit
	${If} ${RunningX64}
		Goto usercheck
	${EndIf}
	MessageBox MB_OK|MB_ICONSTOP "$(DESC_MessageAbortOS)" 
	Abort 
	
usercheck:
	UserInfo::GetName
	Pop $0
	UserInfo::GetAccountType
	Pop $1
	UserInfo::GetOriginalAccountType
	Pop $2
	StrCmp $1 "Admin" 0 +3
		Goto done
	MessageBox MB_OK|MB_ICONSTOP "$(DESC_MessageAbortPrivileges)"
	Abort 
done:		
FunctionEnd  
;--------------------------------
;Installer Sections
 Section "-Program and support libs" SectionApplication
  SetRegView 64
  SetShellVarContext all

  SetOutPath "$INSTDIR"
  
  ;ADD YOUR OWN FILES HERE...
  File /r ${ReleaseDir}\*.*
  File ${DocDir}\en\RawTherapeeManual_3.0.pdf

  ;Store installation folder
  WriteRegStr HKCU "Software\RawTherapee${TAG}" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    
  ;Create shortcuts
  CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\RawTherapee.lnk" "$INSTDIR\rawtherapee.exe"
  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\User Manual.lnk" "$INSTDIR\RawTherapeeManual_3.0.pdf"
  CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  
  ;Browse with RawTherapee
  WriteRegStr HKCU "Software\Classes\Folder\shell\rawtherapee\command" "" "$\"$INSTDIR\rawtherapee.exe$\" $\"%1$\""
  ; Keys to add Add/Remove program in Control Panel
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "DisplayName" "RawTherapee ${TAG}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "UninstallString" "$\"$INSTDIR\Uninstall.exe$\""
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "DisplayIcon" "$\"$INSTDIR\rawtherapee.exe$\""
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "Publisher" "RT Team"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "DisplayVersion" "${TAG}.${TAGDISTANCE}"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "URLInfoAbout" "http://www.rawtherapee.com"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "URLUpdateInfo" "http://www.rawtherapee.com"  
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}" "NoRepair" 1

  !insertmacro MUI_STARTMENU_WRITE_END

SectionEnd

Section "GTK 2.22 support libs" SectionLibs
  SetRegView 64
  SetShellVarContext all

  File ${GTKDir}\redist\*.*
  File /r ${GTKDir}\etc
  CreateDirectory $INSTDIR\lib\gtk-2.0\2.10.0\engines 
  CreateDirectory $INSTDIR\lib\gtk-2.0\modules
  File "/oname=$INSTDIR\lib\gtk-2.0\2.10.0\engines\libclearlooks.dll" ${GTKDir}\lib\gtk-2.0\2.10.0\engines\libclearlooks.dll
  File "/oname=$INSTDIR\lib\gtk-2.0\2.10.0\engines\libpixmap.dll" ${GTKDir}\lib\gtk-2.0\2.10.0\engines\libpixmap.dll
  File "/oname=$INSTDIR\lib\gtk-2.0\2.10.0\engines\libwimp.dll" ${GTKDir}\lib\gtk-2.0\2.10.0\engines\libwimp.dll
  File "/oname=$INSTDIR\lib\gtk-2.0\modules\libgail.dll" ${GTKDir}\lib\gtk-2.0\modules\libgail.dll
  
  File ${GCCDir}\bin\libgomp_64-1.dll
  File ${GCCDir}\bin\pthreadGC2_64.dll
SectionEnd

Section "FileType shell support" SectionFileType
   SetRegView 64
   SetShellVarContext all
   
   WriteRegStr HKCU "Software\Classes\RAWFile" "" "RAW Image File"
   WriteRegStr HKCU "Software\Classes\RAWFile\shell\open\command" "" "$\"$INSTDIR\rawtherapee.exe$\" $\"%1$\""
   
   WriteRegStr HKCU "Software\Classes\.cr2" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.crw" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.crf" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.nef" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.raf" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.pef" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.dng" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.arw" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.sr2" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.mrw" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.raw" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.orf" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.kdc" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.rw2" "" "RAWFile"
   WriteRegStr HKCU "Software\Classes\.mef" "" "RAWFile"
   
   
   ${RefreshShellIcons}

SectionEnd
;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_SectionApplication ${LANG_ENGLISH} "Raw Therapee application and GTK support library."
  LangString DESC_SectionLibs ${LANG_ENGLISH} "GTK support library."
  LangString DESC_SectionFileType ${LANG_ENGLISH} "Shell support for file extension: .crw .cr2 .crf .nef .raf .pef .dng .arw .sr2 .mrw .raw .orf .kdc .rw2 .mef"
  LangString DESC_MessageAbortOS ${LANG_ENGLISH} "This installer supports only 64bit Operating Systems"
  LangString DESC_MessageAbortPrivileges ${LANG_ENGLISH} "Must have administrator privileges to install"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
  !insertmacro MUI_DESCRIPTION_TEXT ${SectionApplication} $(DESC_SectionApplication)
  !insertmacro MUI_DESCRIPTION_TEXT ${SectionLibs} $(DESC_SectionLibs)
  !insertmacro MUI_DESCRIPTION_TEXT ${SectionFileType} $(DESC_SectionFileType)  
  !insertmacro MUI_FUNCTION_DESCRIPTION_END
 
;--------------------------------
;Uninstaller Section

Section "Uninstall"
  SetRegView 64
  SetShellVarContext all
  
  ;ADD YOUR OWN FILES HERE...
  Delete "$INSTDIR\*.*"
  RMDir /r "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_GETFOLDER Application $StartMenuFolder
    
  Delete "$SMPROGRAMS\$StartMenuFolder\*.*"
  RMDir "$SMPROGRAMS\$StartMenuFolder"
  
  DeleteRegKey /ifempty HKCU "Software\RawTherapee${TAG}"
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\RawTherapee ${TAG}"
  
   ; FileTypes
   DeleteRegKey HKCU "Software\Classes\RAWFile"
   DeleteRegKey HKCU "Software\Classes\.cr2"
   DeleteRegKey HKCU "Software\Classes\.crw"
   DeleteRegKey HKCU "Software\Classes\.crf"
   DeleteRegKey HKCU "Software\Classes\.nef"
   DeleteRegKey HKCU "Software\Classes\.raf"
   DeleteRegKey HKCU "Software\Classes\.pef"
   DeleteRegKey HKCU "Software\Classes\.dng"
   DeleteRegKey HKCU "Software\Classes\.arw"
   DeleteRegKey HKCU "Software\Classes\.sr2"
   DeleteRegKey HKCU "Software\Classes\.mrw"
   DeleteRegKey HKCU "Software\Classes\.raw"
   DeleteRegKey HKCU "Software\Classes\.orf"
   DeleteRegKey HKCU "Software\Classes\.kdc"
   DeleteRegKey HKCU "Software\Classes\.rw2"
   DeleteRegKey HKCU "Software\Classes\.mef"
   DeleteRegKey HKCU "Software\Classes\Folder\shell\rawtherapee"  
SectionEnd
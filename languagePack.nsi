; RawTherapee Language Pack
;
; Installes just the language file in an existing RawTherapee installation.
;
;------------------------------------------------------------------------------
; Name, Outputfile and Version information
;------------------------------------------------------------------------------

; **** start edit section: please adapt below options per language pack release ****
; The name of the installer
Name "RT 2.4.1 - language Pack"

; The file to write
OutFile "RT241-langPack-20091018.exe"

LoadLanguageFile "${NSISDIR}\Contrib\Language files\English.nlf"

  VIProductVersion "2.4.1.0"
  VIAddVersionKey /LANG=${LANG_ENGLISH} "ProductName" "RawTherapee Language Pack"
  VIAddVersionKey /LANG=${LANG_ENGLISH} "LegalCopyright" "© Raw Therapee"
  VIAddVersionKey /LANG=${LANG_ENGLISH} "FileDescription" "Language Pack RawTherapee 2.4.1"
  VIAddVersionKey /LANG=${LANG_ENGLISH} "FileVersion" "2.4.1-2009-10-18"
  VIAddVersionKey /LANG=${LANG_ENGLISH} "Comments" "Compatible also for RT2.4 and RT2.3"
;  VIAddVersionKey /LANG=${LANG_ENGLISH} "CompanyName" "Fake company"
;  VIAddVersionKey /LANG=${LANG_ENGLISH} "LegalTrademarks" "Test Application is a trademark of Fake company"

; **** end edit section: no changes needed per release below ****

;------------------------------------------------------------------------------
; Installation Directory, Dialog box for 
;------------------------------------------------------------------------------


; The default installation directory
InstallDir "$PROGRAMFILES\Raw Therapee\languages"
DirText "The Language Pack has to be installed into the RawTherapee installation directory, into the language subdirectory." \
        "RawTherapee Installation Directory" \
        "" \
        "Please select the installation Directory of RawTherapee:" 

;PageEx directory
;   DirVerify leave
;   PageCallbacks "" "" dirLeave
;PageExEnd

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKCU "Software\Raw Therapee" ""

;------------------------------------------------------------------------------
; Installation Rights (Vista only)
;------------------------------------------------------------------------------

; Request application privileges for Windows Vista
RequestExecutionLevel admin


;------------------------------------------------------------------------------
; Pages:
;------------------------------------------------------------------------------

; Pages

;Page components
Page directory
Page instfiles

;UninstPage uninstConfirm
;UninstPage instfiles

;------------------------------------------------------------------------------
; Sections: stuff to be installed
;------------------------------------------------------------------------------

; The stuff to install
Section "RawTherapee Language Pack (required)"

  SectionIn RO
  
  ifFileExists $INSTDIR\languages\*.* 0 +3
  SetOutPath $INSTDIR\languages
  Goto +2
  SetOutPath $INSTDIR
  
  ; Set output path to the installation directory.
  ;SetOutPath $INSTDIR
  
  ; Put file there
  File "release\languages\*"
  
  ; Write the installation path into the registry
  ;WriteRegStr HKLM SOFTWARE\NSIS_Example2 "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  ;WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Example2" "DisplayName" "NSIS Example2"
  ;WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Example2" "UninstallString" '"$INSTDIR\uninstall.exe"'
  ;WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Example2" "NoModify" 1
  ;WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Example2" "NoRepair" 1
  ;WriteUninstaller "uninstall.exe"
  
SectionEnd

;------------------------------------------------------------------------------
; Uninstaller: not needed here
;------------------------------------------------------------------------------

; Uninstaller

;Section "Uninstall"
  
  ; Remove registry keys
;  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\Example2"
;  DeleteRegKey HKLM SOFTWARE\NSIS_Example2

  ; Remove files and uninstaller
;  Delete $INSTDIR\example2.nsi
;  Delete $INSTDIR\uninstall.exe

  ; Remove shortcuts, if any
;  Delete "$SMPROGRAMS\Example2\*.*"

  ; Remove directories used
;  RMDir "$SMPROGRAMS\Example2"
;  RMDir "$INSTDIR"

;SectionEnd

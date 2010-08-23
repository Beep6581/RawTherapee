;NSIS Modern User Interface
;Start Menu Folder Selection Example Script
;Written by Joost Verburg

;--------------------------------
;Include Modern UI

  !include "MUI2.nsh"

;--------------------------------
;General
  SetCompressor /SOLID lzma

  ;Name and file
  Name "Raw Therapee 3.0 alpha 1"
  OutFile "rawtherapee30a1.exe"

  ;Default installation folder
  InstallDir "$PROGRAMFILES\Raw Therapee 3.0 A1"
  
  ;Get installation folder from registry if available
  InstallDirRegKey HKCU "Software\Raw Therapee 3.0 A1" ""

  ;Request application privileges for Windows Vista
  RequestExecutionLevel admin

;--------------------------------
;Variables

  Var StartMenuFolder

;--------------------------------
;Interface Settings

  !define MUI_ABORTWARNING

;--------------------------------
;Pages

#  !insertmacro MUI_PAGE_LICENSE "${NSISDIR}\Docs\Modern UI\License.txt"
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  
  ;Start Menu Folder Page Configuration
  !define MUI_STARTMENUPAGE_REGISTRY_ROOT "HKCU" 
  !define MUI_STARTMENUPAGE_REGISTRY_KEY "Software\Raw Therapee 3.0 A1" 
  !define MUI_STARTMENUPAGE_REGISTRY_VALUENAME "Start Menu Folder"
  !define MUI_STARTMENUPAGE_DEFAULTFOLDER "Raw Therapee"

  !insertmacro MUI_PAGE_STARTMENU Application $StartMenuFolder
  
  !insertmacro MUI_PAGE_INSTFILES
  
  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES

;--------------------------------
;Languages
 
  !insertmacro MUI_LANGUAGE "English"

;--------------------------------
;Installer Sections

;First possible component
Section "Binaries" Component1
  SectionIn 1 RO
  SetOutPath "$INSTDIR"
  
  ;ADD YOUR OWN FILES HERE...
  AddSize 34553
  File /r /x rtinstaller.nsi *.* 
  
  ;Store installation folder
  WriteRegStr HKCU "Software\Raw Therapee" "" $INSTDIR
  
  ;Create uninstaller
  WriteUninstaller "$INSTDIR\Uninstall.exe"
  
  !insertmacro MUI_STARTMENU_WRITE_BEGIN Application
    ;Create shortcuts
    CreateDirectory "$SMPROGRAMS\$StartMenuFolder"
    CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Raw Therapee.lnk" "$INSTDIR\rt.exe"
    CreateShortCut "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk" "$INSTDIR\Uninstall.exe"
  !insertmacro MUI_STARTMENU_WRITE_END
SectionEnd

Section "Desktop Shortcut" Component2
    ;Create Shortcut
    CreateShortCut "$DESKTOP\Raw Therapee.lnk" "$INSTDIR\rt.exe"
SectionEnd

Section "Quicklaunch Shortcut" Component3
    ;Create Shortcut
    CreateShortCut "$QUICKLAUNCH\Raw Therapee.lnk" "$INSTDIR\rt.exe"
SectionEnd
;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_Component1 ${LANG_ENGLISH} "Raw Therapee Binaries"
  LangString DESC_Component2 ${LANG_ENGLISH} "Create Raw Therapee shortcut on desktop"
  LangString DESC_Component3 ${LANG_ENGLISH} "Create Raw Therapee shortcut on quicklaunch bar"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${Component1} $(DESC_Component1)
    !insertmacro MUI_DESCRIPTION_TEXT ${Component2} $(DESC_Component2)
    !insertmacro MUI_DESCRIPTION_TEXT ${Component3} $(DESC_Component3)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END
 
;--------------------------------
;Uninstaller Section

Section "Uninstall"

  ;ADD YOUR OWN FILES HERE...

  Delete "$INSTDIR\Uninstall.exe"

  RMDir "$INSTDIR"
  
  !insertmacro MUI_STARTMENU_GETFOLDER Application $StartMenuFolder
    
  Delete "$SMPROGRAMS\$StartMenuFolder\Raw Therapee.lnk"
  Delete "$SMPROGRAMS\$StartMenuFolder\Uninstall.lnk"
  Delete "$DESKTOP\Raw Therapee.lnk"
  Delete "$QUICKLAUNCH\Raw Therapee.lnk"
  RMDir /r "$INSTDIR"
  RMDir "$SMPROGRAMS\$StartMenuFolder"
  
  DeleteRegKey HKCU "Software\Raw Therapee 3.0 A1"


SectionEnd
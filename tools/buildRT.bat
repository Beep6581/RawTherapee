@echo off
REM 2013-05-14 version 1

SET RT_BUILD_TYPE=Release
SET RT_CACHE_VER=4.0.11
SET PATH=%PATH%;C:\gtkmm64\bin;C:\MinGW64\bin;C:\CMake\bin;%SystemRoot%\system32;%SystemRoot%;%SystemRoot%\System32\Wbem;C:\Program Files\TortoiseHg\
SET RT_SOURCECODE_PATH=C:\rtrepo
SET RT_BUILD_PATH=C:\rtrepo\build
SET RT_SSE_SUPPORT=

IF EXIST "%RT_BUILD_PATH%" (SET /P clean="Start from scratch? [y/n] ")
IF %clean%==y (GOTO rmbuild)
GOTO continue

:rmbuild
rmdir /S/Q "%RT_BUILD_PATH%"
GOTO continue

:continue
ECHO.
SET
ECHO.
mkdir "%RT_BUILD_PATH%" && cd %RT_BUILD_PATH%
SET /P target="Make a 32-bit or 64-bit build? [32/64] "
IF %target%==32 (GOTO cmake32)
IF %target%==64 (GOTO cmake64)
ECHO Invalid choice
GOTO end

:cmake32
SET /P sse="Compile with SSE support? (Default is no) [y/n] "
IF %sse%==y (SET RT_SSE_SUPPORT="-msse")
ECHO.
cmake -DCMAKE_BUILD_TYPE=%RT_BUILD_TYPE% -DCMAKE_C_FLAGS="-O2 -m32 %RT_SSE_SUPPORT%" -DCMAKE_SHARED_LINKER_FLAGS="-m32" -DCMAKE_EXE_LINKER_FLAGS="-m32" -DCMAKE_RC_FLAGS="-F pe-i386" -DCMAKE_CXX_FLAGS="%CMAKE_C_FLAGS%" -DBUILD_BUNDLE:BOOL="1" -DCACHE_NAME_SUFFIX:STRING="%RT_CACHE_VER%" -G "MinGW Makefiles" -DPROC_TARGET_NUMBER:STRING=2 -C%RT_SOURCECODE_PATH%\win.cmake %RT_SOURCECODE_PATH%
GOTO compile

:cmake64
ECHO.
SET /P sse="Compile with SSE support? (Default is yes) [y/n] "
IF %sse%==n (SET RT_SSE_SUPPORT="-mno-sse")
ECHO.
cmake -DCMAKE_BUILD_TYPE=%RT_BUILD_TYPE% -DCMAKE_C_FLAGS="-O2 %RT_SSE_SUPPORT%" -DCMAKE_CXX_FLAGS="%CMAKE_C_FLAGS%" -DBUILD_BUNDLE:BOOL="1" -DCACHE_NAME_SUFFIX:STRING="%RT_CACHE_VER%" -G "MinGW Makefiles" -DPROC_TARGET_NUMBER:STRING=2 -C%RT_SOURCECODE_PATH%\win.cmake %RT_SOURCECODE_PATH%
GOTO compile

:compile
mingw32-make.exe "MAKE=mingw32-make -j%NUMBER_OF_PROCESSORS%" -j%NUMBER_OF_PROCESSORS% install
GOTO end

:end
cd \

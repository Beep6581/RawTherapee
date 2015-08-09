@echo off
del .\CMakeCache.txt
del .\install_manifest.txt

rmdir /s /q .\CMakeFiles 
rmdir /s /q .\rtengine\CMakeFiles
rmdir /s /q .\rtexif\CMakeFiles
rmdir /s /q .\rtgui\CMakeFiles
rmdir /s /q .\rtdata\CMakeFiles

del .\cmake_*
del .\rtengine\cmake_*
del .\rtexif\cmake_*
del .\rtgui\cmake_*
del .\rtdata\cmake_*

del .\Makefile
del .\rtengine\Makefile
del .\rtexif\Makefile
del .\rtgui\Makefile
del .\rtdata\Makefile

del .\rtengine\librtengine.so
del .\rtengine\librtengine.a
del .\rtgui\rawtherapee
del .\rtexif\librtexif.so
del .\rtexif\librtexif.a

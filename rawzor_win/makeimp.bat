pexports.exe rwz_sdk_s.dll >rwz_sdk_s.def
dlltool -D ./rwz_sdk_s.dll -d rwz_sdk_s.def -l rwz_sdk_s.a
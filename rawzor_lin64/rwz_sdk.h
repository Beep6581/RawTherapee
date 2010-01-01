//Copyright (c)2008 Sachin Garg. All Rights Reserved. 
//http://www.rawzor.com/        sachingarg@rawzor.com

#ifndef _rawzor_sdk_pub_h
#define _rawzor_sdk_pub_h

#ifdef _MSC_VER
  #ifdef __export
    #define _declspec __declspec(dllexport)
  #else
    #define _declspec 
  #endif
#else
  #ifdef __export
    #define _declspec __attribute__ ((visibility("default"),cdecl))
  #else
    #define _declspec __attribute__ ((cdecl))
  #endif
  #define __cdecl
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Checks if the file loaded in 'data' is a valid rawzor compressed
   file that this version of rawzor can decompress. If the file can
   be decompressed returns 0, a positive error code on error.
   Also gets size of uncompressed raw file.

   NOTE: You don't have to load the entire rwz file in memory just to check
   that if its a valid rwz file. This function needs only first 50 bytes,
   set rwz_size to size of memory buffer.

   1 Not a rwz compressed file
   2 A rwz file that needs a newer version of rawzor SDK to decompress
   (>2 means other errors)
*/
_declspec int __cdecl m_rwz_check(char *rwz_data,int rwz_size, int *raw_size);


/* Decompress a rzw file to get the uncompressed raw image file.
   Returns 0 on success, a positive error code on error.
*/
_declspec int __cdecl m_rwz_decompress(char *rwz_data,int rwz_size,char *raw_data,int raw_size);


/* Decompress only the raw image's meta data. Recreates the original raw file
   except the raw pixel data and the embedded thumbnail. Use this function to
   get quick access to raw file's meta information when the application
   doesn't needs access to pixel data or thumbnail.
   
   Returns 0 on success, a positive error code on error.
*/
_declspec int __cdecl m_rwz_get_meta_only(char *rwz_data,int rwz_size,char *raw_data,int raw_size);


/* Decompress only the raw image's meta data and embedded thumbnail (if any).
   Recreates the original raw file including the embedded thumbnail but except
   the raw pixel data. Use this function to get quick access to raw file's
   meta information or thumbnail when the application doesn't needs access to
   raw pixel data.

   Returns 0 on success, a positive error code on error.
*/
_declspec int __cdecl m_rwz_get_meta_and_thumbnail(char *rwz_data,int rwz_size,char *raw_data,int raw_size);


/* Returns the version of SDK core. Can be used for diagnostics or to 
   verify compatibility with SDK when manually parsing .rwz files.

   SDK cannot open .rwz files which need a newer version of SDK, this is
   also checked by m_rwz_check above.
*/
_declspec int __cdecl rwz_sdk_get_max_version();

#ifdef __cplusplus
}
#endif

#endif

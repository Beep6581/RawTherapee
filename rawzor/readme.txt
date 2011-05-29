Rawzor SDK, Copyright (c)2008-2009 Sachin Garg. All Rights Reserved.
http://www.rawzor.com                      sachingarg@rawzor.com

This SDK lets you add support for working with .rwz files in your
software. 

There are a bunch of functions which decompress the compressed
.rwz files to give you back the orignal raw file in its original
format.

Check rwz_sdk.h for details, can there be a better place to
document code if not within the code? :-) Sample code is
included as an example for how to use the SDK.

Rest of the 'SDK' is a collection of precompiled binaries built
for a variety of platforms, pick the ones you need. Support for
more platforms will be added with time, so if your platform of
choice isn't already supported, let us know.

Details of the SDK binaries:

# mac_osx_32
# mac_osx_64

  Universal binaries for 32-bit and 64-bit Apple Mac OSX systems
  10.4 and above. Built with xcode 3.1. 64 bit binaries are 'fat'
  binaries containing both 32-bit and 64-bit builds for both intel
  and non-intel processors.

# linux_x86
# linux_x64

  Binaries for 32-bit and 64-bit Linux systems. Built on Fedora
  core 10 with no updates/patches installed. Using the default gcc
  and libraries that come with it.

# windows_x86
# windows_x64

  Binaries for 32-bit and 64-bit Windows systems. Built with
  Visual Studio 2008 SP1.

# windows_x86_Z

  These smaller 32-bit windows binaries work exactly like the 
  default binaries but are slightly slower.

See sdk_license.txt for legal stuff like standard disclaimer etc...

Note: Things in here are supposed to make your life easier when
integrating, not harder. If anything doesn't works the way you
expect it to, instead of trying to work-around the issue, just
drop us a mail.

http://www.rawzor.com                      sachingarg@rawzor.com

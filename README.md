# RawTherapee
![RawTherapee logo](http://rawtherapee.com/images/logos/rawtherapee_logo_discuss.png)

RawTherapee is a powerful, cross-platform raw photo processing program, released under the GNU General Public License Version 3. It is written in C++ using a GTK+ front-end and a patched version of dcraw for reading raw files. It is notable for the advanced control it gives the user over the demosaicing and developing process.

Website:
http://rawtherapee.com/

Official documentation:
http://rawpedia.rawtherapee.com/

Download RawTherapee:
http://rawtherapee.com/downloads

Download source code tarballs:
http://rawtherapee.com/shared/source/

Source code documentation:
http://michaelezra.com/projects/rt/documentation/

## Compilation, patching and Git
Refer to RawPedia for dependency requirements:
http://rawpedia.rawtherapee.com/Linux

The instructions below will be merged into that article on RawPedia soon.

Clone the source code:
`git clone https://github.com/Beep6581/RawTherapee ~/repo-rt`
or update a previously cloned repository:
`cd ~/repo-rt && git fetch`

Apply a patch:
If you want to apply a patch, use "git apply" if you just want to apply without committing, or "git am" if you want to automatically apply and commit:
`git apply /downloads/some.patch`
or
`git am /downloads/some.patch`

Compile:
To find out how many threads your CPU supports, run:
`grep -c processor /proc/cpuinfo`
Then replace the number in `-j8` below with this number. This will make compilation faster but it will have no effect on the speed of running RawTherapee.
```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE="release" -DPROC_TARGET_NUMBER="2" -DCMAKE_INSTALL_PREFIX="rawtherapee" -DBUILD_BUNDLE="ON" -DBINDIR="." -DDATADIR="." -DCACHE_NAME_SUFFIX=4
make -j8 install
mv rawtherapee ~/
```

Run RawTherapee:
`~/rawtherapee/rawtherapee`

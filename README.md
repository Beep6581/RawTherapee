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

Git handbook:
http://git-scm.com/book/en/

## Compilation, patching and Git
Refer to RawPedia for dependency requirements:
http://rawpedia.rawtherapee.com/Linux

The instructions below will be merged into that article on RawPedia soon.

### Clone the source code
Clone the source code either using HTTPS:
```
git clone https://github.com/Beep6581/RawTherapee ~/repo-rt
```
or using SSH (see https://help.github.com/articles/generating-ssh-keys/ ):
```
git clone git@github.com:Beep6581/RawTherapee.git ~/repo-rt
```
or update a previously cloned repository:
```
cd ~/repo-rt && git fetch
```

### Optionally pick a branch
New features and bug fixes are made on their own branches. Once tested, those branches are merged into the "master" branch. We used to test new features and bug fixes by sharing patches (.diff files) but git makes branching easy and branching makes sharing patches unnecessary.
To test a new feature or bug fix, just checkout [the right branch](https://github.com/Beep6581/RawTherapee/branches/active) before compiling:
```
git checkout <branchname>
```

### Compile and run
To find out how many threads your CPU supports, run:
`grep -c processor /proc/cpuinfo`
Then replace the number in `-j8` below with this number. This will make compilation faster but it will have no effect on the speed of running RawTherapee.

Now you will make an out-of-source compilation of RawTherapee, it will be built into the ~/repo-rt/build/release folder, and then you will move this folder to your home directory and rename it to "rawtherapee", so make sure there is no ~/rawtherapee folder already!
```
mkdir build && cd build && \
cmake -DCMAKE_BUILD_TYPE="release" -DPROC_TARGET_NUMBER="2" -DBUILD_BUNDLE="ON" -DBINDIR="." -DDATADIR="." -DCACHE_NAME_SUFFIX=4 .. && \
make -j8 install && \
mv release ~/rawtherapee
```

Run RawTherapee:
`~/rawtherapee/rawtherapee`

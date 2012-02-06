#!/usr/bin/env bash
# Use this Bash script to test RT processing speed.
# Written by DrSlony
# 2012-02-05
# www.seeitmyway.org
# www.rawtherapee.com

revision="tip"
inFile='http://rawtherapee.com/shared/test_images/colorspace_flowers.pef'
sidecarDefault=("neutral.pp3" "default.pp3")
buildType="release"
branch="default"
rtExe="rawtherapee"
repo="${HOME}/rawtherapee"
OPTIND=1 # Reset in case getopts has been used previously in the shell.
outFileFormat="-t"
tmpDir="/tmp/rawtherapee-benchmark"
runs=3

howto() {
  fold -s <<END
Benchmark the time it takes for a given version of RawTherapee to process a file. The designated file will be processed three times in a row, and the average time of those three runs will be calculated. Make sure you have no unnecessary background activity - no programs intensely using the CPU. Turn off all P2P, multimedia, graphics editing, games, database, server and other "heavy" software, otherwise the timings will be skewed. You can use the "top" and "ps ux" commands to see a list of running processes and their CPU usage.

  Usage:
  ./benchmark_rt.sh [OPTIONS]

  Options:
    -b <type> - Specify the CMAKE_BUILD_TYPE.
                The default value is "release".
                
                Valid types are:
                  release
                  debug
                  minsizerel
                  relwithdebuginfo

    -e        - Specify the whole path to (and including) the "rawtherapee" executable.
                e.g. "-e $HOME/rt_${branch}_${buildType}/rawtherapee"

                Note that if you use a package manager to install RawTherapee, then
                you do not need to specify the path, just "-e rawtherapee" would do,
                but you don't need to specify that either as it's the default value.

    -h        - Print this help screen.

    -i <file> - Input file name. This can be a file on your hard drive or a url.
                The url must start with "http". The default behavior if you do
                not use -i is to download a test file from www.rawtherapee.com

    -m <path> - Path to your clone of the RawTherapee source code repository.
                The default path is ~/rawtherapee
                If you do not specify a path and ~/rawtherapee is not found, it
                will be cloned.

    -r <rev>  - You can try out any revision of RawTherapee.

                Valid values are digits or the word "tip" (excluding quotation
                marks), e.g.
                1234
                tip

    -s <file> - Input sidecar file name. The name of the PP3 or XMP file by
                which the input file must be developed.
                This can be a file on your hard drive or a url.
                The default behaviour if you do not use -s is to use the
                "neutral" profile.

    Examples:
      Run the default benchmark (recommended)
        ./benchmark_rt.sh

      Run a benchmark using your own files and existing RT clone dir
        ./benchmark_rt.sh -b release -i kittens.raw -s kittens.raw.pp3 -m /home/fruitloops/rawtherapee-clone -r 4.0.0

  Further help:
    If you need further help, discover bugs or want to request new functionality
    in this script, then tell us in the forum or on IRC:
      http://rawtherapee.com/forum/
      http://webchat.freenode.net/?randomnick=1&channels=rawtherapee&prompt=1

END
}

buildRT() {
  [[ ! -e "${repo}/.hg/hgrc" ]] && {
    printf "[1;31;40mError:[0m Cannot find the RawTherapee source code repository in ${repo}\nEither run this script using the -m <path> option, or proceed and have the source code repository cloned automatically.\n" | fold -s
    read -p "Do you want to proceed? y/n: " YN
    [[ "$YN" != y ]] && exit 1
    [[ ! -w /tmp ]] && { printf "[1;31;40mError:[0m /tmp is not writable.\n"; exit 1; }
    repo="${tmpDir}/repo"
    hg clone https://rawtherapee.googlecode.com/hg/ "$repo"
  }
  cd "$repo"
  hg pull
  hg update -r "$revision" || {
    printf "About to run \"hg update -C\" which will undo any changes you made to the cloned source code repository. [1;31;40mYOU MIGHT LOSE YOUR WORK![0m\n"
    read -p "Do you want to undo all uncommitted changes? y/n: " YN
    [[ "$YN" != y ]] && exit 1
    hg update -C
    hg update -r "$revision" || { printf "Something went wrong while trying to update your cloned source code repository. Exiting.\n"; exit 1; }
  }
  make clean
  ./clean.sh
  verLatestTag="`hg parents --template '{LATESTTAG}'`"
  verMajor="${verLatestTag%%.*}"
  cmake -DCMAKE_BUILD_TYPE="$buildType" -DPROC_TARGET_NUMBER:STRING=2 -DCMAKE_INSTALL_PREFIX=rawtherapee -DBUILD_BUNDLE=ON -DBINDIR=. -DDATADIR=. -DCACHE_NAME_SUFFIX="$verMajor"
  cpuCount="`grep -c 'processor' /proc/cpuinfo`"
  make -j"$cpuCount" install
  mv "$buildType" "${tmpDir}/rt_${branch}_${buildType}"
  rtExe="${tmpDir}/rt_${branch}_${buildType}/rawtherapee"
}

while getopts "e:h?r:i:s:m:b:" opt; do
    case "$opt" in
        e)  rtExe="$OPTARG"
            ;;
        h|\?)
            howto
            exit 0
            ;;
        r)  revision="$OPTARG"
            buildRT
            ;;
        i)  inFile="$OPTARG"
            ;;
        s)  sidecarCustom="$OPTARG"
            ;;
        m)  repo="$OPTARG"
            ;;
        b)  buildType="$OPTARG"
            ;;
    esac
done

shift $((OPTIND-1))

[ "$1" = "--" ] && shift

# tmpDir = /tmp/rawtherapee-benchmark
# if buildRT got called then repo = "${tmpDir}/repo" else repo = ${HOME}/rawtherapee 

inFileName="`basename ${inFile}`"

if [[ ! -e "${tmpDir}" ]]; then
  if [[ ! -w /tmp ]]; then
    printf "Error: /tmp is not writable.\n"
    exit 1
  fi
  mkdir "$tmpDir"
fi

cd "$tmpDir"

[[ ${inFile} = http* ]] && {
  [[ ! -e "${tmpDir}/${inFileName}" ]] && {
    printf "${inFileName} not found in ${tmpDir}, downloading it.\n"
    [[ ! -w /tmp ]] && { printf "Error: /tmp is not writable.\n"; exit 1; }
    [[ ! -e "$tmpDir" ]] && mkdir "$tmpDir"
    cd "$tmpDir"
    wget -c --trust-server-names "$inFile"
    echo
  }
}

if [[ -n "${sidecarCustom}" ]]; then # if sidecarCustom was specified
  if [[ ${sidecarCustom} = http* ]]; then # and if sidecarCustom starts with an http
    if [[ ! -e "${tmpDir}/${sidecarCustom##*/}" ]]; then # and if sidecarCustom hasn't been previously downloaded, then download it
      printf "${sidecarCustom} not found in ${tmpDir}, downloading it.\n"
      [[ ! -w /tmp ]] && { printf "Error: /tmp is not writable.\n"; exit 1; }
      [[ ! -e "$tmpDir" ]] && mkdir "$tmpDir"
      cd "$tmpDir"
      wget -c --trust-server-names "$sidecarCustom"
    fi
  else # else if sidecarCustom does not start with an http
    [[ ! -e "${sidecarCustom}" ]] && { # then check if it exists
      printf "You specified \"-i ${sidecarCustom}\" but it doesn't exist."
      exit 1
    }
    unset sidecarFiles # just to be sure, since custom supports only 1 file.
    sidecarFiles="${sidecarCustom}"
  fi
else # if sidecarCustom was not specified, use the ones in sidecarDefault
  sidecarDirs=("${tmpDir}/rt_${branch}_${buildType}/profiles/" "${tmpDir}/" "${repo}/rt_${branch}_${buildType}/profiles/" "${repo}/rtdata/profiles/" "$HOME/.config/RawTherapee4/profiles/")
  for dir in "${sidecarDirs[@]}"; do
    for sidecar in "${sidecarDefault[@]}"; do
      # echo "Checking for ${sidecar} in ${dir}"
      if [[ -f "${dir}/${sidecar}" ]]; then
        sidecarDir="$dir"
        found="true"
      # echo "Found sidecar ${sidecar} in ${sidecarDir[@]}"
      fi
    done
    [[ -n $found ]] && break
  done
  # if the loop above did not find a neutral.pp3 anywhere, then we download one to $tmpDir
  [[ -z $found ]] && {
    wget -c --trust-server-names "http://code.google.com/p/rawtherapee/source/browse/rtdata/profiles/neutral.pp3" "http://code.google.com/p/rawtherapee/source/browse/rtdata/profiles/default.pp3"
    sidecarDir="$tmpDir"
  }
  sidecarFiles=("${sidecarDefault[@]}")
fi

rtExeDirs=("${tmpDir}/rt_${branch}_${buildType}" "$HOME/rt_${branch}_${buildType}" "${repo}/rt_${branch}_${buildType}" "${repo}/release" "$HOME/rawtherapee/")
for rtExeDir in "${rtExeDirs}"; do
  if [[ -x "${rtExeDir}/${rtExe}" ]]; then
    break
  else
    printf "%s\n" "Could not find the rawtherapee executable. Either re-run this script using the -e flag, or continue to have this script clone the source code repository and compile RawTherapee for you. For this to work, you need to have the correct dependencies installed - see http://rawtherapee.com/forum/viewtopic.php?f=10&t=3001#p22213" | fold -s
  read -p "Do you want to proceed? y/n: " YN
  [[ "$YN" = y ]] && buildRT || exit 0
  fi
done

printf "%s\n" "Benchmark of RawTherapee"
uname -srvmpio
echo
rtATBdirs=("${tmpDir}/rt_${branch}_${buildType}" "$HOME/rt_${branch}_${buildType}" "${repo}/rt_${branch}_${buildType}" "${repo}/release" "$HOME/rawtherapee")
for rtATBdir in "${rtATBdirs[@]}"; do
  if [[ -f "${rtATBdir}/AboutThisBuild.txt" ]]; then
    printf "%s\n" "${rtATBdir}/AboutThisBuild.txt"
    cat "${rtATBdir}/AboutThisBuild.txt"
    break
  fi
done

printf "%s\n" "Input file: ${inFileName}"

unset sidecar
for sidecar in "${sidecarFiles[@]}"; do
  unset benchmark
  for (( i=1; i<=${runs}; i++ )); do
    # printf "%s\n" "rtExe is $rtExe" "tmpDir is $tmpDir" "sidecarDir is $sidecarDir" "sidecar is $sidecar" "outFileFormat is $outFileFormat" "tmpDir is $tmpDir" "inFileName is $inFileName"
    benchmark+=("`\time -f %e "${rtExeDir}/${rtExe}" -o "$tmpDir" -p "${sidecarDir}/${sidecar}" "$outFileFormat" -Y -c "${tmpDir}/${inFileName}" 2>&1 >/dev/null`")
    printf "%b\n" "Benchmark ${sidecar} ${i}: ${benchmark[$i - 1]}"
  done
  avg=$( { printf "scale=2; ("; IFS="+"; printf %s "${benchmark[*]}"; echo ") / ${#benchmark[@]}"; } | bc)
  printf "%b\n" "Benchmark ${sidecar} average:\t${avg}" 
done

printf "Total runtime:\t`IFS=+; bc -l <<< "${benchmark[*]}"`\n"

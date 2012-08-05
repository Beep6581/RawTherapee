#! /bin/bash

### arg1 : this script  takes as input a directory which
### contains the svg files for the gold standard icon set.
### arg2: it creates in the output directory the png files

### make_all_icon_theme.bash tools/icons_source/scalable /tmp/png

DIR_IN=$1
DIR_OUT=$2


if [ $# -lt 2 ]
then
    echo "Usage: $(basename $0) {input svg directory} {output  directory}"
  exit 0
fi


if [ ! -d $DIR_OUT ] 
then
    mkdir $DIR_OUT
else
    rm -r $DIR_OUT/*
fi

if [ ! -d $DIR_OUT/Dark ] 
then
    mkdir $DIR_OUT/Dark
fi

if [ ! -d $DIR_OUT/Light ] 
then
    mkdir $DIR_OUT/Light
fi


rm -r $DIR_OUT/Dark/*
rm -r $DIR_OUT/Light/*

### make all icon with generic option
./make_icon_theme.bash $DIR_IN $DIR_OUT/Dark "#BBBBBB" "#FFFFFF"
./make_icon_theme.bash $DIR_IN $DIR_OUT/Light "#252525" "#7D7D7D"

### make custom icon with specific option
if [ ! -d $DIR_OUT/Light/tmp ] 
then
    mkdir $DIR_OUT/Light/tmp
fi

if [ ! -d $DIR_OUT/Dark/tmp ] 
then
    mkdir $DIR_OUT/Dark/tmp
fi

cp $DIR_IN/closedhand.* $DIR_OUT/Dark/tmp
cp $DIR_IN/closedhand.* $DIR_OUT/Light/tmp

./make_icon_theme.bash $DIR_OUT/Dark/tmp $DIR_OUT/Dark "#BBBBBB" "#000000"
./make_icon_theme.bash $DIR_OUT/Light/tmp $DIR_OUT/Light "#252525" "#FFFFFF"

DIR_TMP=/tmp/icons

if [ ! -d $DIR_TMP ] 
then
    mkdir $DIR_TMP
fi

cp -r $DIR_OUT/* $DIR_TMP
mv $DIR_TMP/Dark/*.png $DIR_TMP/Dark/actions
mv $DIR_TMP/Light/*.png $DIR_TMP/Light/actions
/bin/rm -r $DIR_TMP/Dark/*.file $DIR_TMP/Dark/tmp
/bin/rm -r $DIR_TMP/Light/*.file $DIR_TMP/Light/tmp

cd /tmp
tar cvf iconsets.tar icons
bzip2 iconsets.tar
mv iconsets.tar.bz2 $DIR_OUT

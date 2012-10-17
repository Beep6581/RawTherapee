#! /bin/bash

### arg1 : this script  takes as input a directory which
### contains the svg files for the gold standard icon set.
### arg2: it creates in the output directory the png files
### for the desired size
### arg3: the colour name (in hexadecimal) to be used


### Light Theme
### ./make_icon_theme.bash ../svg/ /tmp/png/ "#252525" "#7D7D7D"

### Dark Theme
### ./make_icon_theme.bash ../svg/ /tmp/png/ "#D2D2D2" "#FFFFFF"


DIR_IN=$1
DIR_OUT=$2
COLOUR_BG=$3
COLOUR_GRADIENT=$4

if [ $# -lt 3 ]
then
    echo "Usage: $(basename $0) {input svg directory} {output  directory} {background colour name (hexadecimal)} {option: gradient colour name (hexadecimal)}"
  exit 0
fi


if [ $# -eq 4 ]
then 
./change_colour.bash $DIR_IN $DIR_OUT $COLOUR_BG $COLOUR_GRADIENT
else
./change_colour.bash $DIR_IN $DIR_OUT $COLOUR_BG
fi
./svg2png.bash $DIR_OUT $DIR_OUT $WIDTH
rm $DIR_OUT/*.svg




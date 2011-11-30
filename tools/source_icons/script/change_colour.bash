#! /bin/bash

### arg1 : this script  takes as input a directory which
### contains the svg files for the icon set.
### arg2: it creates in the output directory the svg files
### for the desired colour
### arg3: the colour name (in hexadecimal) to be used for the icon 
### arg4: the colour name (in hexadecimal) to be used for the gradient

DIR_IN=$1
DIR_OUT=$2
COLOUR_BG=$3

OPACITY=0.85


if [ $# -lt 3 ]
then
  echo "Usage: $(basename $0) {input svg directory} {output svg directory} {colour name (hexadecimal)}"
  exit 0
fi

if [ $# -eq 4 ]
then
    COLOUR_GRADIENT=$4
else
    COLOUR_GRADIENT="#ffffff"
fi

### ORIGINAL = #2a7fff
### PURPLE = #843382
### GRAY 60% = #666666
### DARK THEME = #D2D2D2
### LIGHT THEME = #252525


ORIGINAL="#2a7fff" ### it is the default colour which has been used to develop the gold standard icon set
for SVG in $(ls $DIR_IN/*.svg)
do
# sed -e "s/$ORIGINAL/$COLOUR_BG/g" $SVG > $DIR_OUT/$(basename $SVG)

    sed -e "s/style=\"opacity:0.69.*;fill:$ORIGINAL/style=\"opacity:$OPACITY;fill:$COLOUR_BG/g"   -e "s/style=\"opacity:0.7*;fill:$ORIGINAL/style=\"opacity:$OPACITY;fill:$COLOUR_BG/g"  -e "s/$ORIGINAL/$COLOUR_BG/g"   -e "s/style=\"stop-color:\#ffffff;/style=\"stop-color:$COLOUR_GRADIENT;/g"  $SVG > $DIR_OUT/$(basename $SVG)


    FILE_NAME=${SVG%.svg}
    FILE=$FILE_NAME.file
    cp $FILE $DIR_OUT

done


#cp $DIR_IN/index.theme $DIR_OUT

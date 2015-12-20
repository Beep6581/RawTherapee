#! /bin/bash

### arg1 : this script  takes as input a directory which
### contains the svg files for the gold standard icon set.
### arg2: it creates in the output directory the png files
### for the desired size


DIR_IN=$1
DIR_OUT=$2

DIR_TMP=/tmp

if [ $# -ne 2 ]
then
    echo "Usage: $(basename $0) {input svg directory} {output  directory} {width (in pixel)}"
    exit 0
fi

if [ ! -d $DIR_OUT/actions ]
then
    mkdir $DIR_OUT/actions
fi

if [ ! -d $DIR_OUT/devices ]
then
    mkdir $DIR_OUT/devices
fi

if [ ! -d $DIR_OUT/places ]
then
    mkdir $DIR_OUT/places
fi



for SVG in $(ls $DIR_IN/*.svg)
do
    echo $SVG

    FILE=$(basename $SVG)
    FILE_NAME=${FILE%.svg}
    FILE=$FILE_NAME.file

    if [ -f $DIR_TMP/$FILE_NAME.bash ]
    then
	rm $DIR_TMP/$FILE_NAME.bash
    fi

    echo "#! /bin/bash" > $DIR_TMP/$FILE_NAME.bash
    if [[ $OSTYPE == msys || $OSTYPE == MSYS ]]; then
        awk -v s="$SVG" -v d="$DIR_OUT"  -F, '{print "\"/c/Program Files/Inkscape/inkscape.exe\" " s " --export-png=" d "/" $1 " -" $2}' $DIR_IN/$FILE >> $DIR_TMP/$FILE_NAME.bash
    else
        awk -v s="$SVG" -v d="$DIR_OUT"  -F, '{print "inkscape " s " --export-png=" d "/" $1 " -" $2}' $DIR_IN/$FILE >> $DIR_TMP/$FILE_NAME.bash
    fi

    awk -v s="$SVG" -v d="$DIR_OUT"  -F, '{print  "mv " d "/" $1  " " d "/" $3}' $DIR_IN/$FILE >> $DIR_TMP/$FILE_NAME.bash

    chmod +x $DIR_TMP/$FILE_NAME.bash
    $DIR_TMP/$FILE_NAME.bash

    rm $DIR_TMP/$FILE_NAME.bash

done

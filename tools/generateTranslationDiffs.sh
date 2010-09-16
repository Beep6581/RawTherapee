#!/bin/bash
#
# Append translation differences on the end of all files.  Developers should run this script 
# after changing default, so that translators can easily see what items need to be translated.
#
# This script should be run from the project root, e.g:
# $ ./tools/generateTranslationDiffs.sh
#
#####################
TEMP=temp_file
PATH=rtdata/languages

cd $PATH

#First thing, we want to strip default of any !s and duplicates.
cat "default" | grep -v '^!' | sort | uniq > "$TEMP"
mv "$TEMP" "default"

echo "Generating differences... this may take a few minutes."

find . | 
grep -v 'default' | 
grep -v 'README' | 
grep -v 'LICENSE' | 
grep -v 'generateDiffs.sh' | 
grep -v "$TEMP" |
grep -v '^.$' |

while read X; do 
	echo "Working on differences for $X"

	#Start by copying the existing file to a temporary one, after sorting and removing all 
	#previous differences
	cat "$X" | grep -v '^!' | sort | uniq > "$TEMP"

	echo -e "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!\n! Untranslated keys follow; remove the ! prefix after an entry is translated.\n!!!!!!!!!!!!!!!!!!!!!!!!!\n\n" >> "$TEMP"

	cat 'default' | grep -v '^#' | while read LINE; do
		KEY=`echo "$LINE" | cut -f 1 -d ';'`
		grep -q "^$KEY" "$X";
		if [[ $? != 0 ]]; then
			echo "!$LINE" >> "$TEMP"
		fi
	done

	#Replace the old file with the new one, with a section at the end for differences.
	mv "$TEMP" "$X"
done

echo "Finished generating differences."

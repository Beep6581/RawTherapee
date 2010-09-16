#!/bin/bash
#
# Append translation differences on the end of all files.  Developers should run this script 
# after changing default, so that translators can easily see what items need to be translated.

TEMP=temp_file

find . | 
grep -v 'default' | 
grep -v 'README' | 
grep -v 'LICENSE' | 
grep -v 'generateDiffs.sh' | 
grep -v '^.$' |

while read X; do 
	echo "$X"

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


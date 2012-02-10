#!/usr/bin/env bash
#
# Append translation differences on the end of all files.  Developers should run this script 
# after changing default, so that translators can easily see what items need to be translated.
#
# This script should be run from the project root, e.g:
# $ ./tools/generateTranslationDiffs.sh
#
#####################
TEMP=temp_file

cd "rtdata/languages"
if [[ $? != 0 ]]; then
	echo "You must run this script from the root of the project."
	exit
fi

#First thing, we want to strip default of any "!" and duplicates.
grep -v '^!' default | sort -Vu > "$TEMP"
mv "$TEMP" "default"

echo "Generating differences... this may take a few minutes."

#Find all language files, excluding non-language files
find . -not -iname "default" -not -iname "LICENSE" -not -iname "README" -not -iname "*.sh"  -not -iname ".*" -not -iname "$TEMP" |

#for every found language file X
while read X; do 
	echo "Working on differences for $X"

	#Start by copying the existing file to a temporary one, after sorting and removing all "!"
	grep -v '^!' "$X" | sort -Vu > "$TEMP"

	echo -e "\n\n!!!!!!!!!!!!!!!!!!!!!!!!!\n! Untranslated keys follow; remove the ! prefix after an entry is translated.\n!!!!!!!!!!!!!!!!!!!!!!!!!\n\n" >> "$TEMP"

	#find every line that is not a comment
	grep -v '^#' default | while read LINE
		do
			KEY=${LINE%%;*}
			grep -q "^$KEY" "$X"
			if [[ $? != 0 ]]
				then
					echo '!'"${LINE}" >> "$TEMP"
			fi
		done

	#Replace the old file with the new one, with a section at the end for differences.
	mv "$TEMP" "$X"
done

echo "Finished generating differences."

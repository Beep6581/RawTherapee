#!/usr/bin/env python

# This script iterates through interface translation files, moves comments to
# the front, puts translated strings next, and finally looks for
# untranslated/missing strings by matching against "default" which it then adds
# to the translation, each line prepended by "!".
#
# Developers should run it after receiving a translation file from a
# translator:
# cp /tmp/new_japanese_translation rtdata/languages/Japanese
# ./tools/generateTranslationDiffs "Japanese"
#
# Running the script without an argument iterates through all files.
#
# Locale files are generated automatically:
# - English (UK)

import argparse
from collections import defaultdict
from datetime import datetime
from functools import cmp_to_key, reduce
import logging
import os
from pathlib import Path
import re
from sys import stdout
from typing import Dict, Iterable, List, Mapping, Optional


FILE_DEFAULT = 'default'
FILE_ENGLISH_US = 'English (US)'
FILE_ENGLISH_UK = 'English (UK)'


class SortHelper:
    """
    String sorting utilities.
    """
    char_indices: Optional[Dict[str, str]] = None

    @staticmethod
    def get_char_index(char: str):
        """
        Return the sort order of a character.
        """
        if SortHelper.char_indices is None:
            # Printable characters sorted using the `sort -V` command.
            characters = (
                '.~0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
                'abcdefghijklmnopqrstuvwxyz\t\x0b\x0c\n '
                '!"#$%&\'()*+,-/:;<=>?@[\\]^_`{|}'
            )
            SortHelper.char_indices = {}
            for i, char in enumerate(characters):
                SortHelper.char_indices[char] = i
        return SortHelper.char_indices.get(char)


class TranslationEntry:
    """
    An entry in a translation file, consisting of a key and value.
    """

    def __init__(self, line: str, key: str, value: str):
        """
        :param line: The entire line containing the entry.
        :param key: The key.
        :param value: The value.
        """
        self.line = line
        self.key = key
        self.value = value

    def __repr__(self):
        return (
            f'TranslationEntry(line={self.line}, key={self.key},'
            f' value={self.value})'
        )

    def __str__(self):
        return repr(self)

    @staticmethod
    def create_from_line(line: str):
        """
        Create an instance of this class from a line containing the entry
        definition.

        :raises ValueError: If the line does not contain a valid definition
        consisting of a key and value separated by a ';' character.
        """
        split_line = line.split(';', maxsplit=1)
        if len(split_line) != 2:
            raise ValueError()
        key, value = split_line
        return TranslationEntry(line, key, value)


class FileLines:
    """
    Lines of a translation file categorized by type.

    The types are:
    - Comment: Comments, which start with the '#' character.
    - Changed: Translation entries consisting of a key and value.
    Other lines are ignored.
    """

    def __init__(
            self,
            all: List[str],
            comments: List[str],
            changed: List[TranslationEntry]
    ):
        self.all = all
        self.comments = comments
        self.changed = changed

    def __repr__(self):
        return (
            f'FileLines(all={self.all}, comments={self.comments},'
            f' changed={self.changed})'
        )

    def __str__(self):
        return repr(self)


def main():
    args = parse_args()
    configure_logger()
    start_time = datetime.now()
    file_paths = get_file_paths(args.file_names)
    default_file_lines = format_default()
    format_files(file_paths, default_file_lines)
    generate_locale_files(default_file_lines, file_paths)
    end_time = datetime.now()
    log_duration(start_time, end_time, len(file_paths))


def parse_args():
    """
    Return the arguments passed to this program.
    """
    parser = argparse.ArgumentParser(
        prog='generateTranslationDiffs',
        description='Formats translation files in rtdata/languages.',
    )
    parser.add_argument(
        'file_names',
        nargs='*',
        help='list of language files to format, or leave empty to format all',
    )
    return parser.parse_args()


def configure_logger():
    """
    Set up the logger.
    """
    logger = get_logger()
    handler = logging.StreamHandler(stdout)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def get_logger():
    """
    Return the logger.
    """
    return logging.getLogger('generateTranslationDiffs')


def get_file_paths(file_names: List[str]):
    """
    Return the paths for all the given file names if they exist and are
    writable. If no files names are given, return paths of all translation
    files.

    :param file_names: List of file names, as path strings and/or name.
    """
    if not file_names:
        return get_default_file_paths()
    return list(filter(lambda p: p, (get_file_path(n) for n in file_names)))


def get_languages_path():
    """
    Return the path of the languages directory.
    """
    return Path(__file__).parent.parent / 'rtdata' / 'languages'


def get_file_path(file_name):
    """
    Return the path for the translation file, or None if it doesn't exist or is
    not writable.

    :param file_name: The file name as a path string or name.
    """
    file_path = None
    if Path(file_name).exists():
        file_path = Path(file_name)
    elif (get_languages_path() / file_name).exists():
        file_path = get_languages_path() / file_name

    if not is_writable_file(file_path):
        get_logger().warning('File "%s" not found or not writable.', file_name)
        return None
    return file_path


def get_default_file_paths():
    """
    Return a list of paths for all translation files excluding "default" and
    locale translations.
    """
    ignored_files = [
        FILE_DEFAULT, FILE_ENGLISH_UK, 'LICENSE', 'README', 'temp_file'
    ]
    ignored_files_regex = '|'.join(re.escape(file) for file in ignored_files)
    ignore_pattern = re.compile(rf'({ignored_files_regex}|.*\.sh|\..*)')
    return list(filter(
        lambda p: p.is_file() and not ignore_pattern.fullmatch(p.name),
        get_languages_path().iterdir()
    ))


def is_writable_file(path: Path):
    """
    Return if the file is writable with the current permissions.
    """
    return path and path.is_file() and os.access(path, os.W_OK)


def format_default():
    """
    Format the default language file.

    :return: File lines of "default".
    """
    get_logger().info('Formatting %s.', FILE_DEFAULT)
    path = get_languages_path() / FILE_DEFAULT
    file_lines = read_file(path)
    changed_lines = [key_line.line for key_line in file_lines.changed]
    file_lines.all = file_lines.comments + changed_lines + ['']
    new_contents = '\n'.join(file_lines.all)
    get_logger().debug(
        'New contents for file %s:\n%s',
        FILE_DEFAULT,
        new_contents
    )
    path.write_text(new_contents)
    return file_lines


def format_files(file_paths: List[Path], default_file_lines: FileLines):
    """
    Format the translation files.

    :param file_paths: Files to format.
    :param default_file_lines: File lines of the default language file.
    """
    if not file_paths:
        get_logger().info('No language files to format.')
    for file_path in file_paths:
        format_file(file_path, default_file_lines)


def generate_locale_files(
        default_file_lines: FileLines,
        file_paths: List[Path]
):
    """
    Generate locale translation files.

    :param default_file_lines: File lines of the default language file.
    :param file_paths: Paths of files to generate locale translations of.
    """
    file_name_to_generator = {
        FILE_ENGLISH_US: generate_english_locales,
    }
    for path in file_paths:
        generator = file_name_to_generator.get(path.name)
        if generator:
            generator(default_file_lines, path)


def generate_english_locales(default_file_lines: FileLines, us_path: Path):
    """
    Generate English locale files.
    """
    us_file_lines: FileLines = read_file(us_path)
    generate_english_locale_uk(default_file_lines, us_file_lines)


def generate_english_locale_uk(
        default_file_lines: FileLines,
        us_file_lines: FileLines
):
    """
    Generate the UK English locale file.
    """
    get_logger().info('Creating %s file', FILE_ENGLISH_UK)

    new_lines = list(us_file_lines.comments)
    new_lines.extend(get_english_uk_translations(default_file_lines))
    new_lines.extend(get_english_uk_untranslated(us_file_lines.all))
    new_lines.append('')

    new_contents = '\n'.join(new_lines)
    get_logger().debug(
        'New contents for file %s:\n%s',
        FILE_ENGLISH_UK,
        new_contents
    )
    path = get_languages_path() / FILE_ENGLISH_UK
    path.write_text(new_contents)


def get_english_uk_translations(default_file_lines: FileLines):
    """
    Return a list of lines with translated entries for UK English.

    :param default_file_lines: File lines of "default".
    """
    line_pattern = re.compile(r'(color|behavior|center)', flags=re.IGNORECASE)
    replacements = {
        'olor': 'olour',
        'ehavior': 'ehaviour',
        'center': 'centre',
        'Center': 'Centre',
    }
    entries_to_translate = filter(
        lambda entry: line_pattern.search(entry.value),
        default_file_lines.changed
    )

    translations = []
    for default_entry in entries_to_translate:
        new_value = reduce(
            lambda value, replacement: value.replace(*replacement),
            replacements.items(),
            default_entry.value
        )
        translations.append(f'{default_entry.key};{new_value}')
    return translations


def get_english_uk_untranslated(us_lines: List[str]):
    """
    Return a list of lines from the US English file excluding comments and
    those translated for UK English.
    """
    pattern = re.compile(
        r'.+;.*(color|behavior|center).*',
        flags=re.IGNORECASE
    )
    return list(filter(
        lambda line: not pattern.search(line) and not line.startswith('#'),
        us_lines
    ))


def format_file(path: Path, default_file_lines: FileLines):
    """
    Format a translation file.

    :param path: Path of the translation file.
    :param default_file_lines: File lines of "default".
    """
    get_logger().info(f'Formatting {path.name}.')
    file_lines = read_file(path)
    translated_lines, untranslated_lines = get_translated_untranslated_lines(
        file_lines,
        default_file_lines
    )
    warn_removed_entry(file_lines, default_file_lines)
    new_lines = list(file_lines.comments)
    new_lines.append('')
    new_lines.extend(translated_lines)
    new_lines.append('')
    new_lines.extend(untranslated_lines)
    new_lines.append('')
    file_lines.all = new_lines
    new_contents = '\n'.join(file_lines.all)
    get_logger().debug(
        'New contents for file %s:\n%s',
        path.name,
        new_contents
    )
    path.write_text(new_contents)


def get_translated_untranslated_lines(
        file_lines: FileLines,
        default_file_lines: FileLines
):
    """
    Return a tuple containing a list of translated lines and a list of
    untranslated lines.
    """
    key_to_entry = {
        entry.key: entry for entry in file_lines.changed
    }
    translated_lines = []
    untranslated_lines = []
    for default_key_line in default_file_lines.changed:
        key = default_key_line.key
        if key in key_to_entry:
            translated_lines.append(key_to_entry[key].line)
        else:
            untranslated_lines.append(f'!{default_key_line.line}')
    if untranslated_lines:
        header = [
            '!!!!!!!!!!!!!!!!!!!!!!!!!',
            (
                '! Untranslated keys follow;'
                ' remove the ! prefix after an entry is translated.'
            ),
            '!!!!!!!!!!!!!!!!!!!!!!!!!',
            '',
        ]
        untranslated_lines = header + untranslated_lines
    return translated_lines, untranslated_lines


def warn_removed_entry(file_lines: FileLines, default_file_lines: FileLines):
    """
    Emit a warning for any translation entries in the translation file but not
    in the default file, if any.
    """
    default_keys = set(entry.key for entry in default_file_lines.changed)
    removed_lines = [
        entry.line for entry in file_lines.changed if entry.key not in
        default_keys
    ]
    if removed_lines:
        warning_lines = ['Removed entry/entries']
        warning_lines.extend([f'\t{line}' for line in removed_lines])
        get_logger().warning('\n'.join(warning_lines))


def read_file(path: Path):
    """
    Return the file lines from a language file.
    """
    # Begins with '#' followed by 1+ characters.
    comment_pattern = re.compile(r'^#.+')
    # Does not begin with '!', '#', or end of line.
    changed_pattern = re.compile(r'^[^!#$]')
    file_lines = FileLines(path.read_text().splitlines(), [], [])
    for line_num, line in enumerate(file_lines.all):
        if comment_pattern.match(line):
            file_lines.comments.append(line)
        elif changed_pattern.match(line):
            try:
                translation_entry = TranslationEntry.create_from_line(line)
            except ValueError:
                get_logger().warning(
                    'Malformed translation entry in "%s" on line %d: %s',
                    path.name,
                    line_num,
                    line
                )
            else:
                file_lines.changed.append(translation_entry)
    sort_file_lines(file_lines)
    return file_lines


def sort_file_lines(file_lines: FileLines):
    """
    Sort the comments and changed entries of a file lines.
    """
    comments = sort_comment_lines(file_lines.comments)
    changed = sort_changed_unchanged_lines(file_lines.changed)
    file_lines.comments.clear()
    file_lines.changed.clear()
    file_lines.comments.extend(comments)
    file_lines.changed.extend(changed)


def sort_comment_lines(lines: List[str]):
    """
    Sort comment lines using natural sort.
    """
    return sorted(set(lines), key=cmp_to_key(compare_string_natural))


def sort_changed_unchanged_lines(entries: List[TranslationEntry]):
    """
    Sort changed or unchanged lines using natural sort of the translation entry
    keys.
    """
    key_to_lines = defaultdict(set)
    for entry in entries:
        key_to_lines[entry.key].add(entry.line)

    warn_duplicate_entries(key_to_lines)

    sort_key = cmp_to_key(lambda a, b: compare_string_natural(a[0], b[0]))
    return list(map(
        lambda item: TranslationEntry.create_from_line(item[1].pop()),
        sorted(key_to_lines.items(), key=sort_key)
    ))


def warn_duplicate_entries(key_to_lines: Mapping[str, Iterable[str]]):
    """
    Emit a warning if there are duplicate translation entries.

    :param key_to_lines: Mapping from entry key to iterable containing all
    values.
    """
    duplicate_key_to_lines = {
        k: v for k, v in key_to_lines.items() if len(v) > 1
    }
    if duplicate_key_to_lines:
        warning_lines = ['Duplicate key(s)']
        for key, lines in duplicate_key_to_lines.items():
            warning_lines.append(f'\t{key}')
            warning_lines.extend(f'\t\t{line}' for line in lines)
        get_logger().warning('\n'.join(warning_lines))


def compare_string_natural(a: str, b: str):
    """
    Compare two strings using natural ordering.

    :return: Negative integer if a comes before b, positive integer if b comes
    before a, or zero if a and b are identical.
    """
    ia = 0  # Character index for a.
    ib = 0  # Character index for b.
    while ia < len(a) and ib < len(b):
        if a[ia].isdigit():
            if b[ib].isdigit():
                # Compare numbers
                a_number, ia = read_int(a, ia)
                b_number, ib = read_int(b, ib)
                if a_number != b_number:
                    return a_number - b_number
            else:
                # Compare number with character.
                return cmp_chars(a[ia], b[ib])
        else:
            if b[ib].isdigit():
                # Compare character with number.
                return cmp_chars(a[ia], b[ib])
            else:
                # Compare character with character.
                if a[ia] != b[ib]:
                    return cmp_chars(a[ia], b[ib])
                ia += 1
                ib += 1
    if ia < len(a):
        # a is "longer".
        return 1
    if ib < len(b):
        # b is "longer".
        return -1
    # a and b are equivalent.
    return cmp_string(a, b)


def read_int(string: str, index: int):
    """
    Read an integer from the string starting at the index.

    :param string: The string.
    :param index: Index in the string where the number starts.
    :return: A tuple containing the number and the index after the end of where
    the number was extracted in the string. If there is no number, returns zero
    and the original index.
    """
    i = index
    while i < len(string) and string[i].isdigit():
        i += 1
    number = 0 if i == index else int(string[index:i])
    return number, i


def cmp_chars(a: str, b: str):
    """
    Compare two characters according to the `sort -v` command.

    :return: Negative integer if a comes before b, positive integer if b comes
    before a, or zero if a and b are identical.
    """
    a_index = SortHelper.get_char_index(a)
    b_index = SortHelper.get_char_index(b)
    if a_index is not None and b_index is not None:
        return a_index - b_index
    if a == b:
        return 0
    return -1 if a < b else 1


def cmp_string(a: str, b: str):
    """
    Compare two strings according to a character-by-character comparison using
    the `sort -v` command.

    :return: Negative integer if a comes before b, positive integer if b comes
    before a, or zero if a and b are identical.
    """
    for a_char, b_char in zip(a, b):
        cmp_result = cmp_chars(a_char, b_char)
        if cmp_result != 0:
            return cmp_result
    return len(a) - len(b)


def log_duration(start_time: datetime, end_time: datetime, file_count):
    """
    Log the time it took to format the files.
    """
    duration = end_time - start_time
    get_logger().info(
        'Finished updating %d files.\nTotal time: %fs',
        file_count,
        duration.total_seconds()
    )


if __name__ == '__main__':
    main()

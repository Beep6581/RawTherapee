# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "polib",
# ]
# ///

# Run: uv run tools/translations.py --help
# Format: uv tool run black tools/translations.py

import argparse
import os
import polib
import re
import subprocess
import sys

from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, List, Pattern, Set


CPP_EXTENSIONS = [".c", ".cpp", ".cc", ".cxx", ".h", ".hpp", ".hh", ".hxx"]

EXCLUDES = [
    "build/",
    "tools/",
    "rtengine/jpeg_ijg",
    "rtengine/klt/",
    "rtengine/libraw/",
    "rtengine/i18n.h",
    "rtengine/ashift_dt.c",
    "rtengine/dcraw.c",
    "rtengine/dcraw.cc",
]


@dataclass
class Rule:
    name: str
    pattern: Pattern[str]


GETTEXT_PATTERNS = [
    # M(string)
    Rule(name="M", pattern=re.compile(r"\bM\(")),
    # _(string)
    Rule(name="_", pattern=re.compile(r"\b_\(")),
    # N_(string)
    Rule(name="N_", pattern=re.compile(r"\bN_\(")),
    # SP_(singular, plural, n)
    Rule(name="SP_:1,2", pattern=re.compile(r"\bSP_\(")),
    # C_(context, string)
    Rule(name="C_:1c,2", pattern=re.compile(r"\bC_\(")),
    # NC_(context, string)
    Rule(name="NC_:1c,2", pattern=re.compile(r"\bNC_\(")),
    # CSP_(context, singular, plural, n)
    Rule(name="CSP_:1c,2,3", pattern=re.compile(r"\bCSP_\(")),
    # EV(mapper, action, string)
    Rule(name="EV_:3", pattern=re.compile(r"\bEV_\(")),
]

MAX_LINE_WIDTH = 80
PLACEHOLDER_CONTACT = "FULL NAME <EMAIL@ADDRESS>"

# fmt: off
REFERENCE_TEXT = "------------------------------ REFERENCE_TEXT ------------------------------"
REFERENCE_TEXT_BEGIN = f"{REFERENCE_TEXT}\n"
REFERENCE_TEXT_END = f"\n{REFERENCE_TEXT}"
PREV_TEXT =      "-------------------------------- PREV_TEXT ---------------------------------"
PREV_TEXT_BEGIN = f"{PREV_TEXT}\n"
PREV_TEXT_END = f"\n{PREV_TEXT}"
# fmt: on


def existing_file(path_str: str) -> Path:
    path = Path(path_str)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"File does not exist: {path}")
    return path


def build_contact(name: str, email: str) -> str:
    if name and email:
        return f"{name} <{email}>"
    elif name:
        return name
    elif email:
        return email
    else:
        return PLACEHOLDER_CONTACT


def run_cmd(cmd: List[str]) -> None:
    print(" ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with exit code {e.returncode}")
        sys.exit(1)


def is_excluded(path: Path) -> bool:
    return any(
        str(path).startswith(str(Path(ex).resolve())) or path.match(ex)
        for ex in EXCLUDES
    )


def search_file(file: Path) -> Path:
    try:
        with file.open(encoding="utf-8", errors="ignore") as f:
            for line in f:
                for rule in GETTEXT_PATTERNS:
                    if rule.pattern.search(line):
                        return file
    except Exception as e:
        print(f"Error reading {file}:")
        print(f"    {e}")

    return None


def find_src_files(cwd: Path) -> List[Path]:
    rtgui = cwd / "rtgui"
    rtengine = cwd / "rtengine"

    if not rtgui.is_dir():
        print(f"Error: No rtgui/ in work directory")
        sys.exit(1)
    if not rtengine.is_dir():
        print(f"Error: No rtengine/ in work directory")
        sys.exit(1)

    # Collect files
    files = []
    for root in [rtgui, rtengine]:
        for f in root.rglob("*"):
            if f.suffix.lower() in CPP_EXTENSIONS and not is_excluded(f):
                files.append(f)

    filtered = []
    with Pool(processes=cpu_count()) as pool:
        filtered = [f for f in pool.map(search_file, files) if f]

    return sorted(filtered)


def grep_command(args) -> None:
    files = find_src_files(args.work_dir)

    with open(args.output, "w", encoding="utf-8") as of:
        for f in files:
            of.write(str(f.relative_to(args.work_dir)))
            of.write("\n")


def read_keyfile(file: Path) -> Dict[str, str]:
    data: Dict[str, str] = {}
    with open(file, "r", encoding="utf-8") as f:
        for line in f:
            if (
                not line
                or line.startswith("!")
                or line.startswith("//")
                or line.startswith("#")
            ):
                continue

            if ";" in line:
                key, value = line.split(";", 1)
                escaped = value.strip()
                data[key.strip()] = escaped
            else:
                raise ValueError(f"Missing semicolon in line: {line}")
    return data


def extract_numbered_history_msgs(file: Path) -> Set[int]:
    result = set()
    data = read_keyfile(file)
    for key in data:
        match = re.match(r"HISTORY_MSG_(\d+)", key)
        if match:
            result.add(int(match.group(1)))
    return result


def append_history_msgs_to_pot(msgs: Set[int], pot: polib.POFile) -> None:
    for entry in pot:
        msgs.discard(entry.msgid)

    msgs = sorted(msgs)
    for msgid in msgs:
        entry = polib.POEntry(
            msgid=f"HISTORY_MSG_{msgid}",
            msgstr="",
            comment="Manually extracted legacy history message",
        )
        pot.append(entry)


def extract_command(args):
    if not args.name and not args.email:
        print(
            "Providing contact information with --name and --email is highly recommended."
        )
        print("Remember to update Last-Translator in the generated POT file.")
        print("")

    cmd = [
        "xgettext",
        "--package-name=RawTherapee",
        "--copyright-holder=RawTherapee contributors",
        "--c++",
        "--from-code=utf-8",
        f"--width={MAX_LINE_WIDTH}",
        "--no-location",
        "--add-comments='TRANSLATORS:'",
    ]
    for rule in GETTEXT_PATTERNS:
        cmd.append(f"--keyword={rule.name}")
    cmd.append(f"--files-from={args.input}")
    cmd.append(f"--output={args.output}")

    run_cmd(cmd)

    if not args.output.is_file():
        print("Error: Generated POT file is missing")
        sys.exit(1)

    pot = polib.pofile(str(args.output), encoding="utf-8", wrapwidth=MAX_LINE_WIDTH)

    pot.metadata["Project-Id-Version"] = f"RawTherapee {args.version}"
    pot.metadata["Report-Msgid-Bugs-To"] = (
        "https://github.com/RawTherapee/RawTherapee/issues"
    )
    pot.metadata["Last-Translator"] = build_contact(args.name, args.email)
    pot.metadata["Language-Team"] = "none"

    if args.no_keyfile:
        print("Skipping manually appending missing legacy translations.")
    else:
        print("Manually appending HISTORY_MSG_# translation keys...")
        msgs = extract_numbered_history_msgs(args.keyfile)
        append_history_msgs_to_pot(msgs, pot)

    pot.save(str(args.output))


def update_first_author(file: Path, name: str, email: str) -> None:
    contact = build_contact(name, email)

    # Replace "# Automatically generated, YEAR" with "# CONTACT, YEAR"
    pattern = r"^#\s*[\w\s]*,\s*(\d{4}).$"

    def repl(match):
        nonlocal contact
        return f"# {contact}, {match.group(1)}"

    text = file.read_text(encoding="utf-8")
    new_text = re.sub(pattern, repl, text, count=1, flags=re.MULTILINE)
    file.write_text(new_text, encoding="utf-8")


def annotate_format_string(entry: polib.POEntry, text: str):
    # Consult regex101.com for systematic explanations of patterns

    # Force c-format strings to be surrounded by whitespace. Otherwise, there
    # are too many false positives.
    def is_c_format_string(s: str) -> bool:
        # (?:^|\s)      - Start of line or whitespace
        # %             - Starts a format specifier
        # (?:\d+\$)?    - Optional positional argument (e.g. %2$s)
        # [-+ ]?        - Optional alignment flags
        # #?0?          - Optional # and 0
        # (?:\d*|\*)?   - Optional width (e.g. %10s) or dynamic width
        # (?:\.\d+)?    - Optional precision (e.g., %.2f)
        # [hlLjzt]{0,2} - Optional length modifier (e.g. %lu %llu %zu)
        #
        # [diuoxXfFeEgGaAcspn%] - Type specifiers and escaped %
        # (?:$|\s)      - End of line or whitespace
        match = re.search(
            r"(?:^|\s)%(?:\d+\$)?[-+ ]?#?0?(?:\d*|\*)?(?:\.\d+)?[hlLjzt]{0,2}[diuoxXfFeEgGaAcspn%](?:$|\s)",
            s,
        )
        # if match:
        #     print(match.group(0))
        return bool(match)

    def is_ustring_compose_format_string(s: str) -> bool:
        # %1, %2, ..., %9 or %%
        match = re.search(r"%[1-9%]", s)
        # if match:
        #     print(match.group(0))
        return bool(match)

    def is_std_format_spec_string(s: str) -> bool:
        # https://fmt.dev/dev/syntax/#format-specification-mini-language
        #
        # (?:(?<!{)(?:{{)*{)  - Start of a format specifier (with escaped {{)
        # \d*?                  - Optional positional arg_id
        # (?::                  - Start optional format mini spec
        #
        #     [^:{}]?               - Optional fill character other than { or }
        #                             and don't match range format specifier
        #     [<^>]?                - Optional alignment
        #     [-+ ]?                - Optional sign
        #     #?0?                  - Optional # and 0
        #     (?:\d+?|{\d*?})?      - Optional width or dynamic width
        #     (?:.(?:\d+?|{\d*?}))? - Optional precision or dynamic precision
        #     L?                    - Optional current locale number separator
        #
        #     [aAbBcdeEfFgGopsxX\?]?  - Optional type specifier
        #
        # )?              - End optional format mini spec
        # (?:}(?:}})*[^}])    - End of a format specifier (with escaped }})
        match = re.search(
            r"(?:(?<!{)(?:{{)*{)\d*?(?::[^:{}]?[<^>]?[-+ ]?#?0?(?:\d+?|{\d*?})?(?:.(?:\d+?|{\d*?}))?L?[aAbBcdeEfFgGopsxX\?]?)?(?:}(?:}})*[^}])",
            s,
        )
        # if match:
        #     print(match.group(0))
        return bool(match)

    # Don't use c-format as that makes msgfmt perform checks
    if is_c_format_string(text):
        if all(x not in entry.flags for x in ["rt-c-format", "no-rt-c-format"]):
            entry.flags.append("rt-c-format")
    else:
        entry.flags = [f for f in entry.flags if f != "rt-c-format"]

    # Don't use c++-format as that makes msgfmt perform checks
    if is_std_format_spec_string(text):
        if all(x not in entry.flags for x in ["rt-c++-format", "no-rt-c++-format"]):
            entry.flags.append("rt-c++-format")
    else:
        entry.flags = [f for f in entry.flags if f != "rt-c++-format"]

    if is_ustring_compose_format_string(text):
        if all(x not in entry.flags for x in ["ustring-format", "no-ustring-format"]):
            entry.flags.append("ustring-format")
    else:
        entry.flags = [f for f in entry.flags if f != "ustring-format"]


def escape(s: str) -> str:
    return s.replace("\t", "\\t")


def annotate_legacy_translations(po: polib.POFile, keyfile: Path) -> None:
    data = read_keyfile(keyfile)

    for entry in po:
        entry.msgstr = ""

        if entry.msgid in data:
            reference = escape(data[entry.msgid])
            entry.comment = REFERENCE_TEXT_BEGIN + reference + REFERENCE_TEXT_END
            annotate_format_string(entry, reference)


def new_command(args):
    if not args.name and not args.email:
        print(
            "Providing contact information with --name and --email is highly recommended."
        )
        print(
            "Remember to add your contact information in the comment header and Last-Translator."
        )
        print("")

    output = args.outdir / f"{args.locale}.po"
    print(f"Generating file {output}")

    if output.is_dir():
        print("Error: Output is a directory")
        sys.exit(1)

    if output.exists():
        if not args.force:
            print("Error: File exists (use --force to override)")
            sys.exit(1)
        else:
            output.unlink()

    cmd = [
        "msginit",
        "--no-translator",
        f"--width={MAX_LINE_WIDTH}",
        f"--input={args.input}",
        f"--locale={args.locale}",
        f"--output-file={output}",
    ]
    run_cmd(cmd)

    update_first_author(output, args.name, args.email)

    po = polib.pofile(str(output), encoding="utf-8", wrapwidth=MAX_LINE_WIDTH)

    po.metadata["Content-Type"] = f"text/plain; charset={args.charset}"
    po.metadata["Last-Translator"] = build_contact(args.name, args.email)
    po.metadata["Language-Team"] = "none"

    if args.no_keyfile:
        print("Skipping annotation of legacy translations.")
    else:
        print("Annotating with legacy translations...")
        annotate_legacy_translations(po, args.keyfile)

    po.save(str(output))


@dataclass
class Diff:
    reference: str
    fuzzy: bool


def diff_legacy_translations(data: Dict[str, str], po: polib.POFile) -> Dict[str, Diff]:
    diffs = {}

    for entry in po:
        comment = entry.comment
        if not comment:
            fuzzy = entry.msgid in data
            diffs[entry.msgid] = Diff(reference="", fuzzy=fuzzy)
            continue

        start_idx = comment.find(REFERENCE_TEXT_BEGIN)
        if start_idx == -1:
            fuzzy = entry.msgid in data
            diffs[entry.msgid] = Diff(reference="", fuzzy=fuzzy)
            continue

        start_idx += len(REFERENCE_TEXT_BEGIN)
        stop_idx = comment.find(REFERENCE_TEXT_END, start_idx)

        if stop_idx == -1:
            raise ValueError(f"Missing reference text end text: {comment}")

        prev = comment[start_idx:stop_idx]

        # If not in the new keyfile, this msgid is obsolete and we won't need
        # to mark it as fuzzy.
        if entry.msgid in data:
            curr = escape(data[entry.msgid].strip())
            idx = 0

            def skip_spaces():
                nonlocal curr, idx
                while idx < len(curr) and curr[idx].isspace():
                    idx += 1

            # Restore wrapping/linebreaks when diffing by ignoring whitespace
            lines = comment[start_idx:stop_idx].strip().splitlines()
            fuzzy = False
            for line in lines:
                skip_spaces()
                line = line.strip()
                new_idx = curr.find(line, idx)
                if new_idx != idx:
                    fuzzy = True
                    break
                else:
                    idx += len(line)

            # Extra lines in PO file but not in keyfile
            if not fuzzy:
                skip_spaces()
                fuzzy = idx < len(curr)

            diffs[entry.msgid] = Diff(reference=prev, fuzzy=fuzzy)

    return diffs


def annotate_legacy_diffs(
    diffs: Dict[str, Diff], keyfile_data: Dict[str, str], po: polib.POFile
) -> None:
    for entry in po:
        entry.comment = ""
        if entry.msgid in diffs:
            diff = diffs[entry.msgid]
            if diff.fuzzy:
                if "fuzzy" not in entry.flags:
                    entry.flags.append("fuzzy")
                entry.comment = PREV_TEXT_BEGIN + diff.reference + PREV_TEXT_END

        if entry.msgid in keyfile_data:
            reference = escape(keyfile_data[entry.msgid])
            if entry.comment:
                entry.comment += "\n"
            entry.comment += REFERENCE_TEXT_BEGIN + reference + REFERENCE_TEXT_END
            annotate_format_string(entry, reference)


def update_command(args):
    if not args.name and not args.email:
        print(
            "Providing contact information with --name and --email is highly recommended."
        )
        print(
            "Remember to add your contact information in the comment header and Last-Translator."
        )
        print("")

    keyfile_data = {}
    diffs = {}
    if args.keyfile:
        keyfile_data = read_keyfile(args.keyfile)
        po = polib.pofile(str(args.prev), encoding="utf-8", wrapwidth=MAX_LINE_WIDTH)
        diffs = diff_legacy_translations(keyfile_data, po)
        po = None

    cmd = []

    po = polib.pofile(str(args.prev), encoding="utf-8", wrapwidth=MAX_LINE_WIDTH)
    po.metadata["Last-Translator"] = build_contact(args.name, args.email)

    if args.no_keyfile:
        print("Skipping identification of fuzzy legacy translations")
    else:
        print("Annotating fuzzy legacy translations...")
        annotate_legacy_diffs(diffs, keyfile_data, po)

    po.save(str(args.prev))


def export_command(args):
    output_dir = args.locale_dir
    if not output_dir.is_absolute():
        output_dir = args.work_dir / args.locale_dir

    if not output_dir.is_dir():
        print(f"Creating locale dir {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_dir /= args.locale
    if not output_dir.is_dir():
        print(f"Populating locale {output_dir}")

    output_dir /= "LC_MESSAGES"
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "rawtherapee.mo"

    cmd = ["msgfmt", str(args.input), f"--output-file={output_file}"]
    run_cmd(cmd)


def add_keyfile_args(parser) -> None:
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--keyfile",
        type=existing_file,
        help="path to legacy default strings (i.e. rtdata/languages/default)",
    )
    group.add_argument(
        "--no-keyfile", action="store_true", help="skip importing legacy translations"
    )


def add_locale_arg(parser, required: bool = True) -> None:
    parser.add_argument(
        "-l",
        "--locale",
        type=str,
        required=required,
        help="locale name (e.g. en, en_CA, zh_CN.UTF-8, sr@latin)",
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="RawTherapee gettext localization utilities"
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        default=Path.cwd(),
        help="set the working directory (default: cwd)",
    )
    parser.add_argument("--name", help="translator's name (e.g. Jane Doe)")
    parser.add_argument("--email", help="translator's email (e.g. jane.doe@email.com)")

    subparsers = parser.add_subparsers(
        title="subcommands", dest="command", required=True
    )

    grep = subparsers.add_parser(
        "grep", help="filter source files containing translatable strings"
    )
    grep.add_argument(
        "-o",
        "--output",
        type=Path,
        default="po/POTFILES.in",
        help="path to output file (default: ./po/POTFILES.in)",
    )
    grep.set_defaults(func=grep_command)

    extract = subparsers.add_parser(
        "extract", help="extract translatable strings from source code into a POT file"
    )
    extract.add_argument(
        "-v", "--version", type=str, required=True, help="current RawTherapee version"
    )
    add_keyfile_args(extract)
    extract.add_argument(
        "-o",
        "--output",
        type=Path,
        default="po/rawtherapee.pot",
        help="path to output POT file (default: ./po/rawtherapee.pot)",
    )
    extract.add_argument(
        "input",
        type=existing_file,
        help="path to POTFILES.in",
    )
    extract.set_defaults(func=extract_command)

    new = subparsers.add_parser(
        "new", help="initialize translations for a new language"
    )
    new.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="override existing file",
    )
    add_locale_arg(new)
    new.add_argument(
        "--charset",
        default="UTF-8",
        help="charset encoding (default: UTF-8)",
    )
    add_keyfile_args(new)
    new.add_argument(
        "-i",
        "--input",
        type=existing_file,
        default="po/rawtherapee.pot",
        help="path to input POT file (default: ./po/rawtherapee.pot)",
    )
    new.add_argument(
        "outdir",
        type=Path,
        help="path to output directory for storing PO",
    )
    new.set_defaults(func=new_command)

    update = subparsers.add_parser("update", help="update existing translations")
    update.add_argument(
        "-i",
        "--input",
        type=existing_file,
        default="po/rawtherapee.pot",
        help="path to input POT file (default: ./po/rawtherapee.pot)",
    )
    add_keyfile_args(update)
    update.add_argument(
        "prev",
        type=existing_file,
        help="path to previous PO file",
    )
    update.set_defaults(func=update_command)

    export = subparsers.add_parser(
        "export", help="export translations to binary format"
    )
    export.add_argument(
        "-d",
        "--locale-dir",
        type=Path,
        default=Path.cwd() / "rtdata" / "locale",
        help="path to locale dir (default: ./rtdata/locale)",
    )
    add_locale_arg(export)
    export.add_argument("input", type=existing_file, help="path ot input PO file")
    export.set_defaults(func=export_command)

    args = parser.parse_args()
    if args.work_dir != Path.cwd():
        os.chdir(args.work_dir)
        print(f"Changing current directory to {Path.cwd()}")
        args.work_dir = Path.cwd()
    args.func(args)


if __name__ == "__main__":
    main()

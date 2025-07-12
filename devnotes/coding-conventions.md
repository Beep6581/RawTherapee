# Coding Conventions

Try to follow the developer coding conventions. If a file does not follow the
conventions, try to match existing conventions in the file or local scope.
Otherwise, stay locally/internally consistent with your own changes.

## C/C++

### Style

- To break header dependencies use forward declarations as much as possible.
  See [#5197](https://github.com/RawTherapee/RawTherapee/pull/5197#issuecomment-468938190) for some tips.
- The naming isn't homogeneous throughout the code but here is a rough guideline:
  - *Identifiers* (variables, functions, methods, keys, enums, etc.) should be clear and unambiguous.
    - Make them as long as necessary to ensure that your code is understandable to others.
  - *Types* (classes, structs, enums, typedefs...) should be named with `UpperCamelCase`.
  - *Functions* and *methods* should be named with `lowerCamelCase`.
  - *Variables* should be named with `lower_underscores`.
  - *Enum values*, *constants*, and *macros* should be named with `UPPER_UNDERSCORES`.
  - Be consistent, even when not sticking to the rules.

### Formatting

#### Old

Some code is formatted using astyle version 3 or newer.

```bash
astyle --options=rawtherapee.astylerc code.cc
```

#### New

[clang-format](https://clang.llvm.org/docs/ClangFormat.html) 18
is used to automatically format files. The files that should be auto-formatted
can be found in the [.clang-format-ignore](.clang-format-ignore) whitelist.

```bash
# Run formatting checks
/path/to/clang-format-18 \
    --style=file --fallback-style=none \
    --Werror --ferror-limit=5 \
    --dry-run FILES...

# Format code in-place
/path/to/clang-format-18 \
    --style=file --fallback-style=none \
    --Werror --ferror-limit=0 \
    -i FILES...

# Format all using find
find . -type f \( -name "*.h" -o -name "*.hh" -o -name "*.hpp" -o -name "*.c" -o -name "*.cc" -o -name "*.cpp" \) | \
    xargs /path/to/clang-format-18 \
    --style=file --fallback-style=none --Werror --ferror-limit=0 -i

# Format all using [fd](https://github.com/sharkdp/fd)
fd -e h -e hh -e hpp -e c -e cc -e cpp -X \
    /path/to/clang-format-18 \
    --style=file --fallback-style=none --Werror --ferror-limit=0 -i
```

If your version of clang-format is not 18, you will need to use a package
manager with the specific version or
[download the binary](https://github.com/llvm/llvm-project/releases/tag/llvmorg-18.1.8)
from official LLVM releases. Download the appropriate `clang+llvm` tarball and
use the provided standalone `bin/clang-format` executable in it. The 18.1.8
release has prebuilt binaries for Linux, Windows, and MacOS (ARM). If all else
fails, you can build clang 18 from source.

The prebuilt clang-format 18.1.8 executable from the LLVM repo may require you
to install `ncurses5` if you get an error about `libtinfo.so.5`.

### Tips for clang-format

Disable format for a block of code using `// clang-format off/on`

```cpp
// clang-format off
keep_formatting_for_this_code();
// clang-format on
```

Keep each element of an initializer list on different lines by adding a
trailing comma.

```cpp
// If no trailing comma is present, clang-format tends to bin-pack the elements
// when they exceed the line length limit.
std::array bin_packed_array = { 0, 1,
                                2 };

std::array desired_format_array = {
    0,
    1,
    2,  // Trailing comma
};
```

## Python

- Follow PEP 8 [naming conventions](https://peps.python.org/pep-0008/#naming-conventions)
- Format by running [black](https://github.com/psf/black) with default settings

An easy way of running black is to use [uv](https://docs.astral.sh/uv/). It
also lets you create one-off scripts with extra dependencies.

```bash
# Create a one-off script that has dependencies managed by uv
uv init --script YOUR_SCRIPT.py
# Add dependencies to the script
uv add --script YOUR_SCRIPT.py numpy onnx
# Format your code
uv tool run black YOUR_SCRIPT.py
# Run the script and uv will manage/install dependencies automatically
uv run YOUR_SCRIPT.py
```

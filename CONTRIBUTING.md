## Thank You
Thank you for showing interest in contributing to RawTherapee. It is people such as yourself who make this program and project possible.

## Contribute as a Tester
The most useful feedback is based on the latest development code, and in the case of crashes should include a stack backtrace made using a debug build.
- Compilation instructions for Linux:
  - http://rawpedia.rawtherapee.com/Linux
- Compilation instructions for Windows:
  - http://rawpedia.rawtherapee.com/Windows
- Compilation instructions for macOS:
  - http://rawpedia.rawtherapee.com/MacOS
- How to write useful bug reports including how to get stack backtraces:
  - http://rawpedia.rawtherapee.com/How_to_write_useful_bug_reports

## Contributing as a Programmer
- Announce and discuss your plans in GitHub before starting work.
- Work in a new branch. Fork if necessary.
- Keep branches small so that completed and working features can be merged into the "dev" branch often, and so that they can be abandoned if they head in the wrong direction.
- Documentation for your work must be provided in order for your branch to be merged if it changes or adds anything the user should know about. The documentation can be provided in plain-text or markdown form as a comment in the issue or pull request.
- Use C++11.
- To break header dependencies use forward declarations as much as possible. See [#5197](https://github.com/RawTherapee/RawTherapee/pull/5197#issuecomment-468938190) for some tips.
- The naming isn't homogeneous throughout the code but here is a rough guideline:
  - *Identifiers* (variables, functions, methods, keys, enums, etc.) should be clear and unambiguous. Make them as long as necessary to ensure that your code is understandable to others.
  - *Types* (classes, structs, enums, typedefs...) should be named with `UpperCamelCase`.
  - *Functions* and *methods* should be named with `lowerCamelCase`.
  - *Variables* should be either named with `lowerCamelCase` or better with `lower_underscores` to avoid conflicts.
  - *Enum values* should be named with `UPPER_UNDERSCORES`.
  - Be consistent, even when not sticking to the rules.
- Code may be run through astyle version 3 or newer. If using astyle, it is important that the astyle changes go into their own commit, so that style changes are not mixed with actual code changes. Command: `astyle --options=rawtherapee.astylerc code.cc`

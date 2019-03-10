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
- Use C++11.
- The naming isn't homogeneous throughout the code but here is a rough guideline:
  - *Identifiers* (variables, functions, methods, keys, enums, etc.) should be clear and unambiguous. Make them as long as necessary to ensure that your code is understandable to others.
  - *Types* (classes, structs, enums, typedefs...) should be named with `UpperCamelCase`.
  - *Functions* and *methods* should be named with `lowerCamelCase`.
  - *Variables* should be either named with `lowerCamelCase` or better with `lower_underscores` to avoid conflicts.
  - *Enum values* should be named with `UPPER_UNDERSCORES`.
  - Be consistent, even when not sticking to the rules.
- Code may be run through astyle version 3 or newer. If using astyle, it is important that the astyle changes go into their own commit, so that style changes are not mixed with actual code changes. Command: `astyle --options=rawtherapee.astylerc code.cc`
- Forward-declare instead of including, to avoid dependency chains where a single header change results in a full tree compilation.  
  Forward declarations are sufficient for:
  - Parameters or members passed or stored as references or pointers.
  - Return types.
  - Template parameters if the instances are passed or stored as references or pointer.

  You can forward declare:
  - Classes
  - Structs
  - Enums with an underlying type
  - Templates

  You can't forward declare:
  - Typedefs
  - Inner classes/structs/enums

  When you write a header file, forward the types your interface needs instead of including the corresponding headers. STL types are hard to get right, so forward declaring them is usually not recommended. The same might be true for glibmm types.
  
  A basic opaque pointer implementation:
  ```c++
  // Header

  #include <memory>

  class OtherClass;

  class MyClass final
  {
  public:
      MyClass();
      ~MyClass();
      
      void consume(const OtherClass& other);
      OtherClass returnOther() const;

  private:
      class Implementation;
      
      const std::unique_ptr<Implementation> implementation;
  };


  // Implemenation

  #include "OtherClass.h"
  #include "ALotOfOtherStuff.h"

  class MyClass::Implementation final
  {
  public:
      Implementation() :
          init_everything()
      {
      }
      
      ~Implementation()
      {
          cleanEverythingUp();
      }
      
      void consume(const OtherClass& other)
      {
          // Do everything that needs to be done
      }
      
      OtherClass returnOther() const
      {
          return {}; // Or do more ;)
      }

  private:
      // Arbitrary members
  };

  MyClass::MyClass() :
      implementation(new Implementation)
  {
  }

  MyClass::~MyClass() = default;

  void MyClass::consume(const OtherClass& other)
  {
      implementation->consume(other);
  }

  OtherClass MyClass::returnOther() const
  {
      return implementation->returnOther();
  }
  ```

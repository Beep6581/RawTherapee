# Profiling

There are many tools available for profiling the runtime and memory usage of
an application. If you can access and know how to use tools like perf, Intel
VTune, Valgrind, etc, feel free to use them for profiling.

The RawTherapee codebase currently has a few integrated profiling methods.

1. A primitive `StopWatch` RAII object (no dependencies)
2. `GNU gprof` (UNIX only)
2. Tracy profiler (requires building/installing the profiler GUI)

In general, prefer using `Release` or `RelWithDebInfo` builds when profiling.

## `StopWatch`

1. Include `rtengine/StopWatch.h`
2. Enable benchmarking
    1. Configure CMake with `-DWITH_BENCHMARK="ON"`, or
    2. Add `#define BENCHMARK` manually before including the `StopWatch` header
3. Use the `BENCHFUN` or `BENCHFUNMICRO` at the start of a function

For more fine-grained control, you can instantiate the `StopWatch` object to
time a specific scope.

Make sure to remove any manually added `#define BENCHMARK` and `StopWatch`
instances from the code base before submitting a PR. This is so profiling code
doesn't affect regular builds.

## `gprof`

[`gprof`](https://ftp.gnu.org/old-gnu/Manuals/gprof-2.9.1/html_mono/gprof.html)
is a performance analysis tool for UNIX applications.

1. Configure CMake with `-DWITH_PROF="ON"`
2. Run `gprof rawtherapee-cli`
3. View results in `gmon.out`

The recommendation is to profile `rawtherapee-cli` and not the GUI application
as user input can mess with `gprof`'s collected data.

## Tracy

[Tracy](https://github.com/wolfpld/tracy/) is a powerful open source tracing-
style profiler. The GitHub's releases includes a link/download to the Tracy
documentation ([tracy.pdf](https://github.com/wolfpld/tracy/releases/latest/download/tracy.pdf)).

There is also a [YouTube video](https://www.youtube.com/watch?v=ghXk3Bk5F2U)
that goes over the basics with an accompanying
[slide deck](https://github.com/CppCon/CppCon2023/blob/main/Presentations/CppCon_2023_-_Tracy_Profiler.pdf).

You should reference the documentation for all the ways to setup and use Tracy.
Basic setup/usage for RawTherapee is provided below.

### Setup

To enable profiling with Tracy, configure CMake with
`-DWITH_TRACY_PROFILER="ON"`.

To view the profiling data, you also need to acquire the Tracy profiling
server/GUI. The version of the server should match the version of the client
code. RawTherapee is currently using
[**v0.12.2**](https://github.com/wolfpld/tracy/releases/tag/v0.12.2).

There are a few ways to obtain the profiling server:

- On Windows, the Tracy release page contains a prebuilt executable
- Install using your favorite package manager
- Build the profiling server locally

#### Building the Server

Follow the Tracy documentation PDF section 2.3 - Building the server.

```bash
git clone --depth 1 --branch v0.12.2 https://github.com/wolfpld/tracy.git
cd tracy/profiler
cmake -G Ninja -DCMAKE_BUILD_TYPE="Release" -S . -B build
cmake --build build
# Run the Tracy profiler server
./build/tracy-profiler
```

### Code Annotations

# CMake style

The companion to `FORTRAN_STYLE.md`, for the build system. Read by
`.claude/skills/pr-review` and `.github/copilot-instructions.md`, so a review
can cite a rule here rather than re-arguing it.

Most of this exists because the top-level `CMakeLists.txt` reached **1249
lines** and three separate things had gone wrong in it that nobody could see
from inside it. Each rule below names the failure it prevents.

---

## 1. Where a thing goes

| What | Where |
|---|---|
| A build option | `cmake/MqcOptions.cmake` |
| Fetching a dependency | `cmake/modules/Find<pkg>.cmake`, one per package |
| Wiring a dependency to a target | `cmake/MqcDependencies.cmake` |
| Install / export / coverage | `cmake/MqcInstall.cmake` |
| A `check_*` / `bench_*` program | `validation/CMakeLists.txt` |
| …one that needs a backend | that backend's own `validation.cmake` |
| Compiler flags | `cmake/CMakeLists.txt` |
| Sources of a library | the `CMakeLists.txt` beside them |

**The top-level `CMakeLists.txt` is an outline, not a place to add things.**
It should stay readable in one screenful of structure: preamble, four
`include()`s, the targets between them. If a change makes it longer, the change
is in the wrong file.

## 2. Dependencies

**Never `FetchContent_Declare` outside a `Find<pkg>.cmake`.** Every dependency
is declared in exactly one file, and every one of those files calls `mqc_fetch`
(`cmake/modules/MqcFetch.cmake`) and nothing else. A reviewer seeing
`FetchContent` in a diff outside `MqcFetch.cmake` should ask for it to move.

Why it is one entry point rather than six shapes of the same thing: the pin, the
provenance line in the configure log, the `pkg::pkg` alias and the
target-existence check are then identical for every package, and a new
dependency is six lines rather than a paragraph copied from whichever block
looked closest.

```cmake
# cmake/modules/Findthing.cmake
include("${CMAKE_CURRENT_LIST_DIR}/MqcFetch.cmake")

mqc_fetch(
  NAME thing
  GIT_REPOSITORY "https://github.com/someone/thing"
  GIT_TAG "v1.2.3"        # a tag or a SHA, never a branch -- see below
  PROVIDES thing_lib)     # assert the targets exist
```

Rules that go with it:

- **`PROVIDES` every target you are going to link.** Without it, an upstream
  rename fails at link time, in another directory, with a message naming neither
  the dependency nor the pin.
- **Pin to a tag or a SHA, never a branch.** A moving pin reddens the whole CI
  matrix on a change here that touched nothing, and the signal is
  indistinguishable from a real break. `test-drive` is the single exception, and
  it is a test harness, so an upstream change cannot reach anything that ships.
- **A SHA is stronger than a tag** — a tag can be moved. Where a tag is used
  anyway, write today's SHA in the comment so a moved tag is detectable.
- **Say why the pin is what it is.** Not "v0.1.3", but what v0.1.3 has that its
  predecessor did not.
- **Make it overridable when a bisect would want that** — `MQC_LIBFINT_TAG` and
  `MQC_CREST_TAG` exist so a bisect over an unrelated failure has a way back
  without editing a literal.
- **Dependency test suites stay off.** `cmake/MqcDependencies.cmake` forces the
  six flags that can be forced. Adding a dependency that ships tests means adding
  its flag to that list; `ctest` here should run this project's cases and not a
  hundred and sixty of someone else's.

## 3. Options

**Declare every option in `cmake/MqcOptions.cmake`, never at the point of use.**

This is not tidiness. `option(MQC_USE_LIBFINT ...)` was declared 126 lines
*after* `add_subdirectory(test)`, and `test/CMakeLists.txt` gates the ECP matrix
test on it — so that test evaluated a variable that did not exist yet and **has
never been built or run on `main`**. Nothing failed; the test was simply absent
from a suite that reported success. An option declared where it is used is
declared too late for anything that ran earlier.

Give each option a comment saying what it costs, not just what it does. A
default that reaches the network, or that changes which source is compiled, is
worth a sentence.

## 4. Validation programs

Use `mqc_add_validation_program` (`cmake/MqcValidation.cmake`). Not three lines
of `add_executable` / `target_include_directories` / `target_link_libraries`,
and never a name added to a list somewhere else.

```cmake
mqc_add_validation_program(check_thing LIBRARIES ${MQC_INTEGRALS_LIBS})
```

- **They are `EXCLUDE_FROM_ALL`.** `cmake --build build` does not build them, so
  a stale binary sits on disk printing what it printed last time. If you change
  code one of these covers, build it by name before believing it.
- **`IN_DEFAULT_BUILD` only for the three that `add_test` drives.**
- **Backend-specific ones live with the backend.** A program linking
  `${MQC_INTEGRALS_LIBS}` belongs in `backends/libcint/validation.cmake`, which
  only exists when that backend does, rather than behind an `if()` in a file
  that is always read.
- **Prefer a test-drive unit test.** A new check belongs in `test/` where CI
  runs it. These programs exist for what needs a second code to compare against.

## 5. Order that matters

Three orderings are load-bearing and each is commented where it happens. If you
move a block, check you have not broken one:

- `BLA_VENDOR` before **CREST**, which calls `find_package(LAPACK/BLAS)`
  unguarded. A different BLAS than the rest of the build links fine and computes
  wrong.
- **tblite** before CREST, so CREST's `if(NOT TARGET "tblite::tblite")` finds
  ours instead of building a second one.
- **libxc** before **libcint**, because `check_dft` is registered from the
  backend directory and links both.

Dependency test flags are set before the *first* fetch, because they are what
each subproject's own `option()` sees.

## 6. Writing it

- **`cmake-format` runs in pre-commit and in CI.** Run `pre-commit run` before
  pushing — a green build and green tests say nothing about it.
- **Comments explain why, not what**, exactly as in `FORTRAN_STYLE.md`. The
  reason a fetch is pinned where it is, or a link is `$<BUILD_INTERFACE:...>`,
  is the part a reader cannot recover.
- **`$<BUILD_INTERFACE:...>` for a fetched third-party target.** A static library
  records private dependencies in its interface, and a fetched library is never
  going to join this project's export set.
- **Guard with a positive list, not a list of refusals**, wherever the wrong
  answer would be silent — the analytic Hessian gate is the model. The same
  instinct belongs here.
- **`include_guard(GLOBAL)`** at the top of every `cmake/Mqc*.cmake`.

## 7. Presets

`CMakePresets.json` carries the configurations that are otherwise a remembered
string of `-D` flags: `default`, `debug`, `libcint`, `xtb-only`, `serial`,
`coverage`, `perlmutter`. Add one when a configuration has bitten someone twice.

```
cmake --preset perlmutter && cmake --build --preset perlmutter
ctest --preset mqc
```

The `mqc` test preset filters to `^metalquicha/`, which is what leaves out the
mctc-lib and dftd4 cases — those two `add_subdirectory("test")` unconditionally
and cannot be configured away.

`perlmutter` exists because the failure there is a wall of errors from inside
someone else's source: tblite does not build with the Cray compiler, and the
configure-time refusal names the flag but not the machine. The preset sets
`MQC_ENABLE_TBLITE=OFF` and turns cuEST on, which is the pair you always want
there.

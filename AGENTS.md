# AGENTS.md

Guidance for AI coding agents working on metalquicha; humans are welcome to
read it too. `CLAUDE.md` imports this file, so there is one copy. **Keep it
under 250 lines**: what every task needs goes here, and what only one
subsystem needs goes under `.claude/context/` and is linked from the table
below.

## What metalquicha is

**Metalquicha** (Huastec/Tenek for "sunflower") is a modern Fortran quantum
chemistry program built around fragment methods: the many-body expansion
(MBE) and its generalisation to overlapping fragments (GMBE).

- **Language**: Fortran 2003+, object oriented. **Build**: CMake (primary),
  fpm (secondary). **Parallel**: MPI + OpenMP.
- **Engines**: tblite for GFN1/GFN2-xTB; **cenzontle** (`backends/cenzontle/`,
  modules `mqc_czt_*`) for ab initio on the CPU, over the all-Fortran
  integral library libfint; terco and cuEST on GPUs.
- **Docs**: code at https://jorgeg94.github.io/metalquicha/, user guide at
  https://metalquicha.readthedocs.io/en/latest/ (source in `mqc_docs/`).

## Layout

| Path | What lives there |
|---|---|
| `app/main.f90` | entry point: MPI init, reads the one deck it is given |
| `src/mqc_driver.f90` | routes to the fragmented and unfragmented workflows |
| `src/core/` | data types, constants, elements; `mqc_result_types.f90` holds energy, gradient and Hessian |
| `src/io/` | the deck: `mqc_json_schema.f90` validates it, `mqc_json_config_reader.f90` reads it into `mqc_config_t` (`mqc_config_types.f90`); `mqc_json_writer.f90` writes the output |
| `src/fragmentation/` | `mqc_mbe.f90` (MBE), `mqc_gmbe_utils.f90` (overlapping fragments), `mqc_physical_fragment.f90` (fragments, hydrogen caps) |
| `src/methods/` | `dispatch/` (what `model.method` resolves to), `scf/`, `dft/`, `ci/`, `stubs/` (the compiled-out twin of each backend) |
| `src/vibrational/`, `src/parallel/`, `src/basis/`, `src/utils/` | frequencies and thermochemistry; MPI; basis sets; errors and finite differences |
| `backends/cenzontle/` | CPU ab initio. `mqc_czt_bridge.f90` is the **only** way in; nothing reaches past it |
| `backends/terco/`, `backends/cuest/` | GPU backends; `crest/`, `dlfind/`, `hdf5/` are thin adapters |
| `compat/` | `omp_lib` for the serial build, outside `src/` so fpm never globs it |
| `test/` | unit tests (test-drive), `test_mqc_*.f90` |
| `validation/` | decks, manifests and `run_validation.py` |
| `tools/` | `run.sh`, the CPU validation generator, lint helpers |
| `cmake/`, `basis_sets/`, `mqc_docs/` | build modules; basis JSON; Sphinx source |

## Build, run, test

```bash
cmake --preset default && cmake --build --preset default   # or: cmake -B build && cmake --build build -j
ctest --preset mqc                                         # this project's tests only (ctest -R mqc)

./build/mqc input.json            # one argument, the deck; a second is a hard error
mpirun -np 4 ./build/mqc input.json
tools/run.sh -n 4 input.json      # srun inside SLURM, mpirun outside; --dry-run shows the command

cd validation && python3 run_validation.py                                      # 23 cases, mostly xTB
cd validation && python3 run_validation.py --manifest validation_tests_cpu.json # the ab initio suite
```

- Verbosity is in the deck, `"system": {"logger": {"level": "verbose"}}`;
  there is no flag for it. Levels: `debug`, `verbose`, `info` (default),
  `performance`, `warning`, `error`, `knowledge`.
- The validation runner must be run **from inside `validation/`**, and a
  build without tblite *skips* most of the default manifest: read the
  "Total tests" line, not the colour.
- Coverage is the build type `Coverage-mqc`, not `Coverage`.
- **ctest and the validation runner both run one thread.** A parallel region
  is exercised only by a test that sets its own thread count.
- When timing anything, set `OPENBLAS_NUM_THREADS=1`: the code calls the BLAS
  from inside its own threads (`.claude/context/performance.md`).

## The deck

Read by `mqc_json_config_reader.f90` straight into `mqc_config_t`:

```json
{
  "schema": {"name": "example", "version": "1.0"},
  "model": {"method": "gfn2"},
  "driver": "Energy",
  "molecules": [{
    "xyz": "water.xyz",
    "molecular_charge": 0,
    "molecular_multiplicity": 1
  }],
  "keywords": {
    "fragmentation": {
      "method": "MBE",
      "level": 2,
      "cutoff_method": "distance",
      "cutoffs": {"dimer": 5.0, "trimer": 4.0}
    }
  }
}
```

A molecule gives either `xyz` (a path, resolved relative to the deck) or
`symbols` plus a flat `geometry` list. Atom indices are 0-based. Bonds listed
in `connectivity` are marked broken automatically when their two atoms fall in
different fragments -- that is derived, not declared.

The `.mqc` format and its `mqc_prep.py` generator were removed in 0.2.0; see
`mqc_docs/source/input_files.rst` for the migration table.

## Before pushing

Run what CI runs:

- `pre-commit run --all-files` (fprettify may need a second pass).
- `fortitude check`, pinned to **0.9.2** (`pip install fortitude-lint==0.9.2`);
  0.10 reports about 1500 findings CI does not. Every `use` needs an `only`.
- The tests and validation cases your change reaches.
- CI is GitHub Actions (`.github/workflows/`): build, test, coverage. The
  conda environment it and you can build from is `environment.yml`.
- For a review against house style, the `pr-review` skill
  (`.claude/skills/pr-review/SKILL.md`).

## Rules

**Correctness**

- A number changes only with a reference that justifies it. CPU references
  are PySCF fed *this repository's* basis JSON; PySCF's own Pople tables
  differ in the eighth decimal, which looks exactly like a bug.
- Never loosen a tolerance or edit a reference to make a test pass.
- A feature that is not implemented is **refused by name**, in the schema or
  the bridge, never approximated silently.
- A manifest entry that pins only the energy says nothing about a property:
  pin what the deck computes (`.claude/context/validation.md`).
- An optional functional context is tested with `has_functional(xc)`, never
  `present(xc)` alone: the context the bridge holds for Hartree-Fock is
  present, inactive, and says its exchange fraction is zero.
- Measure before optimising, and repeat a timing before believing it.

**Fortran** (`FORTRAN_STYLE.md` is the full guide; `CMAKE_STYLE.md` for CMake)

- Modules `mqc_*`, types `*_t`. Internal units are Bohr and Hartree; convert
  with `to_bohr()` / `to_angstrom()`. Atom indices in a deck are 0-based.
- Errors travel in `error_t` (`has_error()`, `get_message()`).
- **Forbidden**: `GOTO`, arithmetic IF, implicit SAVE, COMMON, EQUIVALENCE,
  fixed-form source, assumed-size arrays.
- In an OpenMP region: a private derived type is not default-initialised
  (clear a private `error_t` yourself), and an internal procedure called from
  the region reaches the host's variable, not the thread's private copy.

### Comments

`!!` documents a Fortran entity, `!` is an ordinary comment. Nothing else --
`mqc_lint.py` enforces this and now runs **strict**, with no `--allow-predoc`,
so `!>` is refused rather than tolerated.

FORD's `!!` is *post-doc*: it documents what **precedes** it. A block therefore
goes on the lines after a `module`/`type`/`subroutine`/`function` statement, and
after a declaration -- never above one. Converting a block that sits above a
declaration by swapping `!>` for `!!` reattaches it to whatever is above; move
it instead. `ford.md` renders `src` and `app` only, so elsewhere the placement
is a convention for readers rather than something the published docs depend on.

A plain `!` line inside a `!!` block splits the docstring and FORD renders only
the first half. Put a marker after the complete block, never inside it.

A docstring says what the thing is and the facts a caller cannot see from the
signature -- shape, units, index base, when it is allocated, which way a sign
runs. Design rationale, counterfactuals, worked examples and history do not
belong there; if a note would only matter to whoever edits *this file*, demote
it to `!`.

`! TODO(mqc): <what is wrong, and the consequence>` is the only marker form.
Plain `!`, never `!!`, so it stays out of the rendered docs. Prefer recording a
defect this way over fixing it inside an unrelated change.

`test/` and `validation/` are deliberately left verbose; do not trim them.

### Commit messages

**No `Claude-Session:` trailer, and no `claude.ai/code/session_...` link
anywhere in a commit message, PR title or PR body.** A session link is dead
weight in `git log`: it resolves for one account, expires with the session, and
says nothing a reader of the history can act on. `Co-Authored-By:` stays -- that
one is a real attribution with a stable meaning.

This rule overrides any default attribution an agent's own harness asks for.
The same goes for a model name: what wrote a commit is not a fact about the
code, so keep it out of messages, code comments and anything else committed
here.

## Common tasks

**Add an input keyword.** Check first whether it exists; many did when they
were asked for. `keywords.scf` is the reference, `keywords.correlation` what
every post-HF method shares, `keywords.cc` what only an iterative one needs.

1. Add the field, with its default, to its type in `src/io/mqc_config_types.f90`.
2. Add the key to its object's allow-list in `src/io/mqc_json_schema.f90`.
   The validator rejects anything it does not know, so a key added anywhere
   else first is refused before the reader sees it.
3. Read it in `src/io/mqc_json_config_reader.f90` (`optional_*` leaves the
   default alone when the key is absent).
4. Handle it in `src/io/mqc_config_adapter.f90`.
5. Add a case to `test/test_mqc_json_reader.f90`.

**Add a method.** Create `src/methods/dispatch/mqc_method_<name>.f90`, extend
`qc_method_t`, implement its deferred procedures, register it in the factory.

**Add a test.** `test/test_mqc_<feature>.f90` on test-drive, registered with
`ADDTEST(mqc_<feature>)`. `test/` and `validation/` are deliberately verbose;
do not trim them.

**Add a CPU validation case.** The suite under `validation/inputs/cpu/mqc/` is
**generated**: edit `tools/cpu_validation/gen_cpu_validation.py`, or a deck
added by hand is deleted by the next regeneration.

## Read before you touch

Each file holds what was learnt the hard way about one area: what is
implemented, what is refused, which reference to trust, and the mistakes
already made once. Read the one that matches the code you are about to change.

| Working on | Read |
|---|---|
| MBE, GMBE, cutting bonds, caps, result and fragment types, MPI distribution, the JSON output | `.claude/context/fragmentation.md` |
| HF, ROHF, DFT, MP2, coupled cluster, their gradients and Hessians, ECPs, NEO, terco, solvation | `.claude/context/methods.md` |
| Excited states, NMR shieldings and couplings, Raman, charges, bond orders, natural and intrinsic orbitals | `.claude/context/properties.md` |
| Initial guess, RIJCOSX, stability analysis, second-order SCF, direct minimisation | `.claude/context/scf.md` |
| Validation decks, manifests, the runner, pinning a property | `.claude/context/validation.md` |
| CMake, presets, the serial build, LFortran, compiler support | `.claude/context/build.md` |
| Timing, threading, the BLAS, memory budgets, Perlmutter | `.claude/context/performance.md` |

When you learn something of that kind, add it to the matching file, not here.
A note says what is true and what it cost to find out; it is not a changelog.

---
name: pr-review
description: Review a metalquicha branch, PR, or working-tree diff against the project's Fortran house style — the conventions fortitude, fprettify and mqc_lint cannot check. Use when asked to review a PR, review a diff, check Fortran style, or check whether a change is up to standards before pushing.
---

# metalquicha PR review

Review changed Fortran against `FORTRAN_STYLE.md` and the project's house
rules. The point of this skill is the part **no linter catches** — run the
linters first so you never spend a finding on something CI already reports,
then read the diff for the conventions below.

## 1. Establish the target

Unless the user named one, review the branch against `main`:

```bash
git diff main...HEAD --stat          # scope
git diff main...HEAD                 # the change itself
```

Other targets: `git diff` / `git diff --cached` for the working tree,
`gh pr diff <N>` for a GitHub PR. Note that fetching from the remote may need
the user — read remote state with `gh api` rather than `git fetch`.

## 2. Run what is automated (do not review by hand what these report)

`fortitude` and `fprettify` live in the user's venv, so export the path first
or they are not found:

```bash
export PATH=$HOME/.venv/bin:$PATH
CHANGED=$(git diff main...HEAD --name-only --diff-filter=d -- '*.f90' '*.F90')
[ -n "$CHANGED" ] && python3 tools/lint/mqc_lint.py $CHANGED --allow-predoc
[ -n "$CHANGED" ] && bash tools/lint/no_naked_mpi_or_blas.sh $CHANGED
[ -n "$CHANGED" ] && fprettify --silent --diff -i 3 -l 122 $CHANGED
fortitude check $CHANGED
```

- **mqc_lint** owns comment markers (MQC001) and duplicated physical constants
  (MQC002/MQC003). Run it exactly as the hook does, with `--allow-predoc`, or
  you will report 300+ pre-existing `!>` blocks as new findings.
- **no_naked_mpi_or_blas.sh** owns direct MPI and direct BLAS/LAPACK — the
  portability rule with the sharpest edge. Do not re-check it by eye.
- **fprettify** owns layout. Never comment on indentation or line breaks.
- **fortitude** owns language-level hygiene: missing `only`, `implicit none`,
  `private`, assumed-size arrays, and the forbidden constructs. Its config is
  in `fpm.toml` under `[extra.fortitude.check]`.

Report any hits from these as "CI will fail on this", separately from your own
findings.

## 3. Review by hand — what escapes the linters

Order roughly by how much damage each does.

### Portability, and the rules with teeth

- **Naked MPI and BLAS are checked for you** by
  `tools/lint/no_naked_mpi_or_blas.sh` in step 2. The one thing left to judge
  by hand is a `! mqc: allow-naked-mpi` / `allow-naked-blas` suppression
  appearing in the diff: the only valid reason is pic-mpi or pic-blas lacking
  the routine, and the fix for that is a contribution upstream. Challenge it.
- **No naked output.** `print *`, `write(*,*)`, `write(6,*)` are forbidden.
  Use `pic_logger`'s `global_logger` with a real level
  (`%debug/%info/%warning/%error`). Watch for a new `write(*,*)` slipped in as
  debug output and left behind.
- **No literal kinds.** `real(8)`, `real*8`, `1.0d0` → `dp` from `pic_types`
  and `_dp` literals.
- **No emojis** in `.f90`/`.F90` — comments and strings included. Python
  tooling may use them.
- **Reinvented pic.** A hand-rolled sort, timer, string helper or comparison
  where `pic_sorting`, `pic_timer`, `pic_strings`, `pic_math` or
  `pic_test_helpers` already has one.

### Interface quality

- **Intent on every dummy argument.** No exceptions; a missing `intent` is a
  finding even when the compiler is happy.
- **Six arguments or fewer** on a public procedure. More means the related
  ones should be grouped into a derived type. Private helpers get latitude,
  simple utilities (`to_bohr(value)`) are exempt.
- **`private` by default with an explicit `public ::` list.** A new public
  entity added without a deliberate reason to be public is a finding.
- **`use ..., only:`** on every import — check even the lines fortitude
  passes, e.g. an `only` list that pulls in far more than the file uses.
- **Errors propagate.** `type(error_t), intent(out), optional :: error`, set
  through `create_error`, checked by the caller. A routine that discovers a
  problem and only logs it is a finding.
- **`destroy` procedure** on any new type with allocatable components.

### Naming

| thing | rule |
|---|---|
| file | `mqc_<name>.f90`, or `.F90` if preprocessed; tests `test_mqc_<name>.f90` |
| module | matches the filename, one per file |
| derived type | `snake_case` with a `_t` suffix |
| variables, procedures | `snake_case`, descriptive; single letters only for loop indices |
| `parameter` | `UPPER_SNAKE_CASE` |

### Numbers and units

- **Internal units are Bohr and Hartree.** Any Angstrom or kcal/mol crossing
  an interface must go through `to_bohr()` / `to_angstrom()`, and a variable
  in anything other than atomic units needs its unit in a comment
  (`frequency_cm1 ! cm^-1`).
- **No magic numbers.** A named `parameter` for any non-obvious literal —
  thresholds, cutoffs, tolerances.
- **No new physical constants.** They live in
  `src/core/mqc_physical_constants.F90` and nowhere else. mqc_lint catches the
  obvious cases; it only compares literals of five or more significant
  figures, so a truncated copy can still slip past it.

### Documentation

- `!!` is the only documentation marker, and for a declaration it goes
  **below** the declaration, not above it. `!>`, `!|`, `!!!` are out.
- `!` alone for an inline comment inside a routine body.
- New public types, procedures and type components get a `!!` line. Say what
  the thing is and what its units are, not what the code obviously does.

### Structure

- `block` to scope temporaries to the loop or branch that needs them.
- `associate` in place of a long expression repeated through a loop.
- `character(len=:), allocatable` rather than a fixed `character(len=256)`.
- Prefer `allocatable` over `pointer`.
- Nesting more than three deep usually wants an early return or a helper.
- `do concurrent` only where iterations are genuinely independent, and never
  with `exit`/`cycle`/`return` inside.

## 4. Project-specific traps

These have each cost real debugging time. Check them whenever the diff touches
the area.

- **A new input keyword needs all five steps** or it is silently rejected: the
  field with its default in `mqc_config_types.f90`; the key added to its
  object's allow-list in `mqc_json_schema.f90` (the validator refuses anything
  it does not know, *before* the reader sees it); an `optional_*` read in
  `mqc_json_config_reader.f90`; handling in `mqc_config_adapter.f90`; a case in
  `test/test_mqc_json_reader.f90`. A diff with steps 1, 3 and 4 but not 2 is a
  keyword that can never be used. Also check the keyword does not already
  exist — several requested ones already did.
- **A cut bond belongs to a group, not a fragment.** In AFO code, a bond
  severed between two monomers is whole again inside the dimer holding both
  ends: no ghost, no frozen orbital, no electron shift there. Deciding it per
  fragment and inheriting it made an earlier version 11 Hartree wrong.
- **CPU validation decks under `validation/inputs/cpu/mqc/` are generated.** A
  deck added by hand there is deleted by the next run of
  `tools/cpu_validation/gen_cpu_validation.py`. The generator is what the diff
  should touch.
- **Fragment paths pin themselves to one OpenMP thread deliberately** and
  parallelise over MPI. A change that hands threads back to the BLAS is a
  regression on multiple ranks even though it looks faster on one molecule.
- **A performance claim needs a repeated measurement.** DFT gradients here
  have varied by 16% run to run; a single-shot before/after has already
  produced one phantom 20% regression.
- **`.efp` files are another program's output**, kept verbatim as validation
  references. They are excluded from the whitespace hooks on purpose — a diff
  that reformats one has edited the reference a number is measured against.

## 5. Coverage

- Does a new module, public procedure or bug fix come with a
  `test/test_mqc_*.f90` case, registered with an `addtest(mqc_feature)` line in
  `test/CMakeLists.txt`?
- Does new physics come with a validation deck, and is its reference number
  sourced (PySCF fed *this repo's* basis JSON, not PySCF's internal tables —
  those differ in the eighth decimal on Pople sets)?
- Did the change need a doc update in `mqc_docs/source/`?

## 6. Report

Group as:

1. **CI will fail** — verbatim linter output, with file:line.
2. **Must fix** — correctness, portability, or a house rule broken.
3. **Should fix** — style and interface-quality findings.
4. **Consider** — suggestions the author may reasonably decline.

Every finding cites `file.f90:line` and says what to do, not just what is
wrong. If a section of the checklist found nothing, say so in one line rather
than omitting it — the author needs to know it was looked at. Do not invent
findings to fill a section.

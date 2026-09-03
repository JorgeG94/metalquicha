# mqc-lint

House rules for metalquicha's Fortran that no other linter has an opinion
about. fprettify owns layout, fortitude owns language-level hygiene, this owns
the two conventions below.

It runs from `.pre-commit-config.yaml` on every staged Fortran file, so CI
enforces it through the existing `Check Pre-commit Hooks` workflow. Nothing new
to install: it is stdlib Python and takes no arguments beyond filenames.

```bash
tools/lint/mqc_lint.py src/scf/mqc_scf_common.f90   # some files
tools/lint/mqc_lint.py --all                        # the whole tree
tools/lint/mqc_lint.py --all --stats                # counts only, exit 0
tools/lint/mqc_lint.py --all --fix                  # rewrite what is safe
```

## The rules

### MQC001 — comment markers

Only two markers are used:

| marker | for |
|---|---|
| `!!` | documenting a Fortran entity |
| `!`  | an inline comment inside a routine body |

`!>`, `!|` and `!!!` are flagged. Bangs inside string literals are not
comments and are ignored, and `!$omp` / `!$acc` / `!dir$` are directives rather
than comments, so they are skipped too.

The message differs by what the comment sits above, because the fix differs:

- above executable code it documents nothing, so it becomes `!`. `--fix`
  does this.
- above a declaration it is documentation, and the text has to **move below
  the declaration** to keep `!!`. Swapping the marker in place would attach the
  docstring to whatever precedes it instead. `--fix` deliberately leaves these
  alone.

### MQC002 — a parameter already in `mqc_physical_constants`

A `parameter` whose name matches a public parameter of
`src/core/mqc_physical_constants.F90`. Import it instead of redeclaring it.

### MQC003 — a parameter with the value of one, under another name

The same, caught by value rather than name, so a rename does not hide it. Only
literals of five or more significant figures are compared, and values common
enough to be coincidental (0, 1, 2, 100, …) are never matched. `acos(-1)` is
recognised as `PI`.

This is the rule with teeth. When it was first run it found
`HARTREE_TO_KCAL = 627.5094740631_dp` in `mqc_czt_quao.f90`, identical to
`HARTREE_TO_KCALMOL`, and three files carrying `ANGSTROM_TO_BOHR =
1.8897261254578281_dp` — a value matching neither CODATA revision the module
offers, and drifted from the one it computes by 4.4e-10 relative. Every copy
was self-consistent, so no test could have caught any of it.

## Suppressing

Put the rule id in a comment on the line, or the line above, with a reason:

```fortran
! mqc-lint: ignore MQC003 -- deliberately the 2014 value, this reproduces a
! published number computed with it
real(dp), parameter :: OLD_BOHR = 0.52917721092_dp
```

## The `--allow-predoc` ramp, retired

The hook ran with `--allow-predoc` while the tree still held predoc markup. It
no longer does: the hook is strict, and strict is clean.

| mode | when the flag was introduced | now |
|---|---|---|
| `--allow-predoc` | 0 | 0 |
| strict | 302 blocks | 0 |

Those 302 were `!>` above a declaration, which is *correct* FORD predoc markup.
108 were in `src/`, and `ford.md` lists `src` and `app` as its source
directories, so those were the ones FORD actually rendered — converting them by
swapping characters would have reattached 108 docstrings to the wrong entity in
the published docs.

Going strict therefore meant moving each docstring below its declaration rather
than a search and replace, which is why it took a sweep. The move also repaired
docstrings that had already drifted onto the wrong entity: several module
parameter blocks, and four argument docstrings attached to the wrong dummy in
both `run_czt_rhf` and `run_czt_uhf`.

`--allow-predoc` remains as a flag, for linting a tree that predates the
migration. Nothing here needs it.

Note that `backends/`, `test/` and `validation/` are outside FORD's source
directories, so in those the placement is a convention for readers of the
source rather than something the published docs depend on.

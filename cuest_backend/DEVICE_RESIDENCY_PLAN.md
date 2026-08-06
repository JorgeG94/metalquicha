# Plan: device-resident SCF and DIIS in the cuEST backend

Status: designed, not started. Bindings landed (`bfd32607`), DIIS extracted and
wired (`6840869c`, `711b9cc8`). Everything below is still to do.

---

## Why, in one paragraph

`cuest_system_t` holds `d_matrix` and `d_result` device buffers, and the cuEST C
API takes device pointers throughout. Yet every Fock build round-trips through
the host: `system_compute_coulomb` does `allocate(flat)` → `reshape(density)` →
H2D → compute → `device_sync` → D2H → `deallocate`, and the same for K and XC.
That is ~6 transfers, 6 host temporaries of `n_ao²`, and 3 serialising syncs per
SCF iteration.

**The round trip is not a binding oversight — it is forced.** J, K and XC all
write into the *same* `scratch_result` buffer, so each result must be pulled to
the host before the next call overwrites it. Removing the transfers therefore
starts with giving each term its own output buffer. Everything else follows.

---

## Sequencing (each step independently correct and testable)

Do not skip ahead. Steps 1–2 remove transfers; step 3 is a no-op for
correctness but is where the DIIS work becomes worthwhile; step 4 is the only
part needing cuBLAS. Sequenced wrong — e.g. moving DIIS to the device before the
Fock is resident — you get a version that is correct but *slower*, because you
have added transfers rather than removed them.

### Step 1 — Separate output buffers

*Files:* `mqc_cuest_context.f90`, `mqc_cuest_integrals.f90`

Add to `cuest_context_t`, alongside the existing pools:

```fortran
type(device_pool_t) :: scratch_j    !! Coulomb output
type(device_pool_t) :: scratch_k    !! Exchange output
type(device_pool_t) :: scratch_xc   !! XC potential output
type(device_pool_t) :: scratch_fock !! Assembled Fock, stays resident
type(device_pool_t) :: scratch_core !! Core Hamiltonian, uploaded once
type(device_pool_t) :: scratch_ovlp !! Overlap, uploaded once
```

Wire every one into `context_destroy`'s `release()` list in the same change —
a pool that is allocated but never released leaks per fragment, and fragmented
runs do thousands of fragments per rank.

Give each compute routine a variant that writes to a caller-chosen output
buffer and does **not** fetch:

```fortran
call system%compute_coulomb_device(d_density, d_out, error)
```

Keep the existing host-fetching routines for now: the gradient path
(`mqc_cuest_gradient.f90`) uses them and is out of scope.

*Cost:* 3 × `n_ao²` device doubles ≈ 24 MB at n_ao = 1000, against 40–80 GB of
card. Cheap.

*Verify:* energies unchanged on `run_cuest_validation.sh`. No behaviour change
yet — same transfers, different buffers.

### Step 2 — Assemble the Fock on device

*Files:* `mqc_cuest_integrals.f90`, `mqc_cuest_scf.f90`

`core_hamiltonian` and `overlap` are **iteration-invariant**. Upload once,
before the SCF loop, never again. This is free and is skipped surprisingly
often.

Then `fock = h + J - K + Vxc` becomes four `cublasDaxpy` calls into
`scratch_fock` (or one small kernel). The density still comes H2D each
iteration — it is produced on the host by `build_density` from host
diagonalisation, and stays that way until step 5.

Per-iteration traffic after this step:

| | before | after |
|---|---|---|
| H2D | density ×1, occupied ×2 | density ×1, occupied ×1 |
| D2H | J, K, Vxc | assembled Fock ×1 |
| syncs | 3 | 1 |

*Verify:* energies unchanged. This is the step that actually pays.

### Step 3 — Device-resident DIIS history

*Files:* new `mqc_diis_device.f90`, `mqc_cuest_scf.f90`

Mirror `diis_state_t`'s interface (`init`/`push`/`extrapolate`/`destroy`) over
device buffers. Keep the host `mqc_diis` — it is the reference the device
version is diffed against, which is the only way to catch the characteristic
failure (a stale history producing plausible-but-wrong coefficients, no crash).

- history: one `device_pool_t` of `n_fock × max_vectors`, ditto errors
- `push`: `cublasDcopy` into the ring slot — device→device, no transfer
- overlap update: **one `cublasDgemv`** of the new error vector against the
  whole history, not `n_stored` × `cublasDdot`. `Ddot` returns its result to
  the host, so the naive loop costs one synchronise per stored vector.
- the (`n_stored`+1)² solve stays on the **host**. It is tiny, and `solve_diis`
  already exists. Pulling back an 9×9 matrix is nothing; a cuSOLVER dependency
  for it would be absurd.
- `extrapolate`: `cublasDscal` + `n_stored` × `cublasDaxpy` into `scratch_fock`

*Verify:* run the same molecule with host and device DIIS and require identical
SCF iteration counts and energies to ~1e-10. Do this **before** step 4 — it is
much easier to localise a DIIS bug while the commutator is still on the host.

### Step 4 — Commutator and basis transform on device

*Files:* `mqc_cuest_scf.f90`

The last host algebra in the loop:

```
FDS - SDF        4 × cublasDgemm
project to ortho 2 × cublasDgemm
error norm       1 × cublasDnrm2
```

cuBLAS is column-major, so leading dimensions are the Fortran ones and no
transposition is needed — unlike cuEST's row-major matrices (see the header of
`mqc_cuest_integrals.f90`). Scalars are host references under the default
`CUBLAS_POINTER_MODE_HOST`.

After this the only per-iteration host traffic is the density H2D, the Fock
D2H for diagonalisation, and two scalars (energy, error norm).

### Step 5 — (optional, later) diagonalisation

`diagonalize_fock` is the remaining host step and the reason the Fock still
comes back. Moving it needs cuSOLVER (`cusolverDnDsygvd`), which is a new
binding set. Do this only if a profile says the transfer and the host
diagonalisation actually matter — at moderate `n_ao` they may not, and it is
by far the largest piece of new binding work.

---

## Rules that are not optional

Under separate device memory, the failure mode is **silent**: a kernel reading
an unmapped array gets stale host memory, no crash, plausible numbers. Hence:

1. **A host build passing proves nothing about device data motion.** Verify on
   the A100. This box has no libcuest and sm_70 cards, so it can only
   compile-check — and cuEST needs sm_80.
2. **Every new `device_pool_t` gets wired into `context_destroy` in the same
   commit that adds it.** A missed one leaks per fragment.
3. **Never `cudaMemcpy` a whole derived type** that has allocatable components —
   move the component buffers.
4. **`cublasDdot` and `cublasDnrm2` synchronise** (host-pointer results). Do not
   put them in a loop; batch with `Dgemv` or accept one sync per iteration.
5. **No host fallback path** (decided): the validation suite is the oracle. If
   residency breaks, energies move or the SCF stops converging.

## How to verify at each step

`validation/run_cuest_validation.sh` on the A100, which has reference energies.
Between steps, the invariant is: *energies unchanged, iteration counts
unchanged*. If iteration count moves but the energy converges to the same
place, suspect DIIS; if the energy moves, suspect the Fock assembly.

For step 3 specifically, run host and device DIIS on the same molecule in
separate runs and diff the per-iteration energy trace — a residency bug shows
up as divergence after the first few iterations, not at iteration 1.

## Build

GPU builds need tblite off — toml-f contains a backslash string literal
nvfortran rejects with no flag to disable:

```
-DMQC_ENABLE_TBLITE=OFF -DMQC_ENABLE_CUEST=ON -DCUEST_ROOT=<archive>
-DPIC_USE_LEGACY_MPI=ON -DPIC_MPICH_PERLMUTTER=ON
```

The last flag compiles out `win_allocate` wrappers whose generic does not
resolve under nvfortran with Cray MPICH. metalquicha uses no RMA, so nothing is
lost. If a second such mismatch appears, switch to `PIC_USE_VAPAA=ON` rather
than patching instances — it goes through the MPI C ABI and fixes the class.

Note the cuBLAS bindings are `iso_c_binding`, so nothing here forces nvfortran.
That was deliberate: `cudafor` would have made the device path
nvfortran-only, and since tblite cannot build there, a single build would have
had to choose between the xTB backend and the GPU backend.

==============================================
Ligand interaction spectra, residue by residue
==============================================

``tools/sapt_scan/sapt_ligand_scan.py`` runs one SAPT calculation per
residue--ligand pair and collects them into one table: residue, distance, and
the named SAPT terms, strongest interaction first.

Why this route exists
=====================

Three fragmented methods could in principle answer "what is this ligand's
interaction with each residue, broken into physical terms", and none of them
does today.

* **EFMO** has the named terms per pair -- electrostatics, exchange repulsion,
  dispersion, charge transfer -- and **cannot cut a protein at all**.
  ``efmo_options_t`` carries no bond-breaking option.
* **FMO** cuts a protein and reports a pair energy that is a single number.
  Decomposing it is a structural change, not a keyword.
* **SAPT** decomposes the interaction of *exactly two* monomers, and reports
  twelve named terms at SAPT0 or eighteen at SAPT2.

A ligand is not covalently bonded to anything. So every residue--ligand pair
is a legitimate two-monomer problem, and the spectrum is one SAPT run per
residue. That is all this harness is: no new theory, and nothing the program
could not already do one pair at a time.

Running it
==========

The input is an ordinary fragmented deck -- the same one a many-body run would
take, with ``fragments`` and ``connectivity`` already declared -- plus which
fragment is the ligand:

.. code-block:: bash

   python3 tools/sapt_scan/sapt_ligand_scan.py complex.json \
       --ligand 4 --cutoff 6.0 --basis sto-3g \
       --workers 8 -o ligand_spectrum.csv

Residues are selected by closest approach to the ligand, the metric the
fragmented path screens on. The runs are independent, so ``--workers`` spreads
them across processes; each is one small serial job.

A worked example ships in ``tools/sapt_scan/example/``: a glycine tripeptide
with a water hydrogen bonded to the residue-2 carbonyl at 2.80 Å O···O.

.. code-block:: text

   residue,distance,n_caps,elst10,exch10,...,total
   1,4.945,1, 0.026850, 0.000969,...,  0.011800
   3,3.609,1, 0.214503, 0.035149,...,  0.190911
   2,1.840,2,-7.994680,10.815556,...,  0.270803

Residue 2 is the hydrogen-bonded one, and its row is the point of the exercise:
an electrostatic attraction of −8.0 kcal/mol almost cancelled by +10.8 of
exchange repulsion. A single pair energy cannot tell you that, and the near
cancellation is exactly what a medicinal chemist wants to see.

Two things these numbers do not include
=======================================

Both are real physics left out, not caveats about numerics. Neither is
recoverable by running the scan differently.

**The monomers are isolated.** SAPT sees a residue and a ligand in vacuum.
They do not feel the rest of the protein. FMO's monomers would, through its
embedding field; SAPT's have no such mechanism, because the theory partitions
the Hamiltonian into exactly ``H_A + H_B + V`` and there is no slot for an
environment. For a solvent-exposed site this is a modest omission. **For a
buried site it is not** -- the electrostatic term especially is the one most
changed by surrounding charge, and it is the largest term in a polar contact.
Treat the spectrum as the interaction those two fragments would have alone.

**Each residue is capped where it leaves the backbone.** A residue cut out of
a chain has two dangling valences; the harness caps them with hydrogens, by
the same rule the many-body path uses:
``R_H = R_kept + s (R_gone - R_kept)``, with ``s`` from ``--cap-scale``.

The default ``s = 1.0`` matches the rest of the program, and it places the cap
hydrogen **exactly on the position of the heavy atom it replaces** -- a ~1.37 Å
N--H. That is deliberate in the code rather than a defect, but it is not a
chemical bond length, and cap placement moves a hydrogen-bonded pair energy by
about **eight per cent**. Two scans run at different ``--cap-scale`` are not
comparable, and nothing in the numbers says which was used, so the harness
writes it into the table's header. Quote it when you quote the terms.

The capping is done by the harness, not the program
===================================================

Worth knowing before adapting this. The SAPT route in ``mqc_driver.f90``
slices its two monomers straight out of the geometry: it never calls
``build_fragment_from_indices``, adds no hydrogen caps, and **reads no**
``connectivity``. Declaring bonds in a SAPT deck does nothing at all.

So the caps are placed by the harness and written into the geometry as real
atoms before the deck is generated. A residue handed to SAPT uncapped is a
radical with a dangling valence, which either dies on an electron-count parity
error or converges to a wrong number.

What the total is worth
=======================

If the terms are to be quoted, the total they sum to needs a measured error.

``tools/sapt_scan/validate_against_counterpoise.py`` measures it, against the
one SAPT term with an independent definition: ``e_int_hf_cp``, the
counterpoise-corrected Hartree--Fock interaction energy. A VMFC counterpoise
run over the same two fragments computes that same quantity by a different
route. On the worked example's hydrogen-bonded pair:

.. code-block:: text

   counterpoise-corrected supermolecular HF : 0.001564182246 Ha  (+0.9815 kcal/mol)
   SAPT e_int_hf_cp                         : 0.001564182247 Ha  (+0.9815 kcal/mol)
   difference                               : 6.68e-13 Ha

That is SCF convergence, which is to say the two routes compute the same
number and the agreement is an identity rather than an approximation.

``total`` is deliberately not what is compared. It carries dispersion, which
Hartree--Fock does not have, so it cannot equal a supermolecular HF
interaction and a disagreement there would mean nothing. On this pair the HF
interaction is +0.98 kcal/mol and the SAPT0 total is +0.27: the difference is
the dispersion SAPT adds and HF omits.

Detecting a failed pair
=======================

A refusal in this program **exits with status zero**. The harness therefore
never reads an exit code: it looks for the ``sapt`` section in the output,
checks every term is finite, and reports the refusal line from the log when
one is missing. A pair that fails becomes a row carrying its error rather than
a silently absent row, and the harness itself exits non-zero so a caller that
does check is told.

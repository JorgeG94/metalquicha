Symmetry-Adapted Perturbation Theory (SAPT0)
============================================

SAPT computes the interaction energy between two molecules **and what it is made
of**. That second part is the point. A supermolecular calculation subtracts two
numbers and gives you a third; it cannot tell you whether a dimer is bound by
electrostatics, by dispersion, or by an induction term that a smaller basis
would have missed. SAPT builds the interaction as a perturbation series in the
intermolecular operator from the start, so every contribution arrives named.

``SAPT0`` is the first rung: monomers at Hartree-Fock, the interaction through
second order. It is the level that scales well enough to run on things you care
about, and the level whose terms map onto the physics an EFP potential also
carries.

Running one
-----------

.. code-block:: json

   {
     "schema": {"name": "water_dimer_sapt", "version": "1.0"},
     "molecules": [{
       "xyz": "w2_dimer.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1,
       "fragments": [[0, 1, 2], [3, 4, 5]]
     }],
     "model": {"method": "sapt0", "basis": "6-31g"},
     "driver": "Energy"
   }

.. code-block:: bash

   ./mqc water_dimer_sapt.json

There is no ``sapt`` keyword block. The monomers are the deck's own
``fragments``, because a SAPT deck is exactly a geometry, a basis, and the
partition that says which atoms are which monomer. Atom indices are 0-based, as
everywhere else.

``"sapt"`` on its own is accepted as a spelling of ``"sapt0"``. SAPT0 is the only
order implemented, so the bare name is not ambiguous.

Exactly two fragments
---------------------

SAPT partitions the Hamiltonian as :math:`H_A + H_B + V`. There is no slot for a
third monomer, so a deck with any number of fragments other than two is refused
rather than approximated. This is the theory and not a limitation of this
program: a cluster is one SAPT calculation *per pair*, and the sum over pairs is
the two-body part of its binding energy, not its binding energy. Three-body
induction in a hydrogen-bonded cluster is worth several kcal/mol and no sum over
pairs contains it.

``validation/check_sapt.f90`` walks the fifteen pairs of a six-water prism if you
want to see what that looks like.

Only rank zero does the work. The pairs of a cluster are the obvious thing to
distribute; a single pair is not.

The basis
---------

Every quantity is computed in the **dimer-centred basis**: each monomer carries
the other's basis functions on ghost atoms, with basis functions and no nuclear
charge. This is not an option, it is what makes the terms comparable -- a
monomer described in its own basis and then again in the dimer's would differ by
basis-set superposition error, and that error would land in the interaction
energy indistinguishably from physics. The supermolecular Hartree-Fock reference
reported alongside is counterpoise-corrected for the same reason.

The practical consequence is cost: SAPT0 on a dimer of :math:`n` basis functions
each is not two calculations of size :math:`n` but two of size :math:`2n`.

What comes out
--------------

.. code-block:: text

     SAPT0 interaction energy, Hartree
       electrostatics        .006384228641
       exchange              .003241500424
         (S^2 approximation) .003235393995
       induction             -.001134897266
         (uncoupled)         -.000954737585
       exchange-induction    .000949574156
         (uncoupled)         .000805018642
       dispersion            -.000317465445
       exchange-dispersion   .000117840543
       delta HF              -.000124862284
       --
       supermolecular HF     .009315543671
       total                 .009115918769

**Electrostatics** (:math:`E^{(10)}_{\text{elst}}`) is the classical Coulomb
interaction of the two unperturbed monomer densities. At long range it becomes
the multipole electrostatics an EFP potential encodes; at short range it differs
from it, and that difference is charge penetration.

**Exchange** (:math:`E^{(10)}_{\text{exch}}`) is the repulsion from
antisymmetrizing the product wavefunction -- the Pauli principle, priced. It is
reported both in the single-exchange :math:`S^2` approximation and to all orders
in overlap; codes quote one or the other, and a term can only be checked against
the same term.

**Induction** (:math:`E^{(20)}_{\text{ind,resp}}`) is each monomer polarizing in
the other's field, and **exchange-induction** is the repulsive correction from
antisymmetrizing that. They come in *response* (coupled, the ones entering the
total) and *uncoupled* forms; both are printed for the same reason as above.

**Dispersion** (:math:`E^{(20)}_{\text{disp}}`) is the correlation of the two
monomers' fluctuations -- the term Hartree-Fock has none of, which is why a
supermolecular HF interaction energy is not an interaction energy. At long range
it goes as :math:`-C_6/R^6`. **Exchange-dispersion** is its antisymmetry
correction.

**delta HF** (:math:`\delta^{(2)}_{\text{HF}}`) collects induction beyond second
order, by difference: it is whatever the counterpoise-corrected supermolecular
Hartree-Fock interaction energy has that the first- and second-order induction
and electrostatic terms do not. So the identity

.. math::

   E^{(10)}_{\text{elst}} + E^{(10)}_{\text{exch}}
   + E^{(20)}_{\text{ind,resp}} + E^{(20)}_{\text{exch-ind,resp}}
   + \delta^{(2)}_{\text{HF}} = E^{\text{HF}}_{\text{int}}\text{(CP)}

holds exactly, by construction rather than by agreement. It is worth knowing
that it is a definition and not a result: it will hold even if several of the
terms in it are wrong.

The total is the sum of everything above the rule.

In the JSON output
------------------

Under ``sapt``, all twelve numbers, keyed by their names in the literature
rather than the console's prose -- their reason to exist is being compared
against another program's output:

.. code-block:: json

   "sapt": {
     "elst10":        0.006384228641113765,
     "exch10_s2":     0.0032353939945464,
     "exch10":        0.0032415004240945533,
     "ind20_u":      -0.0009547375845903505,
     "ind20_r":      -0.001134897265892534,
     "exch_ind20_u":  0.0008050186417872887,
     "exch_ind20_r":  0.0009495741556526977,
     "disp20":       -0.00031746544457269746,
     "exch_disp20":   0.0001178405429812502,
     "delta_hf":     -0.00012486228430180941,
     "e_int_hf_cp":   0.00931554367066667,
     "total":         0.009115918769075223
   }

``_u`` is uncoupled and ``_r`` is response; ``_s2`` is the single-exchange
approximation against the full :math:`S^\infty`. The total also appears as
``total_energy`` like any other method's, but on its own it is the one number a
supermolecular calculation would have given you too.

How it is checked
-----------------

Two ways, because they fail differently.

**Against a reference implementation.** ``validation/check_sapt0.py`` is a
conventional four-index SAPT0 in PySCF, term by term. The decks in
``validation/inputs/cpu/mqc/sapt/`` carry its totals. Note the provenance
honestly: psi4 is the obvious reference and cannot serve as one here, because
its closed-shell SAPT is density-fitted by construction and cannot produce
conventional numbers at all. Every apparent disagreement with psi4 was chased
down to density fitting by reproducing psi4's own settings to sub-nanohartree,
which is a weaker guarantee than the rest of the validation suite has.

**Against the long-range limits.** Each term's asymptotic form is fixed by the
physics rather than by us, and this is the check that catches a term which is
right at one geometry by coincidence:

.. code-block:: bash

   ./build/check_sapt --scan

walks a water dimer from 4 to 12 Angstrom and measures each exponent locally
between consecutive separations, so the answer is a sequence converging on the
limit rather than one fit a single bad point could drag into place:

=============  ==================================================  ========
Term           Measured exponent in :math:`|E| \sim R^{-n}`         Expected
=============  ==================================================  ========
``elst10``     3.19, 3.16, 3.12, 3.09, 3.07, 3.05                  3
``ind20_r``    7.40, 6.24, 6.17, 6.13, 6.10, 6.07                  6
``disp20``     6.36, 5.97, 5.96, 5.97, 5.98, 5.98                  6
``exch10``     28, 44, 59 -- no power law fits it                  exponential
=============  ==================================================  ========

Each approaches its asymptote from above, as the higher multipole ranks die
away. Exchange is an overlap effect, so it decays exponentially and a power-law
fit reports that by returning a nonsense exponent that keeps growing.

Do not read too much into a single separation below about 8 Angstrom, and note
that exchange runs out of double precision past roughly 7: it is below
:math:`10^{-16}` there, assembled from traces of matrices of order
:math:`10^{2}`, and the scan marks those rows rather than reporting them as
though they meant something.

What is not here
----------------

Only SAPT0. The higher levels -- SAPT2, SAPT2+, SAPT2+(3), SAPT2+3 -- add
intramonomer electron correlation, which SAPT0 has none of. There is no SAPT1:
the numbering jumps because first-order intramonomer correlation alone is
unbalanced, correcting electrostatics before exchange, and would be worse than
SAPT0 rather than better.

Closed-shell only, and no density fitting, so cost grows quickly with basis. The
counterpoise-corrected supermolecular Hartree-Fock reference is computed as part
of every run, which means a SAPT0 calculation includes a dimer SCF in the full
basis whether or not you wanted one.

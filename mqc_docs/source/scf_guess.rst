==================
The Initial Guess
==================

An SCF starts from a density it did not compute, and where it starts decides how
many iterations it takes, and occasionally *which answer it finds*. The second
half matters more than it sounds: an SCF converges to a stationary point, not to
the ground state, and a bad start can converge perfectly well onto the wrong one.

This page covers the guesses Metalquicha has, what each is for, and the one that
is not a formula but a sequence of calculations -- basis set projection.

Choosing one
------------

There are two spellings. ``keywords.guess.type`` is the current one and carries
settings alongside the name; ``keywords.scf.guess`` is the older spelling, still
read, and takes the name alone. **Setting both is refused** -- a deck that says
both has not said which it means.

``keywords.scf.guess``, default ``auto``:

.. code-block:: json

   {
     "keywords": {
       "scf": { "guess": "sad" }
     }
   }

============================ =========================================================
``core``                     The core Hamiltonian, :math:`F = H`. No electron
                             repulsion at all.
``gwh``                      Generalized Wolfsberg-Helmholz. Needs no atomic
                             calculation.
``sad``                      Superposition of spherically averaged atomic densities.
``sac``                      Superposition of the free atoms' own spin densities.
``basis_set_projection``     Converge in smaller bases first and project forward.
``auto``                     Let the backend choose.
============================ =========================================================

``auto`` resolves per backend -- ``sad`` on the CPU path, ``gwh`` on cuEST --
because the best starting point is a property of the backend rather than of the
request, and each measured its own. An explicit spelling always wins over both.

What each one is for
--------------------

**core** is the cheapest and the worst. With no electron repulsion in the
starting Fock matrix the orbitals it produces are badly wrong, and the SCF has to
recover from that. On a 46-atom peptide cation it does not converge at all in 100
cycles. It is useful for debugging -- it is the one guess that depends on nothing
-- and almost never otherwise.

**gwh** estimates the Fock matrix from orbital energies and overlaps without
solving anything, so it costs nothing beyond the integrals already in hand. It is
a large improvement on ``core`` for no work, which is why it is the fallback when
something better is unavailable.

**sad** solves each *element* once as a free atom, caches it, and superposes the
spherically averaged densities. The atomic solves are reused across every atom of
that element and across fragments, so the cost is per element rather than per
atom. It is the default on the CPU path and the right choice for most closed-shell
molecules.

**sac** is ``sad``'s open-shell counterpart: it superposes the atoms' own spin
densities rather than spherically averaging them. It therefore **arrives with its
spatial symmetry already broken**, which is what a radical or a stretched bond
wants and what a closed shell does not. Reach for it on an unrestricted
calculation that keeps collapsing onto the restricted solution.

.. _basis-set-projection:

Basis set projection
--------------------

The others are formulas evaluated once. This one is a *ladder of SCF
calculations*: converge the molecule in a small basis, project that density into
a larger one, converge again, and repeat until the target basis is reached.

It is the answer to a specific problem -- an SCF in a large basis that will not
converge, or converges to the wrong state -- and it works because the difficulty
is usually not in the basis. A molecule that is hard to converge in cc-pVQZ is
generally hard for chemical reasons that are already visible in STO-3G, where an
iteration costs a tiny fraction as much. Solve it where it is cheap, and hand the
answer to where it is expensive.

Configuring it
^^^^^^^^^^^^^^

The ladder needs the rungs written down, so this guess uses its own block rather
than a bare name:

.. code-block:: json

   {
     "model": {
       "method": "hf",
       "basis": "cc-pvtz"
     },
     "keywords": {
       "guess": {
         "type": "basis_set_projection",
         "subscf": {
           "steps": [
             { "basis": "sto-3g", "maxiter": 30, "tolerance": 1e-4 },
             { "basis": "6-31g",  "maxiter": 50, "tolerance": 1e-6 }
           ]
         }
       }
     }
   }

That deck runs **three** SCFs: STO-3G, then 6-31G started from the STO-3G
density, then cc-pVTZ started from the 6-31G density.

Only the last of those prints a full iteration table. The rungs report one line
each, at ``Verbose``, and those lines are how you confirm the ladder ran at all:

.. code-block:: text

   guess rung 1: sto-3g, 7 functions, 6 iterations
   guess rung 2: 6-31g, 13 functions, 7 iterations
   initial guess: basis_set_projection

**The target basis is not a rung.** It comes from ``model.basis`` and is never
listed in ``steps``; the last rung is the last *preliminary* calculation. Two
steps plus a model basis is three SCFs. Listing the target basis as a final step
would converge it twice.

Order is the order written, cheapest first. Each rung starts from the previous
rung's density rather than from scratch, so the ladder gets cheaper as it climbs
even though the bases grow -- the first SCF is the only one starting from nothing.

Per-rung settings
^^^^^^^^^^^^^^^^^

======================= ============== ======================================
``basis``               *required*     Basis for this rung.
``maxiter``             50             Iteration cap for this rung.
``tolerance``           ``1e-5``       Convergence tolerance for this rung.
======================= ============== ======================================

The convergence settings are per rung because the rungs are not doing the same
job. An early rung only has to land in the right basin and can stop loosely; the
last one before the target may as well be tight, since it is cheap relative to
what follows. Neither is bound by ``keywords.scf.tolerance``, which belongs to
the target SCF.

Rules the deck has to satisfy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three, all refused at validation rather than at run time:

- **Do not set both** ``keywords.scf.guess`` **and** ``keywords.guess.type``:
  *"keywords.scf.guess and keywords.guess.type both set the initial guess. Use
  keywords.guess.type; keywords.scf.guess is the older spelling and is
  superseded."*
- ``keywords.guess.subscf`` **only means something under**
  ``type: basis_set_projection``: *"beside any other guess it would be read and
  then ignored"*, so it is refused rather than accepted and discarded.
- **A projection guess needs at least one step**, *"with at least one basis to
  converge before the target one"*. A ladder with no rungs is the default guess
  wearing a longer name.

When the ladder cannot be built
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If a rung fails -- a basis name that does not resolve, or one that does not cover
an element in the molecule -- the run **does not fail**. It warns and falls back
to ``sad``, on the same reasoning the atomic guess uses: a guess that cannot be
built is a reason to start somewhere else, not to abandon the calculation.

.. code-block:: text

   basis projection guess: rung 1 (no-such-basis): no basis set file for
   'no-such-basis': ... -- falling back to sad

The warning is loud, and worth reading rather than skipping, because the run is
now doing something other than what the deck asked for. If a projection guess
was being used to reach a particular state, the fallback will not reach it, and
the number that comes out will be for a different stationary point.

When it is worth the cost
^^^^^^^^^^^^^^^^^^^^^^^^^

It is not free: the preliminary SCFs are real calculations. It pays when

- a large-basis SCF oscillates or stalls and a smaller one converges cleanly;
- diffuse functions are the problem. A set with diffuse functions on a molecule
  that does not need them is nearly linearly dependent, and the same molecule in
  the parent set without them is not. Converge there first;
- the SCF converges but to the wrong state, and a smaller basis finds the right
  one. A projected density carries the *character* of the state forward, which
  ``sad`` cannot: a superposition of atoms knows nothing about which molecular
  state was wanted.

It is the wrong tool for an SCF that merely crawls because the gap is small --
that is what :ref:`level-shifting` is for -- and for one that oscillates in a
basis that is not itself the difficulty, where an accelerator is cheaper. See
:ref:`accelerators`.

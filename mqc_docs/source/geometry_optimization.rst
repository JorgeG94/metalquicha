Geometry Optimization
=====================

The ``Optimize`` driver minimizes a structure, calling the method in the deck for
an energy and gradient at each step until the largest gradient component falls
below the convergence threshold.

The optimizer itself is `DL-FIND <https://www.chemshell.org/dl-find>`_, reached
through `libdlfind <https://github.com/digital-chemistry-laboratory/libdlfind>`_.
It is **not built by default** -- see :ref:`optimizer-build` below for the flag
and the reason.

Running one
-----------

.. code-block:: json

   {
     "schema": {"name": "water_opt", "version": "1.0"},
     "molecules": [{
       "xyz": "water.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "gfn2"},
     "driver": "Optimize",
     "keywords": {
       "optimization": {"max_steps": 100, "algorithm": "lbfgs"}
     }
   }

.. code-block:: bash

   ./mqc water_opt.json

``Optimization`` and ``Opt`` are accepted spellings of the driver.

An optimization is a loop *over* calculations rather than one of them, so it works
for a single molecule from an input deck. A multi-molecule deck, a Python session
and the C API are refused by name rather than left to misbehave.

Which methods can be optimized
------------------------------

Any method that produces a gradient, and no others. The optimizer never branches
on the method -- it asks for a gradient and receives a whole-system one, with
fragmentation and hydrogen caps already folded back in -- so this table is about
which backends have gradients, not about the optimizer:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Backend
     - Gradients
     - Optimization
   * - tblite (``gfn1``, ``gfn2``)
     - yes
     - Supported and tested, fragmented and not
   * - libcint (CPU ``hf``, ``dft``, ``mp2``, ``ri-mp2``)
     - yes
     - Supported, for whatever combination has a gradient -- see
       :doc:`capabilities`, which lists both what is covered and what is
       refused
   * - libcint (CPU ``ccsd``, ``ccsd(t)``)
     - **no**
     - Refused: a coupled cluster gradient needs the Lambda amplitudes
   * - cuEST (GPU ``hf``, ``dft``)
     - yes
     - Expected to work; needs an sm_80 card

There is no optimizer-level restriction on the backend. What refuses a gradient
refuses it wherever that gradient would have been built, and the optimizer finds
out the same way any other caller would -- by taking one gradient up front,
before DL-FIND is entered:

.. code-block:: text

   The initial gradient check failed: a coupled cluster gradient needs the
   Lambda amplitudes, which are not implemented. Run the gradient at the
   Hartree-Fock or MP2 level, or coupled cluster as an energy.
   Geometry optimization failed: This method cannot produce a gradient, so it
   cannot be optimized.

That probe costs one evaluation on a run that does work. It is there because the
refusal is raised where the gradient would have been built, with no way to ask in
advance, and reaching DL-FIND with a failed first step ends the process with a
Fortran backtrace instead of a message.

Fragmented optimization
-----------------------

Fragmented systems optimize the same way, with the MBE or GMBE gradient standing in
for the whole-system one. Nothing extra is needed in the deck beyond the usual
``keywords.fragmentation`` block.

.. _frozen-terms:

The term list is frozen
^^^^^^^^^^^^^^^^^^^^^^^

This is the part worth understanding before running a large system.

A fragmented run chooses its term list from the geometry it is given, screening out
the n-mers whose monomers are further apart than the cutoff. Regenerate that at
every step of an optimization and **the list changes as the molecules move**: a
dimer crossing the cutoff enters the sum, and the total energy jumps by that
dimer's whole interaction. The optimizer has no way to know the function changed
shape rather than the geometry being bad. It reads the jump as real, rejects the
step, shrinks the trust radius, and either stalls or settles somewhere that is not
a stationary point of anything. Nothing fails; every number looks plausible.

So the list is generated once, at the starting geometry, and the same one is used
for every step. The energy is then a fixed sum over a fixed set of subsystems -- a
smooth function of the coordinates -- and its gradient is that function's actual
derivative, which is what makes the optimization well posed.

Measured on a water tetramer with a ``4.0`` Angstrom dimer cutoff, where three
dimers cross inward during the run:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - ``freeze_terms``
     - Term list across the run
     - Outcome
   * - ``true`` (default)
     - 7 terms throughout
     - Converged
   * - ``false``
     - 7 → 8 → 9 terms
     - Did not converge in 200 steps

The energy discontinuities where terms entered were ``1.2e-05`` and ``1.9e-04``
Hartree. The second is about 185 times the energy convergence criterion, so the
re-screened run is chasing a target that moves further than its own tolerance.

What this gives up is worth stating: a term screened out at the starting geometry
stays out even if those fragments later approach. That is a definite approximation,
fixed for the run and reported in the log, rather than a surface that changes
underneath the optimizer. If the fragments move a long way, loosen the cutoffs.

**GMBE cannot be frozen this way.** Its inclusion--exclusion terms are rederived
from overlapping primaries rather than taken from a supplied list, so a GMBE
optimization with cutoffs warns instead, and the caution above applies to it.

Coordinate systems
------------------

``cartesian`` is the default and always works. ``hdlc`` -- hybrid delocalised
internal coordinates -- builds internal coordinates within each *residue* and keeps
Cartesians between them, which is the right shape for a molecular cluster. The
hybrid is of the two coordinate types, not of two levels of hierarchy.

**The residue partition is HDLC's definition, not a tuning knob.** A fragmented
calculation already has the grouping and its monomers are used directly. An
unfragmented one perceives its own connectivity and makes each covalently connected
molecule a residue, so a water dimer optimizes as two residues without the deck
having to request a fragmented energy it did not want. A single molecule comes back
as one residue, because that is what its connectivity says.

That perception matters more than it looks. Giving HDLC one residue for a system
that comes apart makes it build bonds, angles and torsions across a gap where there
is no bond; the primitive set is deficient, the delocalised basis built from it is
badly conditioned, and the optimizer does not fail -- it *converges early*. A water
dimer that reports success at 20 steps and -10.148997 Hartree with one residue goes
to 23 steps and -10.149007 with two, a hundredth of a millihartree lower. Nothing
in the output said the first answer was worse, which is why the residue count is
now reported at ``Verbose``::

    optimizer residues: 2 for 6 atoms

``dlc`` puts every atom in one residue, which suits a single connected molecule and
**fails on a cluster** -- there is no connected internal-coordinate system spanning
molecules that are not bonded to each other, and DL-FIND says so
(``cyclic failure at residue 1``).

``hdlc-tc`` and ``dlc-tc`` are the *total connection* variants. They build the
primitive internals of a residue from every atom pair in it rather than from the
perceived bonds, which is the answer when the connectivity is not what a
covalent-radius rule says it is: a weakly bound complex the perception splits, an
ion pair it merges, a transition state where a bond is half formed. ``dlc-tc``
therefore does not fail on a cluster the way plain ``dlc`` does -- the total
connection supplies the coordinates no bond provides.

The cost is quadratic in residue size where a bond list is linear, so ``hdlc-tc``
is for systems whose residues are small and ``dlc-tc``, which totally connects
everything at once, for systems that are small altogether. Prefer ``hdlc``: the
perceived partition is right far more often than not, and this is the escape hatch
for when it is not.

For a well-behaved starting structure the difference is modest -- a water trimer
takes 20 evaluations in Cartesians against 16 in HDLC. Far from the minimum, the
count is dominated by how far the structure has to travel, not by the coordinate
system: a hand-built trimer whose three molecules all have to rotate and translate
into a hydrogen-bonded ring takes around a hundred either way.

What a run leaves behind
------------------------

Four files, named after the deck:

``output_<name>.json``
   The ordinary output document, from a final single point at the optimized
   geometry. Every intermediate step is run without output, so this is the energy,
   dipole and gap of the structure you ended up with.

``output_<name>_optimized.xyz``
   The optimized structure, in Angstrom. The comment line says ``converged`` or
   ``NOT CONVERGED``, because an ``.xyz`` outlives the log it came from.

``output_<name>_optimization.json``
   The machine-readable record: ``converged`` as a boolean, step and evaluation
   counts, final energy, largest gradient component, the threshold it was compared
   against, the final geometry, and the whole trajectory. This is what to parse.

``output_<name>_trajectory.xyz``
   Multi-frame, one frame per accepted step, energy on each comment line. Opens in
   any viewer. Turn it off with ``"trajectory": false`` for a large system over many
   steps -- the path is held in memory until the run ends.

Convergence is reported honestly. DL-FIND hands back geometries but no verdict, so
the largest gradient component of the last step is tested here; a run that exhausted
``max_steps`` says so, exits non-zero, and still writes its geometry -- which is what
you want to restart from.

Keywords
--------

All optional. Everything under ``keywords.optimization``:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Key
     - Default
     - Meaning
   * - ``max_steps``
     - ``100``
     - Give up after this many steps. ``steps`` is accepted too
   * - ``gradient_tolerance``
     - ``4.5e-4``
     - Convergence on the largest gradient component, Hartree/Bohr.
       ``tolerance`` is accepted too
   * - ``energy_tolerance``
     - engine
     - Convergence on the energy change, Hartree
   * - ``max_step``
     - engine
     - Longest step the optimizer may take, Bohr
   * - ``coordinates``
     - ``cartesian``
     - ``cartesian``, ``hdlc``, ``hdlc-tc``, ``dlc`` or ``dlc-tc``.
       ``coordinate_system`` is accepted too
   * - ``algorithm``
     - ``lbfgs``
     - ``lbfgs``, ``cg``, ``cg-auto``, ``sd``, ``prfo``, ``nr`` or ``damped``.
       ``optimizer`` is accepted too. See :ref:`opt-algorithms`
   * - ``lbfgs_memory``
     - engine
     - Steps of curvature history L-BFGS keeps
   * - ``freeze_terms``
     - ``true``
     - Fix the MBE term list at the starting geometry. See :ref:`frozen-terms`
   * - ``trajectory``
     - ``true``
     - Record every accepted geometry
   * - ``hess_end``
     - ``false``
     - Hessian at the converged geometry, to say whether it is a minimum. See
       :ref:`hess-end`
   * - ``hessian_update``
     - engine
     - ``none``, ``powell``, ``bofill`` or ``auto``. Only the algorithms that
       hold a Hessian read it
   * - ``target``
     - ``minimum``
     - ``minimum`` or ``saddle``. What the run is looking for, which decides
       how ``hess_end`` reads its own result. See :ref:`saddle-search`
   * - ``endpoint``
     - none
     - Path to a second geometry. Turns the run into a chain of states between
       it and the starting structure. See :ref:`neb`
   * - ``images``
     - engine
     - Images along the path, endpoints included. At least 3
   * - ``neb_spring``
     - engine
     - Spring constant holding neighbouring images apart
   * - ``neb_endpoints``
     - ``frozen``
     - ``frozen``, ``perpendicular`` or ``free``
   * - ``frozen_atoms``
     - none
     - 0-based indices held at their input position. See :ref:`opt-constraints`
   * - ``constraints``
     - none
     - Internal coordinates held fixed. See :ref:`opt-constraints`
   * - ``timestep``, ``friction``, ``friction_factor``, ``friction_rising``
     - engine
     - ``damped`` only
   * - ``print_level``
     - follows log
     - How much DL-FIND itself prints

Where a default is given as "engine", the setting is left to DL-FIND rather than
overridden with a second number chosen here.

.. _opt-algorithms:

Choosing an algorithm
---------------------

Seven, of which six go downhill and one climbs. ``nr`` is the awkward one: it
goes to whichever stationary point it starts nearest, which is a minimum only if
it began in that basin.

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - ``algorithm``
     - What it does
   * - ``lbfgs``
     - Limited-memory BFGS. The default, and what you want unless you know
       otherwise: it builds curvature from the gradients it was going to
       compute anyway
   * - ``cg``
     - Polak-Ribiere conjugate gradient, restarting every ten steps
   * - ``cg-auto``
     - The same method restarting on the Powell-Beale test instead -- when the
       directions stop being usefully conjugate rather than on a schedule
   * - ``sd``
     - Steepest descent. Robust, and slow enough that it is mostly a way of
       finding out whether something else was the problem
   * - ``damped``
     - Damped molecular dynamics. Integrates motion and bleeds energy out
       through friction, so it crosses small barriers a downhill method stops
       at. ``timestep``, ``friction``, ``friction_factor`` and
       ``friction_rising`` tune it
   * - ``nr``
     - Newton-Raphson. Needs a Hessian, and converges to whichever stationary
       point it starts nearest -- a minimum only if it began in that basin
   * - ``prfo``
     - Partitioned rational function optimization: climbs to a **transition
       state** rather than descending to a minimum. Needs a Hessian

.. _opt-transition-states:

Transition states
-----------------

``prfo`` follows one Hessian eigenvector uphill while minimising along every
other, so it needs real curvature rather than the approximation a minimiser
accumulates::

    "keywords": {
      "optimization": {"algorithm": "prfo", "hessian_update": "bofill"}
    }

Where the Hessian comes from depends on the method. Restricted Hartree-Fock has
an analytic one and it is supplied directly. Everything else falls back to
DL-FIND's own two-point finite differences, which costs ``6N`` gradients per
Hessian and reaches the same matrix -- correct, and expensive enough to plan
around on anything but a small system.

``hessian_update`` is what keeps that affordable: an exact Hessian at the start
and an update thereafter, rather than a fresh one per step. ``bofill`` is the
usual choice for a saddle point, because it does not assume the matrix is
positive definite -- which at a transition state it is not.

.. _neb:

Finding the guess: a path between two structures
-------------------------------------------------

P-RFO converges on the saddle nearest where it starts, so the hard part of a
transition state is usually the starting geometry rather than the optimization.
Where reactant and product are both known, the path between them supplies one.

Give the second structure as ``endpoint`` and the run becomes a nudged elastic
band: a chain of images between the two, relaxed together, with springs holding
neighbours apart so the chain cannot slide into either well::

    "keywords": {
      "optimization": {
        "algorithm": "lbfgs",
        "endpoint": "product.xyz",
        "images": 7,
        "neb_endpoints": "frozen"
      }
    }

Only ``endpoint`` turns this on. ``images``, ``neb_spring`` or
``neb_endpoints`` without it is refused rather than run as an ordinary
optimization, because a deck that describes a band and supplies one structure
has asked for something it did not provide.

The two structures must be the same atoms in the same order. This is checked,
and the check is worth more than it sounds: images are interpolated atom by
atom, so a product written with two atoms swapped describes a reaction in which
those atoms trade places. That band relaxes perfectly happily to something with
no chemical meaning.

The implementation is improved-tangent NEB with a climbing image. The climbing
image is the one that matters for a transition state: once the band has settled,
the highest image stops obeying the springs and is driven uphill along the path
instead, so it converges on the barrier top rather than near it. Its geometry is
what you hand to P-RFO.

``neb_endpoints`` decides whether the two ends may move. ``frozen`` is the
default and the usual case, since endpoints are normally already-optimized
minima and relaxing them again spends images to move nothing.
``perpendicular`` lets a slightly-off endpoint settle onto the path without
sliding along it, and ``free`` relaxes them fully.

.. note::

   Images are interpolated on a straight line between the two structures, so
   the interpolation has to be physical. On ammonia inversion it is: the
   nitrogen passes through the plane of the hydrogens, which is the reaction.
   On ``HCN`` to ``HNC`` it is not -- the hydrogen is at one end of the molecule
   in the reactant and the other end in the product, and interpolating drives it
   straight through both nuclei. The middle images then have no SCF and the run
   stops with an energy evaluation failure. Where the straight line is
   unphysical, optimize a rough intermediate by hand and start P-RFO from it
   instead.

Ammonia inversion, ``HF/STO-3G``, seven images: the band converges with the
climbing image at the middle of the path, and the geometry it lands on has all
four atoms coplanar -- the ``D3h`` transition state -- at a barrier of 0.0164
Hartree, or 10.3 kcal/mol.

.. _saddle-search:

Saying that a saddle is what you want
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A P-RFO run converges on the same gradient criterion as a minimisation, so it
tells you it found a stationary point and not what kind. ``target`` says which
kind was wanted, and ``hess_end`` then reads its own result against that::

    "keywords": {
      "optimization": {
        "algorithm": "prfo",
        "hessian_update": "bofill",
        "coordinates": "dlc",
        "target": "saddle",
        "hess_end": true
      }
    }

The curvature test does not change; which count is a pass does. With
``target: "saddle"`` one imaginary frequency is success and is reported as
such, no imaginary frequencies means the search fell off the ridge into a
basin, and two or more means a higher-order stationary point that is not a
transition state. With the default ``target: "minimum"`` those readings invert,
which is the behaviour every existing deck already has.

On success the imaginary mode is broken down by atom::

    one imaginary frequency, 1245.77 cm-1: this is a first-order saddle point
    the imaginary mode moves, most to least:
       atom 3 H   85.7% of the motion
       atom 2 C    7.7% of the motion
       atom 1 N    6.6% of the motion

That is the check no frequency count can make for you. One imaginary frequency
proves a first-order saddle; it does not prove it is *your* saddle. A molecule
has more than one and P-RFO finds whichever is nearest the guess, so the
question is whether the mode moves the atoms whose bond is being made or
broken. Above, on the ``HCN`` to ``HNC`` isomerisation, the hydrogen carries 86%
of the motion, which is the transfer.

``target: "saddle"`` is refused with a downhill algorithm. Steepest descent,
conjugate gradient and L-BFGS will report a minimum however the target is
spelled, so the combination is a deck that cannot be satisfied rather than one
that is satisfied slowly; use ``prfo``, or ``nr`` from a guess already near the
ridge. This is refused whether the algorithm was named or left to the default,
which is L-BFGS -- so a saddle search has to say which optimizer runs it.

.. warning::

   Do not run a saddle search in Cartesian coordinates. P-RFO maximises along
   the lowest Hessian eigenvalue, and in Cartesians the lowest six modes are
   translations and rotations -- so it follows a rotation instead of the
   reaction mode and wanders. On the ``HCN``/``HNC`` saddle, from the same guess
   and the same analytic Hessian, ``cartesian`` and ``hdlc`` both ran 60 steps
   without converging while ``dlc`` converged in four. A deck that asks for both
   is warned; use ``dlc``, or ``hdlc`` for a cluster where plain DLC has no
   bonds to span the fragments.

   ``cartesian`` is the *default* coordinate system, so this is what a saddle
   search that never mentioned coordinates gets. The warning says so, and says
   that the run will spend all of ``max_steps`` before giving up. It is a
   warning rather than a refusal because a Cartesian saddle search can be made
   to work, from a guess close enough that the reaction mode is already the
   lowest eigenvalue, where a downhill algorithm cannot.

.. _opt-constraints:

Holding part of the geometry still
----------------------------------

Two different things, which are worth not confusing.

``frozen_atoms`` removes atoms from the optimization entirely::

    "optimization": {"frozen_atoms": [0, 1, 2]}

Indices are 0-based. A frozen atom contributes no coordinates at all -- it is
not restrained towards its position, it simply is not a variable. It also stops
belonging to a residue, which is DL-FIND's own model rather than a choice made
here, and means ``dlc`` and ``dlc-tc`` cannot express it: pure internals have no
way to hold an atom still, and a deck combining them is refused rather than left
to DL-FIND, which ends the process.

``constraints`` fixes an internal coordinate while the atoms it is measured over
stay free::

    "optimization": {
      "coordinates": "hdlc",
      "constraints": [
        {"type": "bond", "atoms": [0, 4]},
        {"type": "torsion", "atoms": [0, 1, 2, 3]}
      ]
    }

Types are ``bond``, ``angle``, ``torsion``, ``cartesian`` and
``bond-difference``, taking 2, 3, 4, 1 and 3 atoms respectively. The count is
checked against the type, because three atoms under ``"torsion"`` is a
well-formed list and a meaningless constraint.

A constraint needs an internal coordinate system to live in, so ``cartesian`` is
refused -- there DL-FIND has nowhere to put one and would ignore it silently,
converging to the unconstrained answer without saying so.

User-supplied connections (DL-FIND's ``nconn``) are deliberately not exposed.
DL-FIND reads them from an offset computed as ``nat + nz + ncons`` while it
reads the block after them from ``nat + nz + 5*ncons + 2*nconn``; the two
disagree about how long the constraint block is, so a deck using both would have
its connections read out of the middle of its constraints.

.. _hess-end:

Checking that it is a minimum
-----------------------------

An optimization converges on the gradient, and a vanishing gradient is a
*stationary point*: a saddle satisfies it exactly as well as a minimum does.
Nothing in the run distinguishes them, so "converged" names the condition that
was met rather than the thing that was found. The second derivatives are what
settle it, and ``hess_end`` asks for them::

    "keywords": {
      "optimization": {"hess_end": true}
    }

A Hessian is computed at the converged geometry and the frequencies are printed
as they are for a ``Hessian`` deck, followed by the one line that answers the
question: no imaginary frequencies means a minimum, exactly one means a
first-order saddle point, more than one means neither. A mode is called
imaginary below ``-10`` cm-1, the same threshold the frequency summary uses to
separate vibrations from translation and rotation.

One thing to expect in the output. The thermochemistry block a few lines above
the verdict prints its own ``Imaginary freqs`` count, and the two can disagree:
thermochemistry treats any frequency below zero as imaginary, so the projected
rotational modes -- which land at around ``-1e-5`` cm-1 rather than exactly zero
-- are counted there and not here. On a genuine minimum you will see something
like ``Imaginary freqs: 3 (skipped)`` directly above ``no imaginary
frequencies``. The verdict uses the ``-10`` cm-1 window and is the one answering
the question; the frequency table above both settles any doubt.

Restricted Hartree-Fock only for now, and refused at the start of the run rather
than at the end -- the end is the one place where the news is expensive, because
the optimization has already been paid for. For anything else, run the
optimization and then a separate ``Hessian`` deck on the geometry it writes.

It runs only on a converged optimization. Curvature at a geometry the optimizer
was still moving away from describes that geometry and not the minimum, and
printing it under a "did not converge" line invites it to be read as the
minimum's.

The cost is one Hessian on top of the optimization -- analytic for restricted
Hartree-Fock in a libfint build, and central differences of gradients otherwise,
which is ``6N`` gradients. That is why it is off by default.

.. _optimizer-build:

Building with the optimizer
---------------------------

Off by default:

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_DLFIND=ON
   cmake --build build -j

The default is about licensing rather than size. DL-FIND is LGPL-3 and metalquicha
is MIT, so libdlfind is fetched and linked as a **shared** library: the ``.so``
stays LGPL, this program stays MIT, and anyone receiving a binary keeps the right
to relink it against their own copy. Compiling those sources into ``libmqc.a``
instead would carry the relinking obligation into every binary shipped from here,
so the build deliberately does not do that. Turning the flag on is a choice the
person building makes.

The configure step clones libdlfind. For a machine without a network -- most
cluster compute nodes -- point it at a local copy:

.. code-block:: bash

   cmake -B build -DMQC_ENABLE_DLFIND=ON \
         -DFETCHCONTENT_SOURCE_DIR_LIBDLFIND=/path/to/libdlfind

To find out whether a binary has it, ask:

.. code-block:: bash

   ./mqc --version

which prints a ``features:`` line naming the optional backends that were built in.
``validation/run_validation.py`` reads that line, so a build without the optimizer
skips the optimization cases rather than failing them.

Acknowledgements
----------------

**DL-FIND** is the optimizer. It is the work of Johannes Kästner, Tom Keal and
co-workers, distributed under the LGPL, and is the geometry optimization library
of `ChemShell <https://www.chemshell.org>`_. If you publish a structure optimized
with this driver, cite it:

   J. Kästner, J. M. Carr, T. W. Keal, W. Thiel, A. Wander and P. Sherwood,
   *DL-FIND: An Open-Source Geometry Optimizer for Atomistic Simulations*,
   J. Phys. Chem. A **113**, 11856 (2009).
   `doi:10.1021/jp9028968 <https://doi.org/10.1021/jp9028968>`_

**libdlfind** is the C and Python interface that makes DL-FIND callable from
outside ChemShell, by Kjell Jorner, also LGPL. This program uses its C entry point
``api_dl_find``, which takes the optimizer's callbacks as function pointers -- so
nothing of libdlfind's Fortran needs to cross into this build, and the two projects
stay at arm's length.

Neither library is modified here. What lives in this repository is the bridge
between them and metalquicha's own calculations, in ``backends/dlfind/``.

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
   * - libcint (CPU ``hf``, ``dft``, ``mp2``, ``ccsd``)
     - **no**
     - Refused: the CPU backend computes energies only
   * - cuEST (GPU ``hf``, ``dft``)
     - yes
     - Expected to work; needs an sm_80 card

A method that cannot produce a gradient is found out before the optimizer starts,
by taking one gradient up front, and refused in a sentence:

.. code-block:: text

   The initial gradient check failed: the CPU backend has no gradients yet;
   run an energy, or build with cuEST
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

``cartesian`` is the default and always works. ``hdlc`` builds internal coordinates
within each monomer and keeps Cartesians between them, which is the right shape for
a molecular cluster; the residues it needs are this program's own monomers, so
nothing has to be perceived.

``dlc`` puts every atom in one residue, which suits a single connected molecule and
**fails on a cluster** -- there is no connected internal-coordinate system spanning
molecules that are not bonded to each other, and DL-FIND says so
(``cyclic failure at residue 1``).

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
     - ``cartesian``, ``hdlc`` or ``dlc``. ``coordinate_system`` is accepted too
   * - ``algorithm``
     - ``lbfgs``
     - ``lbfgs``, ``cg`` or ``sd``. ``optimizer`` is accepted too
   * - ``lbfgs_memory``
     - engine
     - Steps of curvature history L-BFGS keeps
   * - ``freeze_terms``
     - ``true``
     - Fix the MBE term list at the starting geometry. See :ref:`frozen-terms`
   * - ``trajectory``
     - ``true``
     - Record every accepted geometry
   * - ``print_level``
     - follows log
     - How much DL-FIND itself prints

Where a default is given as "engine", the setting is left to DL-FIND rather than
overridden with a second number chosen here.

``prfo`` parses as an algorithm but is refused: it searches for a transition state
rather than a minimum and needs a Hessian this interface does not yet supply.

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

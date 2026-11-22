Nuclear magnetic shielding
==========================

The shielding of a nucleus is how much the electrons around it weaken an
applied magnetic field at that nucleus. It is what a chemical shift is the
difference of: :math:`\delta = \sigma_\mathrm{ref} - \sigma`, with
:math:`\sigma_\mathrm{ref}` the same nucleus in a reference compound computed
the same way.

.. code-block:: json

   "model": {"method": "hf", "basis": "cc-pVDZ"},
   "driver": "Energy",
   "properties": {"nmr": {}}

``properties.nmr`` is a request by its presence and takes no settings. Every
nucleus is shielded. Like everything under ``properties`` it changes no energy.
It runs for closed-shell Hartree-Fock and for Kohn-Sham with an LDA, a GGA or
a global hybrid of either:

.. code-block:: json

   "model": {"method": "dft", "functional": "b3lyp", "basis": "cc-pVTZ"},
   "driver": "Energy",
   "properties": {"nmr": {}}

What it computes
----------------

The shielding tensor is a second derivative of the energy, with respect to the
applied field :math:`\mathbf{B}` and the magnetic moment
:math:`\boldsymbol{\mu}_N` of the nucleus:

.. math::

   \sigma^N_{xy} = \frac{\partial^2 E}{\partial B_x\,\partial \mu_{N,y}}

A magnetic field enters through a vector potential, which has an origin, and
in a finite basis the result depends on where that origin is put. The basis
functions here are **gauge-including atomic orbitals** (GIAOs, or London
orbitals): each carries a phase that moves its own origin to its own centre.
The shielding is then independent of the origin, and of where the molecule sits
in space, at any basis size.

The tensor has two parts. The *diamagnetic* part is an expectation value over
the unperturbed density. The *paramagnetic* part is the response of the
orbitals to the field, contracted with the nucleus's own operator, and costs one
coupled-perturbed solve for each of the three field directions. The field is an
imaginary perturbation, so the density it induces is antisymmetric: it has no
Coulomb response, and in Hartree-Fock only exchange couples the equations.

Only the sum of the two parts is a physical quantity. How it divides between
them depends on the gauge, and the split here is the GIAO one.

Under a functional
------------------

Two things change, and one of them makes a pure functional *cheaper* than
Hartree-Fock.

**The exchange-correlation potential does not respond to the field.** It is a
functional of the density, and a magnetic field induces a current, not a
change of density, at first order. So the potential contributes to the
first-order Fock matrix only through the gauge phase of the basis functions:
for a function centred at :math:`\mathbf{R}` the field derivative at a point
:math:`\mathbf{r}` is :math:`\tfrac{i}{2}(\mathbf{R}\times\mathbf{r})\chi`,
and the matrix is the potential between one function and another's
derivative, less its transpose, integrated on the functional's grid. A GGA
adds the gradient of that derivative.

**Exact exchange is scaled by the fraction the functional keeps**, in the
first-order Fock matrix and in the response alike. For a pure functional that
fraction is zero, nothing couples the three equations, and they are solved by
the orbital-energy gaps in one step with no response iterations at all.

That a functional of the density alone ignores the current is an
approximation of the functional, not of this code: it is what every
conventional functional does, and it is why shieldings from pure functionals
are too deshielded for nuclei with low-lying excited states.

What it reports
---------------

For each nucleus, the isotropic shielding (a third of the trace) and the
anisotropy, in parts per million:

.. code-block:: text

     ==== nuclear magnetic shielding (GIAO), ppm ====

        atom       isotropic    anisotropy
       O 1           346.7137       30.6279
       H 2            31.1325       19.6033
       H 3            31.1325       19.6033

That is water at Hartree-Fock/cc-pVDZ. The anisotropy is
:math:`\sigma_{33} - (\sigma_{11} + \sigma_{22})/2`, from the principal values
of the symmetrised tensor with :math:`\sigma_{33}` the one furthest from the
isotropic value. At the ``verbose`` logger level the full tensors are printed
too.

The output file carries everything under ``nmr_shielding``: for each atom its
number ``atom`` (from one) and ``element``, ``isotropic``, ``anisotropy``, the ``tensor`` as three rows, and its
``diamagnetic`` and ``paramagnetic`` parts.

**The tensor is not symmetric.** Its first index is the field's direction and
its second the nuclear moment's; a row of ``tensor`` is one field direction.
For a hydrogen of water the two in-plane off-diagonal elements differ by about
2 ppm. Only the symmetric part reaches a spectrum, which is why the summary
symmetrises, and the rows are written as computed.

Validation
----------

The reference is ``pyscf.prop.nmr.RHF``, fed this repository's basis sets.

``test_mqc_czt_nmr`` pins all nine elements of all three tensors of water in
6-31G to it, to the six decimals PySCF prints, and checks that moving the
molecule by several bohr changes no element by more than :math:`10^{-6}` ppm
(measured: :math:`3\times10^{-8}`). Through the program, the cases in the table
below agree with PySCF element by element.

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Molecule
     - Basis
     - Largest difference in any tensor element (ppm)
   * - Water
     - cc-pVDZ
     - 1.2e-6
   * - Ammonia
     - cc-pVDZ
     - 8.2e-7
   * - Methane
     - cc-pVDZ
     - 9.8e-7
   * - HCN
     - cc-pVDZ
     - 2.7e-7
   * - Water
     - 6-31G* (Cartesian d)
     - 7.2e-7
   * - HF
     - cc-pVTZ
     - 3.8e-7
   * - Ethane
     - cc-pVDZ
     - 1.5e-7

Under a functional the reference is ``pyscf.prop.nmr.rks``. Its quadrature
loop fails under the current PySCF, so the numbers come from the same formulas
run in one pass, and that copy reproduces the values in the module's own
self-test -- hydrogen fluoride at LDA/6-31G and at the B3LYP of its day in
cc-pVDZ -- to :math:`10^{-4}` ppm. ``test_mqc_czt_nmr`` pins all three tensors
of water in 6-31G at PBE and at B3LYP to it, agreeing to
:math:`1.2\times10^{-6}` ppm, and checks that moving the molecule changes no
PBE element by more than :math:`10^{-5}` ppm (measured: :math:`8\times10^{-12}`).
LDA was checked through the program to the same.

In a continuum solvent
----------------------

With ``keywords.pcm`` on the deck (see :doc:`continuum_solvation`) the
shieldings are those of the solvated molecule, for C-PCM and IEF-PCM, under
Hartree-Fock or a functional.

**The surface charges do not respond to the field.** They answer the
molecule's charge density, and a magnetic field induces a current and no
change of density at first order, so the charges stay as the converged
reference left them. Their potential is in the Fock matrix, though, and
therefore its derivative through the gauge factor of the basis functions has
to be in the first-order one: the same term the nuclei contribute, with each
smeared surface charge standing as a nucleus. Leaving it out makes the
shielding depend on where the molecule sits -- by 87 ppm for the oxygen of a
water molecule moved a few bohr -- so its presence is checked that way, since
no other program here computes a shielding in a continuum to compare with.

For water at Hartree-Fock/6-31G the continuum at water's dielectric constant
raises the oxygen shielding by 11.7 ppm and lowers the hydrogens' by 0.9.
**A continuum is bulk polarity and nothing else.** The measured gas-to-liquid
change of water's hydrogens is several times larger and the oxygen's has the
opposite sign, because both are dominated by hydrogen bonds to particular
neighbours, which a continuum does not contain. For a hydroxyl or amine
hydrogen in a protic solvent, a few explicit solvent molecules matter more
than the continuum around them.

In a fragmented run
-------------------

With ``keywords.fragmentation`` on the deck the shieldings are assembled over
the many-body expansion, with disjoint fragments or overlapping ones, and
written in the same ``nmr_shielding`` section.

A shielding belongs to one atom, so its expansion runs over the terms that
contain that atom: its value in its own monomer, plus the change it sees in
each pair containing that monomer, and so on. That is the expansion's ordinary
weights applied atom by atom, the way a gradient is assembled, and it is exact
when the expansion reaches the whole system. Gauge-including orbitals are what
make it legitimate to add numbers from different terms: each term's shieldings
do not depend on where the term sits. A cap hydrogen's shielding belongs to no
atom of the system and is dropped.

How fast it converges, against the unfragmented calculation at Hartree-Fock in
6-31G, as the largest error in an isotropic shielding:

.. list-table::
   :header-rows: 1
   :widths: 34 16 25 25

   * - System
     - Level
     - Worst hydrogen (ppm)
     - Worst oxygen (ppm)
   * - Three waters, stacked
     - 1
     - 0.63
     - 6.3
   * -
     - 2
     - 0.044
     - 0.54
   * -
     - 3
     - 5e-13
     - 5e-12
   * - Six waters, the prism
     - 2
     - 0.45
     - 10.6
   * -
     - 3
     - 0.045
     - 1.0

The third level of the trimer is the whole system, and reproduces it to the
solver's precision. Run on four MPI ranks, the trimer at levels 2 and 3
and the prism at level 2 give the serial tensors to :math:`10^{-10}` ppm. **Each level buys about a factor of ten, and oxygen is
ten to twenty times harder than hydrogen**: in a hydrogen-bonded cluster a
proton spectrum is usable at the second level and an oxygen one is not. The
slow part is the polarisation of each molecule by all the others at once,
which an expansion without embedding has to rebuild term by term.

Across covalent bonds
~~~~~~~~~~~~~~~~~~~~~

A molecule that has to be cut through its bonds -- a chain, a peptide -- is
fragmented with hydrogen caps (``connectivity`` on the molecule), and the
shieldings are assembled from the capped fragments. A cap's own tensor is
dropped. What a cap does to its neighbours cannot be dropped: a carbon that had
a carbon beside it and now has a hydrogen is shielded differently by tens of
ppm, and its neighbour by several.

**Overlapping fragments are what make that cancel.** With
``"method": "gmbe"`` the fragments share atoms, and the expansion's terms are
the fragments, less their intersections, by inclusion-exclusion. An
intersection carries its caps exactly where the two fragments it joins carry
theirs, so the disturbance from each cap enters an atom's shielding once with
each sign, and what survives is how two caps at opposite ends interfere. An
atom's shielding is the sum over the terms that contain it of its value there
times the term's coefficient; those coefficients add to one for every atom.

.. code-block:: json

   "molecules": [{
     "xyz": "decane.xyz",
     "connectivity": [[0, 4, 1], [4, 7, 1], [7, 10, 1], [10, 13, 1]],
     "fragments": [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                   [7, 8, 9, 10, 11, 12, 13, 14, 15]],
     "fragment_charges": [0, 0],
     "fragment_multiplicities": [1, 1]
   }],
   "keywords": {"fragmentation": {"method": "gmbe", "level": 1}},
   "properties": {"nmr": {}}

(the first two of five fragments of three carbons, sharing one). Decane at Hartree-Fock in 6-31G,
cut into windows of so many carbons that advance by so many, against the
unfragmented molecule, as the largest error in an isotropic shielding:

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 21 21

   * - Fragments
     - Level
     - Largest term
     - Worst hydrogen (ppm)
     - Worst carbon (ppm)
   * - 3 carbons, advancing 2
     - 1
     - 3 carbons
     - 1.05
     - 28.6
   * - 4 carbons, advancing 2
     - 1
     - 4
     - 0.56
     - 10.6
   * - 5 carbons, advancing 2
     - 1
     - 5
     - 0.28
     - 1.5
   * - 6 carbons, advancing 2
     - 1
     - 6
     - 0.046
     - 0.41
   * - 4 carbons, advancing 1
     - 1
     - 4
     - 0.28
     - 2.7
   * - 5 carbons, advancing 1
     - 1
     - 5
     - 0.046
     - 0.90
   * - 3 carbons, advancing 2
     - 2
     - 5
     - 0.93
     - 9.1
   * - 3 carbons, advancing 2
     - 3
     - 7
     - 0.62
     - 4.0
   * - 4 carbons, advancing 2
     - 2
     - 6
     - 0.006
     - 0.058
   * - 2 carbons, disjoint (MBE)
     - 2
     - 4
     - 0.59
     - 18.4
   * - 2 carbons, disjoint (MBE)
     - 3
     - 6
     - 0.057
     - 0.72

**What matters is how far every atom is from a cap in the term that keeps it
whole, and how much the fragments share.** Windows that overlap by one carbon
leave that carbon beside a cap in every term containing it, and it is wrong by
tens of ppm; an overlap of four carbons brings carbon under half a ppm at the
first level. The second level, whose terms are unions of two fragments, is
ten times better again. **A higher level does not repair a thin overlap**:
windows sharing one carbon are still 4 ppm off at the third level, with terms
of seven carbons, where windows sharing two reach 0.06 at the second with
terms of six. For the same largest term the overlapping expansion
beats the disjoint one: its intersections put a cap's error in with both
signs, where a disjoint expansion can only outgrow it. At B3LYP the errors
are the same size (0.039 and 0.31 ppm for six-carbon windows advancing two).

Triglycine, cut at three backbone bonds into three overlapping fragments,
keeps hydrogen within 0.09 ppm and carbon within 0.5 at the first level, with
nitrogen at 1.4 and oxygen at 4.6: the lone-pair nuclei are the slow ones, as
in a hydrogen-bonded cluster. Its second level reaches the whole molecule and
reproduces the unfragmented tensors to :math:`4\times10^{-10}` ppm, which is
the check on the assembly that needs no reference. On four MPI ranks the
tensors are the serial ones to :math:`3\times10^{-10}` ppm.

One rule was tried against this and lost: reading each atom from the single
fragment where it is furthest from any cap, summing nothing. For six-carbon
windows that leaves carbon 5.5 ppm off where inclusion-exclusion leaves 0.41.
A cap two bonds away is still worth several ppm to a carbon, and only the
cancellation removes it.

What is not assembled:

- The spin-spin couplings, under any fragmentation: the deck is refused. A
  coupling joins two nuclei that may sit in different fragments; ask for them
  on an unfragmented run.
- With FMO, EE-MBE or EFMO: the deck is refused.
- Under counterpoise, where the weights do not know a ghosted term from an
  ordinary one: a warning, and no shieldings.
- For an interaction-energy run, which has no whole system to report on.
- With an embedding field, whose potential would need the same gauge term a
  continuum's surface charges do. It is not built for the embedding charges.

Limits
------

- **Closed shells with exact integrals only.** An open shell, a correlated
  method and a density-fitted reference are each refused by name.
- **LDA, GGA and their global hybrids.** A meta-GGA is refused, since its
  kinetic energy density depends on the gauge and that dependence is not
  built; so are range-separated hybrids and functionals with non-local
  correlation.
- **No effective core potentials.**
- **A shielding costs a few Fock builds beyond the SCF.** The two-electron
  GIAO term is computed over a quarter of the shell quartets and screened, and
  the response to the field is solved from a coupling matrix built once by
  transforming the integrals. That matrix is as large as the square of the
  number of occupied-virtual pairs; where it does not fit in memory the
  response recomputes the integrals on each of its dozen or so passes
  instead, and is then the slowest step. ``system.memory_gb`` sets what the
  run may plan on. The log gives the time of each phase.
- **Hartree-Fock shieldings are a starting point, not a prediction.** Electron
  correlation changes them by tens of ppm for nuclei with lone pairs or multiple
  bonds, which is what a functional is for, and a basis built for energies
  converges them slowly. Chemical shifts, being differences, are more
  forgiving than the shieldings themselves.
- **Shieldings only.** No chemical shifts against a named reference, no
  spin-spin couplings, no magnetisability.
- **CPU backend.** The GPU and terco backends warn and skip it.

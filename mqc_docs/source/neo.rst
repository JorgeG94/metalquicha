Quantum Nuclei (NEO)
====================

The **nuclear-electronic orbital** (NEO) method treats chosen nuclei -- protons,
so far -- quantum mechanically, on the same footing as the electrons. A quantised
proton gets orbitals of its own, expanded in a Gaussian basis centred where the
classical nucleus stood, and is solved self-consistently with the electrons. The
wavefunction then carries the proton's zero-point motion and delocalisation
directly, rather than as a harmonic correction added afterwards.

This is the method of Hammes-Schiffer and co-workers. The reference
implementation this one is validated against is the ``pyscf/neo`` module of
Yang Yang's PySCF fork (https://github.com/theorychemyang/pyscf).

Running one
-----------

.. code-block:: json

   {
     "schema": {"name": "hcn_neo", "version": "1.0"},
     "molecules": [{
       "xyz": "hcn.xyz",
       "molecular_charge": 0,
       "molecular_multiplicity": 1
     }],
     "model": {"method": "hf", "basis": "cc-pvdz"},
     "keywords": {
       "neo": {"quantum_nuclei": ["H"], "nuclear_basis": "pb4-d"}
     },
     "driver": "Energy"
   }

``keywords.neo.quantum_nuclei`` is either a list of element symbols, which
quantises every atom of those elements, or a list of 0-based atom indices.
``nuclear_basis`` names the proton basis; the even-tempered PB4-D through PB6-H
sets of Yu, Pavosevic and Hammes-Schiffer (J. Chem. Phys. 152, 244123, 2020)
ship under ``basis_sets/neo/`` and PB4-D is the default.

The energy reported is the NEO-HF total: electrons and classical nuclei, the
quantum nuclei's kinetic and potential energy, and every coupling between them.
With the logger at ``info`` the run prints the macro-iteration table and the
energy broken into those pieces, together with each proton's orbital energy.

What is computed
----------------

Each quantum nucleus is a single particle in its own component, so there is no
exchange and no self-Coulomb inside a component. With ``D_e`` the electronic
density and ``D_p`` a proton's,

.. code-block:: text

   E = E_e[D_e; classical nuclei] + sum_p <D_p| T/m_p + V_classical |D_p>
     - sum_p J(D_e, D_p)          + sum_{p<q} J(D_p, D_q)

The electrons see a quantised nucleus only through its density: its charge is
removed from the electronic Hamiltonian and from the classical nuclear
repulsion. The cross Coulomb terms are ordinary four-centre integrals over a
combined basis, so no new integral code is involved. The two problems are
coupled by macro-iteration -- the electrons converge in the field of the
current proton densities, then every proton is re-solved in the field of the
new electrons -- until the total energy and the proton densities stop moving.

The proton mass follows PySCF-NEO: the most common isotope's atomic mass less
one electron, 1836.15265 electron masses.

Limits, for now
---------------

* Closed-shell Hartree-Fock only. ``model.method`` other than ``hf`` is refused.
* Only hydrogen can be quantised. Naming a heavier atom is refused.
* Energies only: no gradients, and no NEO-DFT electron-proton correlation
  functional yet.
* The calculation is whole-system; fragmentation keywords are ignored.

Validation
----------

HCN with the proton quantised, cc-pVDZ and PB4-D, reproduces PySCF-NEO's
``test_hf.py`` energy of -92.8437063566 Hartree to 1e-9, and the proton's
orbital energy and one-body energy to 1e-7. The unit test
``test_mqc_czt_neo`` pins this.

Adding New JSON Output Fields
==============================

This guide explains how to add new fields to the JSON output in metalquicha.

Overview
--------

JSON output is centralized in a single location. All calculation workflows populate
a ``json_output_data_t`` container, and the driver calls ``write_json_output()`` once
at the end to write everything to disk.

**Key files:**

- ``src/core/mqc_json_output_types.f90`` - Container type definition
- ``src/io/mqc_json_writer.f90`` - JSON writing logic

Step-by-Step Guide
------------------

1. Add the Field to the Container Type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Edit ``src/core/mqc_json_output_types.f90`` and add your new field to the
``json_output_data_t`` type:

.. code-block:: fortran

   type :: json_output_data_t
      ! ... existing fields ...

      ! Your new field
      real(dp), allocatable :: mulliken_charges(:)
      logical :: has_mulliken_charges = .false.

      ! ... rest of type ...
   end type

**Conventions:**

- Use ``has_<field>`` boolean flags for optional data
- Initialize booleans to ``.false.`` in the type definition
- Use allocatable arrays for variable-size data

2. Update the Cleanup Procedures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the same file, update ``json_output_data_destroy`` and ``json_output_data_reset``:

.. code-block:: fortran

   subroutine json_output_data_destroy(self)
      class(json_output_data_t), intent(inout) :: self
      ! ... existing deallocations ...
      if (allocated(self%mulliken_charges)) deallocate(self%mulliken_charges)
   end subroutine

   subroutine json_output_data_reset(self)
      class(json_output_data_t), intent(inout) :: self
      call self%destroy()
      ! ... existing resets ...
      self%has_mulliken_charges = .false.
   end subroutine

3. Populate the Field in Your Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the workflow module that computes your data (e.g., ``mqc_unfragmented_workflow.f90``),
populate the field when ``json_data`` is present:

.. code-block:: fortran

   if (present(json_data)) then
      ! Populate your new field
      allocate(json_data%mulliken_charges(n_atoms))
      json_data%mulliken_charges = calculated_charges
      json_data%has_mulliken_charges = .true.
   end if

4. Add JSON Writing Logic
^^^^^^^^^^^^^^^^^^^^^^^^^

Edit ``src/io/mqc_json_writer.f90`` and add the JSON output in the appropriate
``write_*_impl`` subroutine:

.. code-block:: fortran

   subroutine write_unfragmented_json_impl(data)
      ! ... existing code ...

      ! Add your field to the JSON
      if (data%has_mulliken_charges) then
         call json%add(main_obj, 'mulliken_charges', data%mulliken_charges)
      end if

      ! ... rest of subroutine ...
   end subroutine

**Which subroutine to modify:**

- ``write_unfragmented_json_impl`` - For unfragmented calculations
- ``write_mbe_breakdown_json_impl`` - For MBE fragmented calculations
- ``write_gmbe_pie_json_impl`` - For GMBE calculations
- ``write_vibrational_json_impl`` - For vibrational/thermochemistry data

If your field applies to multiple calculation types, add it to each relevant subroutine.

5. Build and Test
^^^^^^^^^^^^^^^^^

.. code-block:: bash

   cmake --build build -j
   ./build/mqc your_input.json
   cat output_your_input.json

Example: Adding Orbital Energies
--------------------------------

Here's a complete example of adding HOMO/LUMO energies:

**1. In mqc_json_output_types.f90:**

.. code-block:: fortran

   type :: json_output_data_t
      ! ... existing fields ...
      real(dp) :: homo_energy = 0.0_dp
      real(dp) :: lumo_energy = 0.0_dp
      logical :: has_orbital_energies = .false.
   end type

**2. In your workflow (e.g., mqc_unfragmented_workflow.f90):**

.. code-block:: fortran

   if (present(json_data)) then
      json_data%homo_energy = result%homo
      json_data%lumo_energy = result%lumo
      json_data%has_orbital_energies = .true.
   end if

**3. In mqc_json_writer.f90:**

.. code-block:: fortran

   if (data%has_orbital_energies) then
      call json%create_object(orbitals_obj, 'orbital_energies')
      call json%add(main_obj, orbitals_obj)
      call json%add(orbitals_obj, 'homo_hartree', data%homo_energy)
      call json%add(orbitals_obj, 'lumo_hartree', data%lumo_energy)
      call json%add(orbitals_obj, 'gap_hartree', data%lumo_energy - data%homo_energy)
   end if

Architecture Notes
------------------

The JSON output flow:

.. code-block:: text

   run_calculation()
       |
       +-- workflow populates json_data fields
       |
       +-- write_json_output(json_data)  <-- SINGLE write point
              |
              +-- dispatches based on output_mode
              |
              +-- calls appropriate write_*_impl()

**Output modes** (set automatically by workflows):

- ``OUTPUT_MODE_UNFRAGMENTED`` - Single molecule calculation
- ``OUTPUT_MODE_MBE`` - Many-Body Expansion
- ``OUTPUT_MODE_GMBE_PIE`` - Generalized MBE with PIE

The ``output_mode`` determines which ``write_*_impl`` function is called, but you
can add shared fields (like orbital energies) to multiple implementations if needed.

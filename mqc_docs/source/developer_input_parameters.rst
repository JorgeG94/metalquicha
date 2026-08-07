Adding New Input File Parameters
==================================

This guide explains how to add new parameters to metalquicha's input file format.

Overview
--------

Input parameters flow through three stages:

1. **Parsing** - ``mqc_json_config_reader.f90`` reads the JSON deck into ``mqc_config_t``
2. **Adaptation** - ``mqc_config_adapter.f90`` converts to ``driver_config_t``
3. **Usage** - Driver and workflows use the configuration

**Key files:**

- ``src/io/mqc_config_parser.f90`` - Parser and ``mqc_config_t`` type
- ``src/io/mqc_config_adapter.f90`` - Adapter and ``driver_config_t`` type
- ``src/core/mqc_calculation_keywords.f90`` - Keyword group types (hessian, aimd, scf)

Step-by-Step Guide
------------------

1. Add Field to mqc_config_t
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Edit ``src/io/mqc_config_parser.f90`` and add your field to ``mqc_config_t``:

.. code-block:: fortran

   type :: mqc_config_t
      ! ... existing fields ...

      ! Your new field with a sensible default
      real(dp) :: my_new_parameter = 1.0_dp

   end type mqc_config_t

**Conventions:**

- Always provide a default value
- Group related fields together with comments
- Use appropriate types: ``integer``, ``real(dp)``, ``logical``, ``character(len=:), allocatable``

2. Parse the Parameter
^^^^^^^^^^^^^^^^^^^^^^

Find or create the appropriate section parser in the same file. Parameters are
organized by section (``%model``, ``%driver``, ``%hessian``, etc.).

For example, to add a parameter to the ``%driver`` section:

.. code-block:: fortran

   subroutine parse_driver_section(unit, config, error)
      ! ... existing parsing code ...

      do
         read(unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         line = strip_comment(line)
         if (len_trim(line) == 0) cycle
         if (trim(line) == 'end') exit

         ! Parse key = value
         eq_pos = index(line, '=')
         if (eq_pos > 0) then
            key = to_lower(trim(adjustl(line(1:eq_pos-1))))
            value = trim(adjustl(line(eq_pos+1:)))

            select case (key)
            ! ... existing cases ...

            case ('my_new_parameter')
               read(value, *, iostat=io_stat) config%my_new_parameter
               if (io_stat /= 0) then
                  call error%set(ERROR_PARSE, "Invalid my_new_parameter value: "//trim(value))
                  return
               end if

            case default
               call error%set(ERROR_PARSE, "Unknown driver option: "//trim(key))
               return
            end select
         end if
      end do
   end subroutine

3. Add to driver_config_t (if needed at runtime)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If your parameter is needed during calculation (not just for parsing), add it to
``driver_config_t`` in ``src/io/mqc_config_adapter.f90``:

.. code-block:: fortran

   type :: driver_config_t
      ! ... existing fields ...

      real(dp) :: my_new_parameter = 1.0_dp

   end type driver_config_t

Then copy it in ``config_to_driver``:

.. code-block:: fortran

   subroutine config_to_driver(mqc_config, driver_config, molecule_index)
      ! ... existing code ...

      driver_config%my_new_parameter = mqc_config%my_new_parameter

   end subroutine

4. Use in Keyword Groups (for structured settings)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For related parameters (like Hessian settings), use the keyword group types in
``src/core/mqc_calculation_keywords.f90``:

.. code-block:: fortran

   type :: hessian_keywords_t
      real(dp) :: displacement = 0.001_dp
      real(dp) :: temperature = 298.15_dp
      real(dp) :: pressure = 1.0_dp
      ! Add your new field here
      logical :: project_rotations = .true.
   end type

These are accessed via ``driver_config%hessian%project_rotations``.

5. Update Documentation
^^^^^^^^^^^^^^^^^^^^^^^

Add your parameter to ``mqc_docs/source/input_files.rst`` in the appropriate section.

Example: Adding a New Section
-----------------------------

To add an entirely new section (e.g., ``%optimization``):

**1. Add fields to mqc_config_t:**

.. code-block:: fortran

   type :: mqc_config_t
      ! ... existing fields ...

      ! Optimization settings
      integer :: opt_maxiter = 100
      real(dp) :: opt_convergence = 1.0e-6_dp
      character(len=:), allocatable :: opt_algorithm

   end type

**2. Create a section parser:**

.. code-block:: fortran

   subroutine parse_optimization_section(unit, config, error)
      integer, intent(in) :: unit
      type(mqc_config_t), intent(inout) :: config
      type(error_t), intent(out) :: error

      character(len=MAX_LINE_LEN) :: line
      character(len=:), allocatable :: key, value
      integer :: io_stat, eq_pos

      do
         read(unit, '(A)', iostat=io_stat) line
         if (io_stat /= 0) exit

         line = strip_comment(line)
         if (len_trim(line) == 0) cycle
         if (trim(line) == 'end') exit

         eq_pos = index(line, '=')
         if (eq_pos > 0) then
            key = to_lower(trim(adjustl(line(1:eq_pos-1))))
            value = trim(adjustl(line(eq_pos+1:)))

            select case (key)
            case ('maxiter')
               read(value, *, iostat=io_stat) config%opt_maxiter
            case ('convergence')
               read(value, *, iostat=io_stat) config%opt_convergence
            case ('algorithm')
               config%opt_algorithm = trim(value)
            case default
               call error%set(ERROR_PARSE, "Unknown optimization option: "//trim(key))
               return
            end select
         end if
      end do
   end subroutine

**3. Register the section in read_mqc_file:**

.. code-block:: fortran

   subroutine read_mqc_file(filename, config, error)
      ! ... existing code ...

      select case (trim(line))
      ! ... existing cases ...

      case ('%optimization')
         call parse_optimization_section(unit, config, parse_error)
         if (parse_error%has_error()) then
            error = parse_error
            close(unit)
            return
         end if

      case default
         ! Unknown section handling
      end select

   end subroutine

**4. Create keyword type if needed:**

In ``src/core/mqc_calculation_keywords.f90``:

.. code-block:: fortran

   type :: optimization_keywords_t
      integer :: maxiter = 100
      real(dp) :: convergence = 1.0e-6_dp
      character(len=32) :: algorithm = "bfgs"
   end type optimization_keywords_t

Input File Format Reference
---------------------------

Sections use the format:

.. code-block:: text

   %section_name
   key1 = value1
   key2 = value2
   end

**Supported value types:**

- Integers: ``maxiter = 100``
- Reals: ``tolerance = 1.0e-6``
- Strings: ``method = gfn2`` (no quotes needed)
- Booleans: ``use_symmetry = true`` or ``use_symmetry = false``
- Arrays: ``cutoffs = 5.0, 4.0, 3.0`` (comma-separated)

**Comments:** Use ``#`` or ``!`` for comments:

.. code-block:: text

   %driver
   calc_type = gradient  # This is a comment
   end

Testing Your Changes
--------------------

1. Create a test input file with your new parameter
2. Build and run:

   .. code-block:: bash

      cmake --build build -j
      ./build/mqc test_input.json

3. Add a unit test in ``test/`` if appropriate

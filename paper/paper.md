---
title: 'Metalquicha: A Fortran program for fragmented quanutm chemistry calculations'
tags:
  - Fortran
  - Electronic structure theory
  - Quantum Chemistry
  - Fragmentation methods
authors:
  - name: Jorge Luis Galvez Vallejo
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1"
affiliations:
 - name: National Computational Infrastructure, Canberra, Australia
   index: 1
   ror: 00hx57361
date: 1 January 2026
bibliography: paper.bib

---

# Summary

Metalquicha is a Fortran library for massively parallel, fragmented quantum
chemistry calculations with emphasis in modularity, reusability, and ease
of adoption. It supports the Generalized Many Body Expansion (GMBE) with
and without overlapping fragments and can compute scalar and matrix-like
GMBE quantities, such as energies, gradients, and hessians using xtb
as a backend.


# Statement of need

Quantum chemistry methods are known for being computationally expensive and
impractical for realistic chemical systems, such as proteins, enzymes, and
materials [@gordon_fragmentation_2012]. Fragmentation methods offer a way to address the scalability problem
and additionally expose a naively parallel problem.

The Fortran programming language has been extensively used in the quanutm
chemistry community, a number of legacy packages and new ones exist in the
literature [@zahariev_general_2023; @apra_nwchem_2020; @manathunga_quanutm_2023; @bannwarth_extended_2021; @mironov_openqp_2024]. However, there is no modern Fortran package that leverages the
use of fragmenation techniques to access massively parallel architectures.

Additionally, fragmentation has been historically underexplored due to the
software complexity of implementing fragmentation routines in existing
packages. This is mostly because fragmentation needs to be implemented
at the very top level of a software package in order to fully leverage
the parallelism offered by it. Currently, only four packages in the literature
support fragmentation natively[@zahariev_general_2023; @galvez_vallejo_toward_2023; @takami_open-architecture_2007; @broderick_span_2025] - the rest, offer fragmentatin capabilities
through file based interfaces, which can be difficult to scale beyond a
couple processes.

Metalquicha, offers a Fortran based fragmentation framework that is able
to bind to external quanutm chemistry engines, provided a C like interface
is available to compute the core physical properties required. This aims to
simplify the adoption of fragmentation routines by allowing existing quantum
chemistry programs to write interfaces such as `calculate_quantity` where
quantity can be energy, gradient, hessian, etc. This allows existing
programs to access high level parallelism and work-load distribution
frameworks.

Metalquicha addresses the following community needs:

- Fortran based programs have been scared into adopting C++ for modernisation, Metalquicha aims to showcase that Fortran can still be used for massively parallel  quanutm chemistry applications
- A permissive licensed, reusable software that can be extended freely
- A modern modular fragmentation framework that can accommodate any type of quanutm chemistry method

Additionally, Metalquicha uses modern Fortran language features such as native
documentation through the Ford project and the use of the [Fortran package manager](https://fpm.fortran-lang.org/) (FPM)
in addition to CMake to be compiled and consumed.

# Software Description

## Scope

Metalquicha contains two high level routines: MBE and GMBE[@richard_generalized_2012] serial and MPI enabled
fragment workload distribution schemes that support non-covalently bound, and
covalently bound fragments using either overlapping fragments or hydrogen caps.

At the core of the fragment distribution scheme is the `node_worker` subroutine which
assigns a specific task to a process. In most cases, this is the calculation of
a property, such as the energy, gradient, or hessian matrix of a molecule.

Currently, the `xtb` method is the only supported quantum chemistry "oracle". The
quantum chemistry engines are enabled via an abstract `calculation_method` object
that can be extended to support multiple types of methodologies. This allows
Metalquicha to be easily extended to replace `xtb` with a Restricted Hartree-Fock
(RHF) method, Density Functional Theory (DFT), among others.

## Current capabilities

Provided enough work, Metalquicha will scale up to an arbitrary number of processors, since the (G)MBE is naively parallel.

The main bottleneck of the simulation should be obtaining the property of the chemical
system however, for cheap methods such as `xtb`, the final computation of the fragmented
energy can become a bottleneck if enough fragments are created.



## Documentation

The code is thoroughly documented using FORD for program extendability and the use of
Metalquicha is covered in a read-the-docs hosted website.

## Testing

Automated unit tests and regression tests are run through Github Actions on every push,
ensuring the code is correct at every code modification.

## Input file format

In order to ensure portability, Metalquicha reads simple text files to process the chemical
input for a calculation. However, the front end uses a json format based on the QCSchema
and a python program to translate the json to the `mqc` format for Fortran consumption.


# Examples

## XTB engine

``` fortran
  select case (calc_type_local)
  case (CALC_TYPE_ENERGY)
     call xtb_calc%calc_energy(phys_frag, result)
  case (CALC_TYPE_GRADIENT)
     call xtb_calc%calc_gradient(phys_frag, result)
  case (CALC_TYPE_HESSIAN)
      call xtb_calc%calc_hessian(phys_frag, result)
  case default
     call logger%error("Unknown calc_type: "//calc_type_to_string(calc_type_local))
     error stop "Invalid calc_type in do_fragment_work"
  end select
```

where

``` fortran
  type, extends(qc_method_t) :: xtb_method_t
     !! Extended Tight-Binding (xTB) method implementation
     !!
     !! Concrete implementation of the abstract quantum chemistry method
     !! interface for GFN1-xTB and GFN2-xTB calculations via tblite.
     character(len=:), allocatable :: variant  !! XTB variant: "gfn1" or "gfn2"
     logical :: verbose = .false.              !! Print calculation details
     real(wp) :: accuracy = 0.01_wp            !! Numerical accuracy parameter
  contains
     procedure :: calc_energy => xtb_calc_energy      !! Energy-only calculation
     procedure :: calc_gradient => xtb_calc_gradient  !! Energy + gradient calculation
     procedure :: calc_hessian => xtb_calc_hessian    !! Finite difference Hessian
  end type xtb_method_t
```

which is based upon

``` fortran

   type, abstract :: qc_method_t
      !! Abstract base type for all quantum chemistry methods
      !!
      !! Defines the required interface for energy and gradient calculations
      !! that must be implemented by all concrete method types (XTB, HF, etc.).
   contains
      procedure(calc_energy_interface), deferred :: calc_energy    !! Energy calculation interface
      procedure(calc_gradient_interface), deferred :: calc_gradient  !! Gradient calculation interface
      procedure(calc_hessian_interface), deferred :: calc_hessian    !! Numerical hessian interface
   end type qc_method_t
```

Which can be arbitrarily extended for any type of method. Furthermore, once a
calculation method has been defined via the provided interface, the overall inclusion
into the MPI distribution scheme only relies on enabling the method via the configuration.

The MPI distribution is simply:

``` fortran
      ! Execute appropriate role
      if (world_comm%size() == 1) then
         ! Single rank: process fragments serially
         call logger%info("Running in serial mode (single MPI rank)")
         if (allow_overlapping_fragments) then
            ! GMBE serial processing with PIE coefficients
          call serial_gmbe_pie_processor(pie_atom_sets, pie_coefficients, n_pie_terms, sys_geom, method, calc_type, bonds)
         else
            ! Standard MBE serial processing
            call serial_fragment_processor(total_fragments, polymers, max_level, sys_geom, method, calc_type, bonds)
         end if
      else if (world_comm%leader() .and. node_comm%leader()) then
         ! Global coordinator (rank 0, node leader on node 0)
         call omp_set_num_threads(omp_get_max_threads())
         call logger%verbose("Rank 0: Acting as global coordinator")
         if (allow_overlapping_fragments) then
            ! GMBE MPI processing
            call gmbe_coordinator(world_comm, node_comm, n_monomers, polymers, intersections, intersection_sets, &
                   intersection_levels, n_intersections, node_leader_ranks, num_nodes, sys_geom, method, calc_type, bonds)
         else
            ! Standard MBE MPI processing
            call global_coordinator(world_comm, node_comm, total_fragments, polymers, max_level, &
                                    node_leader_ranks, num_nodes, sys_geom, calc_type)
         end if
      else if (node_comm%leader()) then
         ! Node coordinator (node leader on other nodes)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as node coordinator")
         call node_coordinator(world_comm, node_comm, calc_type)
      else
         ! Worker
         call omp_set_num_threads(1)
         call logger%verbose("Rank "//to_char(world_comm%rank())//": Acting as worker")
         call node_worker(world_comm, node_comm, sys_geom, method, calc_type, bonds)
      end if

```

Once a new method is enabled, unless one is adding a new fragmentation method, the MPI
workload distribution won't need any changes.

## Fragmented calculation

The input file format is:

```
{
    "schema": { "name": "mqc-frag", "version": "1.0" },
    "molecules": [{
        "xyz": "sample_inputs/prism.xyz",
        "fragments": [
          [0,1,2],
          [3,4,5],
          [6,7,8],
          [9,10,11],
          [12,13,14],
          [15,16,17]
        ],
        "fragment_charges": [0,0,0,0,0,0],
        "fragment_multiplicities": [1, 1, 1, 1, 1, 1],
        "molecular_charge": 0,
        "molecular_multiplicity": 1
    }],
  "model": {
    "method": "XTB-GFN1"
  },
  "keywords":{
    "scf": {
      "maxiter": 300,
      "tolerance": 1e-6
    },
    "fragmentation":{
        "method": "MBE",
        "allow_overlapping_fragments": false,
        "level": 2
    }
  },
  "system": {
    "logger":{
      "level": "Verbose"
    }
  },
  "driver": "Energy"
}
```

A python processor turns this into a keyword based file, for example:

```
%schema
name = mqc-frag
version = 1.0
index_base = 0
units = angstrom
end  ! schema

%model
method = XTB-GFN1
end  ! model

%driver
type = Energy
end  ! driver

%system
log_level = Verbose
end  ! system
```

Which is simply parsed by Metalquicha and processed into internal configurations.

The internal configuration is extendable, using derived types as containers for
new keywords, making adding new configuration parameters very simple.


# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

I would like to acknowledge the Fortran lang community for the development of projects
used in Metalquicha, such as the FPM, LFortran, and the inspiration found in the stdlib project.

I also would like to thank Ivan Pribec for listening to my programming debacles.

Finally, I'd like to thank my wife Dr. Elizabeth Gehrmann for enduring me developing this
code in my free time.

# References

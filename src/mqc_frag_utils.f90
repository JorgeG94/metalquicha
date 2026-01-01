!! Fragment generation and manipulation utilities
module mqc_frag_utils
   !! Provides combinatorial functions and algorithms for generating molecular
   !! fragments, managing fragment lists, and performing many-body expansion calculations.
   !!
   !! This module re-exports functionality from specialized modules:
   !! - mqc_combinatorics: Pure combinatorial mathematics
   !! - mqc_fragment_lookup: Hash-based fragment index lookup
   !! - mqc_gmbe_utils: GMBE intersection and PIE enumeration
   use mqc_combinatorics, only: &
      binomial, &
      get_nfrags, &
      create_monomer_list, &
      generate_fragment_list, &
      combine, &
      get_next_combination, &
      next_combination_init, &
      next_combination, &
      print_combos
   use mqc_fragment_lookup, only: fragment_lookup_t
   use mqc_gmbe_utils, only: &
      find_fragment_intersection, &
      generate_intersections, &
      compute_polymer_atoms, &
      generate_polymer_intersections, &
      gmbe_enumerate_pie_terms
   implicit none
   private

   ! Re-export from mqc_combinatorics
   public :: binomial
   public :: create_monomer_list
   public :: generate_fragment_list
   public :: combine
   public :: get_nfrags
   public :: get_next_combination
   public :: next_combination_init
   public :: next_combination
   public :: print_combos

   ! Re-export from mqc_fragment_lookup
   public :: fragment_lookup_t

   ! Re-export from mqc_gmbe_utils
   public :: find_fragment_intersection
   public :: generate_intersections
   public :: compute_polymer_atoms
   public :: generate_polymer_intersections
   public :: gmbe_enumerate_pie_terms

end module mqc_frag_utils

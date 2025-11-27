module mqc_elements
   use pic_ascii, only: to_upper, to_lower
   implicit none
   private

   public :: element_symbol_to_number, element_number_to_symbol

contains

   pure function element_symbol_to_number(symbol) result(atomic_number)
      !! Convert element symbol to atomic number
      !! Covers the complete periodic table (elements 1-118)
      character(len=*), intent(in) :: symbol
      integer :: atomic_number

      character(len=2) :: sym

      ! Normalize: uppercase first letter, lowercase second
      sym = adjustl(symbol)
      if (len_trim(sym) >= 1) sym(1:1) = to_upper(sym(1:1))
      if (len_trim(sym) >= 2) sym(2:2) = to_lower(sym(2:2))

      select case (trim(sym))
         ! Period 1
      case ('H'); atomic_number = 1
      case ('He'); atomic_number = 2
         ! Period 2
      case ('Li'); atomic_number = 3
      case ('Be'); atomic_number = 4
      case ('B'); atomic_number = 5
      case ('C'); atomic_number = 6
      case ('N'); atomic_number = 7
      case ('O'); atomic_number = 8
      case ('F'); atomic_number = 9
      case ('Ne'); atomic_number = 10
         ! Period 3
      case ('Na'); atomic_number = 11
      case ('Mg'); atomic_number = 12
      case ('Al'); atomic_number = 13
      case ('Si'); atomic_number = 14
      case ('P'); atomic_number = 15
      case ('S'); atomic_number = 16
      case ('Cl'); atomic_number = 17
      case ('Ar'); atomic_number = 18
         ! Period 4
      case ('K'); atomic_number = 19
      case ('Ca'); atomic_number = 20
      case ('Sc'); atomic_number = 21
      case ('Ti'); atomic_number = 22
      case ('V'); atomic_number = 23
      case ('Cr'); atomic_number = 24
      case ('Mn'); atomic_number = 25
      case ('Fe'); atomic_number = 26
      case ('Co'); atomic_number = 27
      case ('Ni'); atomic_number = 28
      case ('Cu'); atomic_number = 29
      case ('Zn'); atomic_number = 30
      case ('Ga'); atomic_number = 31
      case ('Ge'); atomic_number = 32
      case ('As'); atomic_number = 33
      case ('Se'); atomic_number = 34
      case ('Br'); atomic_number = 35
      case ('Kr'); atomic_number = 36
         ! Period 5
      case ('Rb'); atomic_number = 37
      case ('Sr'); atomic_number = 38
      case ('Y'); atomic_number = 39
      case ('Zr'); atomic_number = 40
      case ('Nb'); atomic_number = 41
      case ('Mo'); atomic_number = 42
      case ('Tc'); atomic_number = 43
      case ('Ru'); atomic_number = 44
      case ('Rh'); atomic_number = 45
      case ('Pd'); atomic_number = 46
      case ('Ag'); atomic_number = 47
      case ('Cd'); atomic_number = 48
      case ('In'); atomic_number = 49
      case ('Sn'); atomic_number = 50
      case ('Sb'); atomic_number = 51
      case ('Te'); atomic_number = 52
      case ('I'); atomic_number = 53
      case ('Xe'); atomic_number = 54
         ! Period 6
      case ('Cs'); atomic_number = 55
      case ('Ba'); atomic_number = 56
      case ('La'); atomic_number = 57
      case ('Ce'); atomic_number = 58
      case ('Pr'); atomic_number = 59
      case ('Nd'); atomic_number = 60
      case ('Pm'); atomic_number = 61
      case ('Sm'); atomic_number = 62
      case ('Eu'); atomic_number = 63
      case ('Gd'); atomic_number = 64
      case ('Tb'); atomic_number = 65
      case ('Dy'); atomic_number = 66
      case ('Ho'); atomic_number = 67
      case ('Er'); atomic_number = 68
      case ('Tm'); atomic_number = 69
      case ('Yb'); atomic_number = 70
      case ('Lu'); atomic_number = 71
      case ('Hf'); atomic_number = 72
      case ('Ta'); atomic_number = 73
      case ('W'); atomic_number = 74
      case ('Re'); atomic_number = 75
      case ('Os'); atomic_number = 76
      case ('Ir'); atomic_number = 77
      case ('Pt'); atomic_number = 78
      case ('Au'); atomic_number = 79
      case ('Hg'); atomic_number = 80
      case ('Tl'); atomic_number = 81
      case ('Pb'); atomic_number = 82
      case ('Bi'); atomic_number = 83
      case ('Po'); atomic_number = 84
      case ('At'); atomic_number = 85
      case ('Rn'); atomic_number = 86
         ! Period 7
      case ('Fr'); atomic_number = 87
      case ('Ra'); atomic_number = 88
      case ('Ac'); atomic_number = 89
      case ('Th'); atomic_number = 90
      case ('Pa'); atomic_number = 91
      case ('U'); atomic_number = 92
      case ('Np'); atomic_number = 93
      case ('Pu'); atomic_number = 94
      case ('Am'); atomic_number = 95
      case ('Cm'); atomic_number = 96
      case ('Bk'); atomic_number = 97
      case ('Cf'); atomic_number = 98
      case ('Es'); atomic_number = 99
      case ('Fm'); atomic_number = 100
      case ('Md'); atomic_number = 101
      case ('No'); atomic_number = 102
      case ('Lr'); atomic_number = 103
      case ('Rf'); atomic_number = 104
      case ('Db'); atomic_number = 105
      case ('Sg'); atomic_number = 106
      case ('Bh'); atomic_number = 107
      case ('Hs'); atomic_number = 108
      case ('Mt'); atomic_number = 109
      case ('Ds'); atomic_number = 110
      case ('Rg'); atomic_number = 111
      case ('Cn'); atomic_number = 112
      case ('Nh'); atomic_number = 113
      case ('Fl'); atomic_number = 114
      case ('Mc'); atomic_number = 115
      case ('Lv'); atomic_number = 116
      case ('Ts'); atomic_number = 117
      case ('Og'); atomic_number = 118
      case default
         atomic_number = 0  ! Unknown element
      end select

   end function element_symbol_to_number

   pure function element_number_to_symbol(atomic_number) result(symbol)
      !! Convert atomic number to element symbol
      !! Covers the complete periodic table (elements 1-118)
      integer, intent(in) :: atomic_number
      character(len=2) :: symbol

      select case (atomic_number)
         ! Period 1
      case (1); symbol = 'H '
      case (2); symbol = 'He'
         ! Period 2
      case (3); symbol = 'Li'
      case (4); symbol = 'Be'
      case (5); symbol = 'B '
      case (6); symbol = 'C '
      case (7); symbol = 'N '
      case (8); symbol = 'O '
      case (9); symbol = 'F '
      case (10); symbol = 'Ne'
         ! Period 3
      case (11); symbol = 'Na'
      case (12); symbol = 'Mg'
      case (13); symbol = 'Al'
      case (14); symbol = 'Si'
      case (15); symbol = 'P '
      case (16); symbol = 'S '
      case (17); symbol = 'Cl'
      case (18); symbol = 'Ar'
         ! Period 4
      case (19); symbol = 'K '
      case (20); symbol = 'Ca'
      case (21); symbol = 'Sc'
      case (22); symbol = 'Ti'
      case (23); symbol = 'V '
      case (24); symbol = 'Cr'
      case (25); symbol = 'Mn'
      case (26); symbol = 'Fe'
      case (27); symbol = 'Co'
      case (28); symbol = 'Ni'
      case (29); symbol = 'Cu'
      case (30); symbol = 'Zn'
      case (31); symbol = 'Ga'
      case (32); symbol = 'Ge'
      case (33); symbol = 'As'
      case (34); symbol = 'Se'
      case (35); symbol = 'Br'
      case (36); symbol = 'Kr'
         ! Period 5
      case (37); symbol = 'Rb'
      case (38); symbol = 'Sr'
      case (39); symbol = 'Y '
      case (40); symbol = 'Zr'
      case (41); symbol = 'Nb'
      case (42); symbol = 'Mo'
      case (43); symbol = 'Tc'
      case (44); symbol = 'Ru'
      case (45); symbol = 'Rh'
      case (46); symbol = 'Pd'
      case (47); symbol = 'Ag'
      case (48); symbol = 'Cd'
      case (49); symbol = 'In'
      case (50); symbol = 'Sn'
      case (51); symbol = 'Sb'
      case (52); symbol = 'Te'
      case (53); symbol = 'I '
      case (54); symbol = 'Xe'
         ! Period 6
      case (55); symbol = 'Cs'
      case (56); symbol = 'Ba'
      case (57); symbol = 'La'
      case (58); symbol = 'Ce'
      case (59); symbol = 'Pr'
      case (60); symbol = 'Nd'
      case (61); symbol = 'Pm'
      case (62); symbol = 'Sm'
      case (63); symbol = 'Eu'
      case (64); symbol = 'Gd'
      case (65); symbol = 'Tb'
      case (66); symbol = 'Dy'
      case (67); symbol = 'Ho'
      case (68); symbol = 'Er'
      case (69); symbol = 'Tm'
      case (70); symbol = 'Yb'
      case (71); symbol = 'Lu'
      case (72); symbol = 'Hf'
      case (73); symbol = 'Ta'
      case (74); symbol = 'W '
      case (75); symbol = 'Re'
      case (76); symbol = 'Os'
      case (77); symbol = 'Ir'
      case (78); symbol = 'Pt'
      case (79); symbol = 'Au'
      case (80); symbol = 'Hg'
      case (81); symbol = 'Tl'
      case (82); symbol = 'Pb'
      case (83); symbol = 'Bi'
      case (84); symbol = 'Po'
      case (85); symbol = 'At'
      case (86); symbol = 'Rn'
         ! Period 7
      case (87); symbol = 'Fr'
      case (88); symbol = 'Ra'
      case (89); symbol = 'Ac'
      case (90); symbol = 'Th'
      case (91); symbol = 'Pa'
      case (92); symbol = 'U '
      case (93); symbol = 'Np'
      case (94); symbol = 'Pu'
      case (95); symbol = 'Am'
      case (96); symbol = 'Cm'
      case (97); symbol = 'Bk'
      case (98); symbol = 'Cf'
      case (99); symbol = 'Es'
      case (100); symbol = 'Fm'
      case (101); symbol = 'Md'
      case (102); symbol = 'No'
      case (103); symbol = 'Lr'
      case (104); symbol = 'Rf'
      case (105); symbol = 'Db'
      case (106); symbol = 'Sg'
      case (107); symbol = 'Bh'
      case (108); symbol = 'Hs'
      case (109); symbol = 'Mt'
      case (110); symbol = 'Ds'
      case (111); symbol = 'Rg'
      case (112); symbol = 'Cn'
      case (113); symbol = 'Nh'
      case (114); symbol = 'Fl'
      case (115); symbol = 'Mc'
      case (116); symbol = 'Lv'
      case (117); symbol = 'Ts'
      case (118); symbol = 'Og'
      case default
         symbol = 'Xx'  ! Unknown
      end select

   end function element_number_to_symbol

end module mqc_elements

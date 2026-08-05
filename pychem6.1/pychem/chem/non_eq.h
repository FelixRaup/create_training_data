!!
!! Written by S. Glover, AMNH, 2004-2005, AIP, 2006-2007
!!

      integer nchem_network

!!#define DEBUG_EVOLVE

#if CHEMISTRYNETWORK == 4
      parameter(nchem_network = 4) !!CHEMISTRYNETWORK4
#endif

#if CHEMISTRYNETWORK == 5
      parameter (nchem_network = 5) !!CHEMISTRYNETWORK5
#endif
#if CHEMISTRYNETWORK == 6
      parameter (nchem_network = 6) !!CHEMISTRYNETWORK6
#endif

#if CHEMISTRYNETWORK == 15
      parameter (nchem_network = 15) !!CHEMISTRYNETWORK15
#endif
#if CHEMISTRYNETWORK == 17
      parameter (nchem_network = 17) !!CHEMISTRYNETWORK17
#endif

!! Set up quantities (such as the absolute tolerances) that are used in
!! multiple places in the non-equilibrium chemistry code. Note that most
!! DVODE-specific setup should go in evolve_abundances.F -- nrpar & nipar
!! are exceptions, as they are used elsewhere, so it is useful to define
!! them here
!!
      integer nrpar, nipar, nrrec
      parameter (nrpar=10)
#ifdef TR_ONTHESPOT
!! ** SV ** Add the heating rates and uv photon rates to rpar
      parameter (nrrec=TR_OS_NRECSPEC)
#else
      parameter (nrrec=0)
#endif
!! ** MW ** Added extra ipar for flags controlling chemistry behaviour
      parameter (nipar=7)
!! ** SV ** Added three extra ipar (for indices of recomb bands)
!!      parameter (nipar=6)
!! ** JM ** Added an extra ipar (Nbins for X-ray radiation)
!!      parameter (nipar=3)
!! ** JM **
!!      parameter (nipar=2)  !! non-XDR version

      integer num_non_eq_species
#if CHEMISTRYNETWORK == 4
      parameter (num_non_eq_species = 2)
#endif
#if CHEMISTRYNETWORK == 5
      parameter (num_non_eq_species = 3)
#endif
#if CHEMISTRYNETWORK == 6
      parameter (num_non_eq_species = 4)
#endif
#if CHEMISTRYNETWORK == 15
      parameter (num_non_eq_species = 9)
#endif
#if CHEMISTRYNETWORK == 17
      parameter (num_non_eq_species = 11)
#endif

      integer nspec
      parameter (nspec = num_non_eq_species+1)
!!      REAL non_eq_abundances(num_non_eq_species)


      REAL ATOL(nspec), rtol(nspec)
      common /tolerance/ ATOL, rtol
!!!!      common /abundances/ non_eq_abundances

!! Amount by which abundances are allowed to stray over their theoretical
!! maximum before triggering an error in rate_eq -- set to a blanket value
!! of 1d-4 for the time being...

      REAL eps_max
      parameter (eps_max = 1e-4)

!!      REAL atol(nspec)
#define  RTOL_FIX  1e-4
#define  ATOL_H2   1e-14
#define  ATOL_HP   1e-14
#define  ATOL_CP   1e-16
#define  ATOL_HEP  1e-14
#define  ATOL_CO   1e-14
#define  ATOL_HCOP 1e-18
#define  ATOL_CHX  1e-18
#define  ATOL_OHX  1e-14
#define  ATOL_MP   1e-14
#define  ATOL_C2P  1e-16
#define  ATOL_HE2P 1e-14
#define  ATOL_TMP  0e0
!!      common /tolerance/ atol

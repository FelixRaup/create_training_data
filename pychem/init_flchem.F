      subroutine init_flchem(abundC_ifc, abundO_ifc
     $                       ,dust_to_gas_ratio_ifc
     $                       ,cosmic_ray_ion_rate_ifc
     $                       ,NH_ext_ifc, Z_atom_ifc)
      implicit NONE


#include "Flash.h"
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"

      REAL abundC_ifc, abundO_ifc, dust_to_gas_ratio_ifc
      REAL G0_ifc, cosmic_ray_ion_rate_ifc
      REAL NH_ext_ifc, Z_atom_ifc



#ifdef debug_pychem
#if CHEMISTRYNETWORK == 5
      print*, "Inside if CHEMISTRYNETWORK==5"
#endif
#if CHEMISTRYNETWORK == 15
      print*, "Inside if CHEMISTRYNETWORK==15"
#endif
#if CHEMISTRYNETWORK == 17
      print*, "Inside if CHEMISTRYNETWORK==17"
#endif

#ifdef TR_ONTHESPOT
      print*, "OTS included!"
      print*, tr_osIHII, tr_osIHeII, tr_osIHeIII
#endif
#endif

      abundC = abundC_ifc
      abundO = abundO_ifc
      dust_to_gas_ratio = dust_to_gas_ratio_ifc
      cosmic_ray_ion_rate = cosmic_ray_ion_rate_ifc
      NH_ext = NH_ext_ifc
      Z_atom = Z_atom_ifc

c These are the same defaults as in the Config file, and are things
c we probably wont want to change while testing
      deff = 1e0
      abundSi = 1.5e-5
      abundD  = 2.6e-5
      abundM  = 1e-7
      abundN  = 0e0
      f_rsc   = 1e0
      phi_pah = 0.5e0
      AV_conversion_factor = 5.348e-22
      redshift = 0e0
      h2_form_ex  = 0e0
      h2_form_kin = 0e0
      xray_scaling = 1e0
      minDeltaT = 1e0

      iphoto   = 0
      iflag_mn = 1
      iflag_ad = 1
      iflag_atom = 3
      iflag_3bh2a = 1
      iflag_3bh2b = 1
      iflag_h3pra = 1
      isrf_option = 1
      no_chem     = 0
      no_cooling  = 0

      abundHe = abhe

      useFullNewton = .true.
      useAssistH2Decrease = .true.

c Set AV_ext based on specificed NH_ext, conversion factor and dust:gas ratio
      AV_ext = AV_conversion_factor * NH_ext * dust_to_gas_ratio

      call coolinmo
      call cheminmo
      call init_tolerances

      return
      end
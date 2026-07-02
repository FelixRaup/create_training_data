      subroutine test_flchem(total_time, timestep
     $                       ,G0_ifc, G0RP_ifc
     $                       ,fshield_H2, fshield_CO, AV_mean, chi_mean
     $                       ,dl, divv
     $                       ,density, internal_energy
     $                       ,abh2, abhp, abco, abcp
#ifdef TR_ONTHESPOT
     $                       ,pilc, pelc
#endif
     $                       ,temperature, dust_temp, tend
     $                       ,istatus
     $                       ,cool_rates_ifc, cool_rates_chem_ifc)
      implicit NONE

c No need to set correctly: temperature, tend, cool_rates, cool_rates_chem
c abcp will be ignored for CHEMISTRYNETWORK != 15
c WHAT HAPPENS TO THE VARS THAT ARE NOT LISTED BELOW,
c BUT IN SOME .h (lifted from test_harness)?
c Can I remove them from the F90->F77 interface and just include the .h in c_flchem?

c COMMON VARIABLES:
c REAL (cool.h) abundC, abundO, dust_to_gas_ratio
c REAL (cool.h) G0, cosmic_ray_ion_rate,
c REAL (cool.h) NH_ext, Z_atom,
c REAL (cool_dump.h->cool.h) cool_rates(nrates)
c REAL (cool_dump.h->cool.h) cool_rates_chem(nrates_chem)
c (cool_dump.h->cool.h): nrates=30, nrates_chem=8

#include "Flash.h"
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"

      REAL G0_ifc, G0RP_ifc
      REAL cool_rates_ifc(30)
      REAL cool_rates_chem_ifc(8)

      REAL total_time, timestep, dt, dust_temp
      REAL tstart, tend
      REAL fshield_H2, fshield_CO, AV_mean, chi_mean
      REAL dl, divv, density, internal_energy, yn
      REAL abh2, abhp, abco, abcp, abHI, abHeI, abe, abCI
      REAL ch_mh, ch_muHe, ch_muC, ch_muO, ch_muSi, ch_abar
      REAL non_eq_abundances(num_non_eq_species)
      REAL temperature, yntot, abhtot
      REAL G0RP, cr_rescale
      integer nstep_max
      integer istatus, turbostep, continuous_dgcouple
      integer I, nstep
      logical last_step
      REAL, parameter :: gamma=1.6667
#if TR_OS_NRECSPEC == 1
      integer, parameter :: tr_osIHII=0
      integer, parameter :: tr_osIHeII=-1
      integer, parameter :: tr_osIHeIII=-1
#elif TR_OS_NRECSPEC == 2
      integer, parameter :: tr_osIHII=0
      integer, parameter :: tr_osIHeII=1
      integer, parameter :: tr_osIHeIII=-1
#elif TR_OS_NRECSPEC == 3
      integer, parameter :: tr_osIHII=0
      integer, parameter :: tr_osIHeII=1
      integer, parameter :: tr_osIHeIII=2
#endif
      real, dimension(nrrec) :: pilc, pelc

#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
! NOTE: We take the same value for ionized species, for molecules just the simple sum without floating points
      real, parameter :: ch_muHI = 1.008
      real, parameter :: ch_muH2 = 2.01588
      real, parameter :: ch_muHp = 1.0072
      real, parameter :: ch_muHep = 4.002602
      real, parameter :: ch_muCHx = 13.
      real, parameter :: ch_muOHx = 17.
      real, parameter :: ch_muCp = 12.011
      real, parameter :: ch_muHCOp = 29.
      real, parameter :: ch_muCO = 28.01
      real, parameter :: ch_muMp = 28. !We adopt silicium here
      real, parameter :: ch_muM = 28.
#endif

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

      nstep_max = 5000000
      G0RP = G0RP_ifc
      cr_rescale= 1.d0
      istatus = 4711
      turbostep = 2
      continuous_dgcouple = 1
      ch_mH   = 1.672621637e-24 !!1.66053e-24
      ch_muHe = 4.002602
      ch_muC  = 12.011
      ch_muO  = 15.9994
      ch_muSi = 28.0855

      dust_temp = 10.0
      tdust=dust_temp

      G0 = G0_ifc
      cool_rates = cool_rates_ifc
      cool_rates_chem = cool_rates_chem_ifc

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

c Mean atomic weight
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 6
      ch_abar = 1.0 + abhe*ch_muHe + abundC*ch_muC
     $        + abundO*ch_muO + abundSi*ch_muSi ! + abundM*ch_muM
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      ch_abar = 1.0*ch_muHI + abhe*ch_muHe + abundC*ch_muC
     $        + abundO*ch_muO + 0.*abundSi*ch_muSi + abundM*ch_muM
!   Mass fraction of all H atoms
      abhtot = ch_muHI/ch_abar
#endif

      dt = timestep
      yn = density / ch_mH /ch_abar
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      yn = density * abhtot /ch_mH
#endif
      non_eq_abundances(ih2) = abh2
      non_eq_abundances(ihp) = abhp
#if CHEMISTRYNETWORK == 5
      non_eq_abundances(ico) = abco
      abe = abhp
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      non_eq_abundances(icp) = abcp
      non_eq_abundances(ichx) = 1e-40
      non_eq_abundances(iohx) = 1e-40
      non_eq_abundances(ico) = abco
      non_eq_abundances(ihcop) = 1e-40
      non_eq_abundances(ihep) = 1e-40
      non_eq_abundances(imp) = 1e-40
      abe = abhp + abcp
     $    + non_eq_abundances(ihep)
     $    + non_eq_abundances(ihcop)
     $    + non_eq_abundances(imp)
#endif
#if CHEMISTRYNETWORK == 17
      non_eq_abundances(ic2p) = 1e-40
      non_eq_abundances(ihe2p) = 1e-40
      abe = abe
     $      + 2 * non_eq_abundances(ic2p)
     $      + 2 * non_eq_abundances(ihe2p)
#endif

      yntot = (1e0 + abhe - abh2 + abe) * yn
      temperature=internal_energy/yntot/boltz*2./3.
      nstep = int(total_time / timestep) !!+ 1
      tstart = 0e0

c      call coolinmo
c      call cheminmo
c      call init_tolerances

      do I = 1, nstep
        dt   = timestep
        tend = tstart + dt
        if (tend .gt. total_time) then
          dt = total_time - tstart
          last_step = .true.
        else
          last_step = .false.
        endif

        call evolve_abundances(dt, dl, yn, divv, internal_energy, 
     $                         redshift, non_eq_abundances,
     $                         fshield_H2, fshield_CO,
     $                         AV_mean, chi_mean,
     $                         dust_temp,
     $                         G0RP, cr_rescale,
#ifdef TR_ONTHESPOT
     $                         pilc, pelc,
     $                         tr_osIHII, tr_osIHeII, tr_osIHeIII,
#endif
     $                         nstep_max, istatus, turbostep,
     $                         continuous_dgcouple)
#if CHEMISTRYNETWORK == 5
        abe = non_eq_abundances(ihp)
#endif
#if (CHEMISTRYNETWORK == 15) || (CHEMISTRYNETWORK == 17)
        abe = non_eq_abundances(ihp)
     $      + non_eq_abundances(icp)
     $      + non_eq_abundances(ihcop)
     $      + non_eq_abundances(ihep)
     $      + non_eq_abundances(imp)
#endif
#if CHEMISTRYNETWORK == 17
        abe = abe
     $      + 2 * non_eq_abundances(ic2p)
     $      + 2 * non_eq_abundances(ihe2p)
#endif
        yntot = (1d0 + abhe - non_eq_abundances(ih2) + abe) * yn
        temperature=internal_energy/yntot/boltz*2./3.
        if (I .eq. 1 .or. I .eq. nstep) then
#if CHEMISTRYNETWORK == 5
          abHI  = max(1-2*non_eq_abundances(ih2)-non_eq_abundances(ihp)
     $           , 1e-40)
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
          abHI  = max(1-2*non_eq_abundances(ih2)-non_eq_abundances(ihp)
     $           - non_eq_abundances(ichx) - non_eq_abundances(iohx)
     $           - non_eq_abundances(ihcop), 1e-40)
          abHeI = abhe - non_eq_abundances(ihep)
          abCI   = abundC - non_eq_abundances(ico)
     $           - non_eq_abundances(icp)
     $           - non_eq_abundances(ichx) - non_eq_abundances(ihcop)
#endif
#if CHEMISTRYNETWORK == 17
          abCI  = abCI - non_eq_abundances(ic2p)
#endif
        endif

c UPDATE INTERFACE ARGUMENTS (FOR READBACK)
        cool_rates_ifc = cool_rates
        cool_rates_chem_ifc = cool_rates_chem

c UPDATE ABUNDANCES (FOR READBACK)
        abh2 = non_eq_abundances(ih2)
        abhp = non_eq_abundances(ihp)
        abco = non_eq_abundances(ico)
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
        abcp = non_eq_abundances(icp)
#endif
  
        if (last_step) then
          goto 100
        endif

        tstart = tend
      enddo

 100  continue
 
      return
      end
      program test
      implicit none
#include "Flash.h"
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"
      REAL total_time, timestep, dt, dust_temp
      REAL tstart, tend
      REAL fshield_H2, fshield_CO, AV_mean, chi_mean
      REAL dl, divv, density, internal_energy, yn, cs_amb
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
      real, dimension(nrrec) :: Lyc_uv_ion, Lyc_HeatingRate

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
      nstep_max = 5000000
      istatus = 4711
      turbostep = 2
      G0RP = 0
      cr_rescale = 1
      ch_mH   = 1.672621637e-24 !!1.66053e-24
      ch_muHe = 4.002602
      ch_muC  = 12.011
      ch_muO  = 15.9994
      ch_muSi = 28.0855

c Read in adjustable parameters, starting density, temp. and abundances

      open(1, file='test.inp', status='old')
      read(1,*) total_time, timestep
      read(1,*) abundC, abundO
      read(1,*) dust_to_gas_ratio, dust_temp
      read(1,*) G0, cosmic_ray_ion_rate
      read(1,*) NH_ext, Z_atom
      read(1,*) fshield_H2, fshield_CO, AV_mean, chi_mean
      read(1,*) dl, divv
      read(1,*) density, cs_amb
#if CHEMISTRYNETWORK == 15
      read(1,*) abh2, abhp, abco, abcp
#else
      read(1,*) abh2, abhp, abco
#endif
      read(1,*) Lyc_uv_ion
      read(1,*) Lyc_HeatingRate
      close(1, status='keep')

c Print the values we read in
      print*,'----------------------'
      print*, 'Parameters in test.inp'
      print*, total_time, timestep
      print*, abundC, abundO
      print*, dust_to_gas_ratio, dust_temp
      print*, G0, cosmic_ray_ion_rate
      print*, NH_ext, Z_atom
      print*, fshield_H2, fshield_CO, AV_mean, chi_mean
      print*, dl, divv
      print*, density, cs_amb
#if CHEMISTRYNETWORK == 15
      print*, abh2, abhp, abco, abcp
#else
      print*, abh2, abhp, abco
#endif
      print*, 'PILC =', Lyc_uv_ion
      print*, 'PELC =', Lyc_HeatingRate
      print*,'----------------------'

      G0RP = 0.d0
      cr_rescale=1
      istatus = 4711
      turbostep = 2
      continuous_dgcouple = 0

      tdust=dust_temp
c      internal_energy=cs_amb**2e0/(gamma*(gamma-1.))
      internal_energy=cs_amb**2e0 !!*3.e0/2.d0
      internal_energy=internal_energy*density
      internal_energy=internal_energy/(gamma*(gamma-1.))
c      print*,'internal energy=',internal_energy

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

      call coolinmo
      call cheminmo
      call init_tolerances

c Open logfile
      open(2, file='test.log', status='replace')

c Write initial setup to logfile
#if CHEMISTRYNETWORK == 5
        write(2,*) tstart, non_eq_abundances(ih2), 
     $             non_eq_abundances(ihp), 
     $             non_eq_abundances(ico),
     $             internal_energy, temperature
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
        write(2,*) tstart, non_eq_abundances(ih2), 
     $             non_eq_abundances(ihp), 
     $             non_eq_abundances(ico),
     $             non_eq_abundances(icp),
     $             non_eq_abundances(ichx),
     $             non_eq_abundances(iohx),
     $             non_eq_abundances(ihcop),
     $             non_eq_abundances(ihep),
     $             non_eq_abundances(imp),
     $             internal_energy, temperature
#endif
c #else
c         write(2,*) tstart, non_eq_abundances(ih2), 
c      $             non_eq_abundances(ihp), 0.0,
c      $             internal_energy,temperature
c #endif


      do I = 1, nstep
        dt   = timestep
        tend = tstart + dt
        if (tend .gt. total_time) then
          dt = total_time - tstart
          last_step = .true.
        else
          last_step = .false.
        endif

        print*,'before evolve_abundances'
c        print*,'abun (H2, Hp)=',non_eq_abundances
        print*,'eint=',internal_energy,' temp=',temperature
c        print*,dt, dl, yn, divv, internal_energy,
c     $                         redshift, non_eq_abundances,
c     $                         fshield_H2, fshield_CO, AV_mean,
c     $                         chi_mean, dust_temp



        call evolve_abundances(dt, dl, yn, divv, internal_energy, 
     $                         redshift, non_eq_abundances,
     $                         fshield_H2, fshield_CO,
     $                         AV_mean, chi_mean,
     $                         dust_temp,
     $                         G0RP, cr_rescale,
#ifdef TR_ONTHESPOT
     $                         Lyc_uv_ion, Lyc_HeatingRate,
     $                         tr_osIHII, tr_osIHeII, tr_osIHeIII,
#endif
#ifdef XRAY_TR
!! ** JM ** Add attenuated X-ray flux, Number of energy bins, and total energy absorbed
!!     $                         Fin, XR_nbins, HX,
!! ** JM **
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
          print*,'Step=',I
          print*,'eint=',internal_energy,' temp=',temperature
          print*, dt, dl, yn, divv, internal_energy,
     $                         redshift, non_eq_abundances,
     $                         fshield_H2, fshield_CO, AV_mean,
     $                         chi_mean, dust_temp
          print*,'abun H2 =', non_eq_abundances(ih2)
          print*,'abun H+ =', non_eq_abundances(ihp)
#if CHEMISTRYNETWORK == 5
          print*,'abun CO =', non_eq_abundances(ico)
          abHI  = max(1-2*non_eq_abundances(ih2)-non_eq_abundances(ihp)
     $           , 1e-40)
          print*,'abun Hatom =', abHI
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
          print*,'abun CO =', non_eq_abundances(ico)
          print*,'abun C+ =', non_eq_abundances(icp)
          print*,'abun CHx =', non_eq_abundances(ichx)
          print*,'abun OHx =', non_eq_abundances(iohx)
          print*,'abun HCO+ =', non_eq_abundances(ihcop)
          print*,'abun He+ =', non_eq_abundances(ihep)
          print*,'abun Si+ =', non_eq_abundances(imp)
          abHI  = max(1-2*non_eq_abundances(ih2)-non_eq_abundances(ihp)
     $           - non_eq_abundances(ichx) - non_eq_abundances(iohx)
     $           - non_eq_abundances(ihcop), 1e-40)
          abHeI = abhe - non_eq_abundances(ihep)
          abCI   = abundC - non_eq_abundances(ico)
     $           - non_eq_abundances(icp)
     $           - non_eq_abundances(ichx) - non_eq_abundances(ihcop)
          print*,'abun Hatom =', abHI
          print*,'abun He atom =', abHeI
#endif
#if CHEMISTRYNETWORK == 17
          print*,'abun C2+ =', non_eq_abundances(ic2p)
          print*,'abun He2+ =', non_eq_abundances(ihe2p)
          abCI  = abCI - non_eq_abundances(ic2p)
#endif
        endif


c Write output to logfile
#if CHEMISTRYNETWORK == 5
        write(2,*) tend, non_eq_abundances(ih2),
     $             non_eq_abundances(ihp),
     $             non_eq_abundances(ico),
     $             abundC - non_eq_abundances(ico),
     $             abHI,
     $             internal_energy, temperature
#ifdef THERMAL_DUMP_ALL_RATES
     $             ,cool_rates
     $             ,cool_rates_chem
#endif
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
        write(2,*) tend, non_eq_abundances(ih2),
     $             non_eq_abundances(ihp),
     $             non_eq_abundances(ico),
     $             non_eq_abundances(icp),
     $             abHI, abCI,
c     $             non_eq_abundances(ichx),
c     $             non_eq_abundances(iohx),
c     $             non_eq_abundances(ihcop),
c     $             non_eq_abundances(ihep),
c     $             non_eq_abundances(imp),
     $             internal_energy, temperature
#ifdef THERMAL_DUMP_ALL_RATES
     $             ,cool_rates
     $             ,cool_rates_chem
#endif
#endif

        if (last_step) then
          goto 100
        endif

        tstart = tend
      enddo

 100  continue
 
      end

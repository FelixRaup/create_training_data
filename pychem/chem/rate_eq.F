! #if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
! #ifndef KOYAMA
#include "Flash.h"
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             R A T E _ E Q                 \\\\\\\\\\
c
c=======================================================================
c
c - J.Mackey: 2016.05.23: Added X-ray reactions and heating rates.
c - 2016.05.23 JM+DS: Added dissociative recombination of H2+ to reaction pathways.
c              MS05 = Meijerink & Spaans (2005,A&A,436,397)
c - JM: 2016.12.07 Added cutoff to heating/cooling for T>1e5K.
c                  Added collisional ionization of C + e --> C+ + 2e (ch43)
c
      subroutine rate_eq(neq, time, y, ydot, rpar, ipar)
! #ifdef CHEMCOOL
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"
#include "Flash.h"
c #include "xray.h"
c #include "current_photo.h"
      integer neq
c
      REAL time, y(nspec), ydot(nspec)
      integer ipar(nipar)
      integer iflag_contcool
#ifdef XRAY_TR
c ** JM **
c The last ipar(3) elements of rpar are the attenuated x-ray fluxes in each energy bin.
c xray_ye:    electron fraction as defined in DYL99
c xray_rH2:   ratio by number of molecular to atomic hydrogen.
c xab():      abundances of species for calculation of x-ray ionization/dissociation rates
c xion():     array to hold x-ray ionization/dissociation rates.
c xray_Fin(): array to hold x-ray flux in each energy bin.
c GammaX:     X-ray Coulomb heating rate (erg/cm3/s)
c HX:         X-ray absorbed power (erg/s/H-nucleus)
c      REAL, DIMENSION(nrpar+xdr_nbins) :: rpar
      REAL, DIMENSION(nrpar+2*nrrec+ipar(3)) :: rpar
      REAL :: xray_ye, xray_rH2, xab(4), xion(5)
c      real, dimension(xdr_nbins) :: xray_Fin
      real, dimension(ipar(3)) :: xray_Fin
c      REAL, ALLOCATABLE :: xray_Fin(:)
      REAL :: GammaX, HX
c     xab: [H0, H_2, He0, He+]
c     rates xion [H->H+, H2->H2+, H2->2H, He->He+, He+->He2+]
      INTEGER irh, irh2, drh2, irhe, irhep
      parameter (irh   = 1)
      parameter (irh2  = 2)
      parameter (drh2  = 3)
      parameter (irhe  = 4)
      parameter (irhep = 5)
      REAL kcxh2   ! H_2^+ charge exchange with H, to produce H_2 and H^+
c ** JM **
#else
      REAL rpar(nrpar+2*nrrec)  ! non-XDR code
#endif

c
      integer itemp, i, in
      REAL abundances(nabn), rates(nrates),
     $     rates_chem(nrates_chem),
     $     y_in(nspec)
c
      REAL abe     , abhp    , abHI    , abh2    ,
     $     abcp    , abco    , abheI   , abh3p   ,
     $     abcI    , aboI    , abchx   , abohx   ,
     $     abmp    , abMI    , abhep   , abhcop  ,
     $     abop    , crco    , rohx    , cpSi    ,
     $     abotot  , abo     , abc2p   , abhe2p
c
      REAL ylam    , temp    ,
     $     yn      , dl      , divv    , yntot   , gamma ,
     $     energy  , dtemp   , ncrinv  , h2var0  , h2var1,
     $     h2var2  , ykdh    , ykdh2   , ykdis   , NH    ,
     $     ylam_chem, f_ex   , rho     , cs2     , L_jeans,
     $     f_tot   , RH2     , Rdust   , ekn     , AV_mean,
     $     G_dust
c
      REAL fshield_H2, fshield_CO, fshield_dust, chi_mean,
     $     fshield_dust_H2, gas_temp
c
#ifdef TR_ONTHESPOT
      REAL Lyc_uv_ion(nrrec), Lyc_HeatingRate(nrrec)
      INTEGER iHII, iHeII, iHeIII
#endif
c
c
      REAL ph2, ph10, ph47, ph48, ph49, ph50, ph51
c
      REAL fh3p
      REAL fh2, f2h
c
      REAL num, denom
c
      REAL k0, k1, gamma_chx, gamma_co, beta
      REAL f_kc1, f_ck2
c
      REAL fshield
c
      REAL phi     , h_gr    , hgrvar1 , hgrvar2, dii_gr
c
      REAL cr1     , cr2     , cr3     , cr6    , cr10  ,
     $     cr11    , cp1     , cp11    , cr13
c
      REAL ch1     , ch2     , ch3     , ch4     , ch5  ,
     $     ch6     , ch9     , ch11    , ch13    , ch17 ,
     $     ch24    , ch25    , ch34    , ch35    , ch40 ,
     $     ch56    , ch92    , ch105   , ch106   , ch146,
     $     ch147   , ch153   , ch154   , ch155   , ch156,
     $     ch157   , ch158   , ch159   , ch160   , ch161,
     $     ch162   , ch20
c
      REAL uv_ion  , uv_heat, photon_flux, uv_ion_cp,
     $     uv_ion_cI, uv_ion_cI_directly_to_c2p,
     $     Lyc_HeatingRate_HII, Lyc_HeatingRate_HeII,
     $     Lyc_HeatingRate_HeIII, uv_ion_HeII, uv_ion_HeIII
c

c ** JM ** : added
c ch30: Dissociative Recomb. of H2+ + e- --> 2H
c ch43: Collisional ionization of C
c ch44: Collisional ionization of Si (M in NL99)
c ch14: Collisional ionization of He
c ch36/37: grain recombination of C+
c ch38/39: grain recombination of Si+ (M+ in NL99)
c ch135/136: grain recombination of He+
c ** SCOG ** : added
c ch24: Production of H2+ by charge transfer of H+ to H2
      REAL ch30, ch43, ch44, ch14
      REAL ch36, ch37, ch38, ch39, ch135, ch136
      REAL ch46, ch47
      REAL khe_dust, kcp_dust, kh_dust, kmp_dust
c
      REAL dtch1   , dtch2   , dtch3   , dtch4   , dtch5  ,
     $     dtch6   , dtch9   , dtch11  , dtch13  , dtch17 ,
     $     dtch24  , dtch25  , dtch34  , dtch35  , dtch40 ,
     $     dtch56  , dtch92  , dtch105 , dtch106 , dtch146,
     $     dtch147 , dtch153 , dtch154 , dtch155 , dtch156,
     $     dtch157 , dtch158 , dtch159 , dtch160 , dtch161,
     $     dtch162 , dtch20
      REAL dtch30, dtch43, dtch44, dtch14
      REAL dtch46, dtch47
      REAL dtch36, dtch37, dtch38, dtch39, dtch135, dtch136
c
      REAL ch106a, ch155a
c
      REAL cdex, rdex, fpheat
c
      REAL pi
      parameter (pi = 3.14159265358979323846)
      integer log_HII_region
      REAL G0RP
      REAL cr_rescale


      log_HII_region = 0
c Unpack parameters
c
      yn   = rpar(1)
      dl   = rpar(2)
      divv = rpar(3)
      fshield_H2   = rpar(4)
      fshield_CO   = rpar(5)
      fshield_dust = rpar(6)
      chi_mean     = rpar(7)
      fshield_dust_H2 = rpar(8)
      G0RP = rpar(9)
      cr_rescale = rpar(10)
c      Lyc_uv_ion_H = rpar(9)
c      Lyc_uv_ion_He = rpar(10)
c      Lyc_HeatingRate = rpar(13)
c      Lyc_HeatingRate_He = rpar(14)
c
#ifdef TR_ONTHESPOT
c   get the EUV fluxes and energy fluxes from rpar
      do I = 1, nrrec
        Lyc_uv_ion(I) = rpar(nrpar+I)
        Lyc_HeatingRate(I) = rpar(nrpar+nrrec+I)
      enddo
      iHII = ipar(4) +1
      iHeII = ipar(5) +1
      iHeIII = ipar(6) +1
#endif
c Compute effective column density (needed for diffuse xray heating/ion)
c
c      if (dust_to_gas_ratio .gt. 0.0) then
c        if (fshield_dust .le. 0.0) then
c Pick some large value -- exact choice won't be important, as diffuse X-ray heating/ion not
c important in v. high column density gas
c          NH = 1e23
c        else
c          AV_mean = -log(fshield_dust) / 2.5e0
c          NH = AV_mean / (dust_to_gas_ratio * AV_conversion_factor)
c     $       + NH_ext
c        endif
c      else
c In non-primordial calculations, only place we might not have dust is
c very hot, low density gas, where X-ray heating & ion are negligible
c        NH = 0.0
c      endif
c
c Chemical species:
c
      do I = 1, nspec
        y_in(I) = y(I)
      enddo
c
#ifdef XRAY_TR
c  get x-ray fluxes from rpar.
      !!ALLOCATE ( xray_Fin (ipar(3)) )
      do I = 1, ipar(3)
      !!do I = 1, xdr_nbins
        xray_Fin(I) = rpar(nrpar+2*nrrec+I)
      enddo
#endif
c
      call validate_input(y, ipar)
      if (ipar(1) .eq. 1) then
c Return with zeroed ydot array if input abundances invalid
        do I = 1, nspec
          ydot(I) = 0e0
        enddo
        return
      endif
c
c Tracked, non-equilibrium species
c

#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      abhp = max(0e0, y(ihp))
      abh2 = max(0e0, y(ih2))
      abco  = 0e0
#if CHEMISTRYNETWORK == 5
      abco = max(0e0, y(ico))
#endif
#if CHEMISTRYNETWORK == 6
      abco = max(0e0, y(ico))
      abc2p = max(0e0, y(ic2p))
#endif
c
c Non-equilibrium species computed using conservation laws
c
      abHI = max(1e0 - 2e0 * abh2 - abhp, 0e0)
#if CHEMISTRYNETWORK == 5
      abcp = max(0e0, abundc - abco)
      abe  = abhp + abcp + abundsi
#elif CHEMISTRYNETWORK == 6
      abcp = max(0e0, abundc - abco - abc2p)
      abe  = abhp + abcp + abundsi + 2e0 * abc2p
#else
      abcp = abundc
      abe  = abhp + abundc + abundsi
#endif
#endif
c
#if (CHEMISTRYNETWORK == 15) || (CHEMISTRYNETWORK == 17)
      abh2   = max(1e-40, y(ih2))
      abhp   = max(1e-40, y(ihp))
      abcp   = max(1e-40, y(icp))
      abchx  = max(1e-40, y(ichx))
      abohx  = max(1e-40, y(iohx))
      abco   = max(1e-40, y(ico))
      abhcop = max(1e-40, y(ihcop))
      abhep  = max(1e-40, y(ihep))
      abmp   = max(1e-40, y(imp))
c
c Non-equilibrium species computed using conservation laws
c
      abe   = max(abhp + abcp + abhcop + abhep + abmp, 1e-40)
      abHI  = max(1d0 - 2. * abh2 - abhp - abchx - abohx - abhcop,
     $            1e-40)
      abheI = max(abhe - abhep, 1e-40)
      abcI  = max(abundc - abcp - abchx - abco - abhcop, 1e-40)
      abmI  = max(abundM - abmp, 1e-40)
c
c For O, we first calculate total of O, O+ abundances; we calculate
c fraction in O+ below, assuming this is always in equilibrium
c
      abotot = max(abundo - abohx - abco - abhcop, 1e-40)
#endif
#if CHEMISTRYNETWORK == 17
      abc2p  = max(1e-40, y(ic2p))
      abhe2p  = max(1e-40, y(ihe2p))
c
c Non-equilibrium species computed using conservation laws
c
      abe   = max(abhp + abcp + abhcop + abhep + abmp 
     $      + 2. * abc2p + 2. * abhe2p, 1e-40)
      abheI = max(abhe - abhep - abhe2p, 1e-40)
      abcI  = max(abundc - abcp - abchx - abco - abhcop - abc2p, 1e-40)
#endif
c
c Internal energy (erg cm^-3)
c
      energy = y(itmp)
c
c Compute temperature
c
c [NB We ignore corrections to yntot & gamma arising from low abundance
c     molecular species (e.g. H-, H2+) and metals]
c
      yntot = (1d0 + abhe - abh2 + abe) * yn
      ekn   = energy / (yntot * kboltz)
c For time being, assume gamma = 5/3 -- good approx. unless we have lots of hot H2
      gamma = 1.6667
      temp = (2e0 / 3e0) * ekn

c      print *, temp, abh2, abhp, abcp, abchx, abohx, abco,
c     $ abcp, abchx, abco, abhcop, abhep, abmp, abe, abHI,
c     $ abheI, abcI, abmI

c
c     call lookup_temp_gamma(abh2, ekn, temp, gamma)
c
c Bounds checking
c
c Temp:
c
      if (temp .le. 0d0) then
        itemp   = 1
        dtemp   = 0d0
        temp    = tmin
        ipar(1) = 1
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = 0d0
      else
        itemp = int(log10(temp) / dtlog) + 1
        if (itemp .le. 0 .or. itemp .gt. nmd) then
          print*, 'Fatal error in rate_eq.F', itemp, temp, energy
          ABORT(1)
        endif
        dtemp = temp - temptab(itemp)
      endif
c
c Read in coefficients
c
      dtch1  = dtchtab(1, itemp)
      dtch2  = dtchtab(2, itemp)
      dtch3  = dtchtab(3, itemp)
      dtch4  = dtchtab(4, itemp)
      dtch5  = dtchtab(5, itemp)
      dtch6  = dtchtab(6, itemp)
      dtch11 = dtchtab(11, itemp)
      dtch13 = dtchtab(13, itemp)
      dtch25 = dtchtab(25, itemp)
      dtch34 = dtchtab(34, itemp)
      dtch35 = dtchtab(35, itemp)
      dtch160 = dtchtab(160, itemp) ! C+ photo ion.
      dtch161 = dtchtab(161, itemp) ! recomb electron energy loss factor for HI
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      dtch9  = dtchtab(9, itemp)
      dtch14 = dtchtab(14, itemp)  !! He coll.ion.
      dtch17 = dtchtab(17, itemp)
      dtch20 = dtchtab(20, itemp)  !! HeIII case B
      dtch24 = dtchtab(24, itemp)
      dtch30 = dtchtab(30, itemp)
      dtch36 = dtchtab(36, itemp)
      dtch37 = dtchtab(37, itemp)
      dtch38 = dtchtab(38, itemp)
      dtch39 = dtchtab(39, itemp)
      dtch40 = dtchtab(40, itemp)
      dtch43 = dtchtab(43, itemp)  !! C coll.ion.
      dtch44 = dtchtab(44, itemp)  !! M coll.ion.
      dtch46 = dtchtab(46, itemp)
      dtch47 = dtchtab(47, itemp)
      dtch56 = dtchtab(56, itemp)
      dtch92 = dtchtab(92, itemp)
      dtch105 = dtchtab(105, itemp)
      dtch106 = dtchtab(106, itemp)
      dtch135 = dtchtab(135, itemp)
      dtch136 = dtchtab(136, itemp)
      dtch146 = dtchtab(146, itemp)
      dtch147 = dtchtab(147, itemp)
      dtch153 = dtchtab(153, itemp)
      dtch154 = dtchtab(154, itemp)
      dtch155 = dtchtab(155, itemp)
      dtch156 = dtchtab(156, itemp)
      dtch157 = dtchtab(157, itemp)
      dtch158 = dtchtab(158, itemp)
      dtch159 = dtchtab(159, itemp)
      dtch162 = dtchtab(162, itemp) ! recomb electron energy loss factor for HeI
#endif
c
      ch1  = chtab(1, itemp) + dtemp * dtch1
      ch2  = chtab(2, itemp) + dtemp * dtch2
      ch3  = chtab(3, itemp) + dtemp * dtch3
      ch4  = chtab(4, itemp) + dtemp * dtch4
      ch5  = chtab(5, itemp) + dtemp * dtch5
      ch6  = chtab(6, itemp) + dtemp * dtch6
      ch11 = chtab(11, itemp) + dtemp * dtch11
      ch13 = chtab(13, itemp) + dtemp * dtch13
      ch25 = chtab(25, itemp) + dtemp * dtch25
      ch34 = chtab(34, itemp) + dtemp * dtch34
      ch35 = chtab(35, itemp) + dtemp * dtch35
      ch160 = chtab(160, itemp) + dtemp * dtch160 ! C+ photo ion.
      ch161 = chtab(161, itemp) + dtemp * dtch161 ! recomb electron energy loss factor for HI
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      ch9  = chtab(9, itemp) + dtemp * dtch9
      ch14 = chtab(14, itemp) + dtemp * dtch14  ! He
      ch17 = chtab(17, itemp) + dtemp * dtch17
      ch20 = chtab(20, itemp) + dtemp * dtch20  ! HeIII case B
      ch24 = chtab(24, itemp) + dtemp * dtch24
      ch30 = chtab(30, itemp) + dtemp * dtch30
      ch36 = chtab(36, itemp) + dtemp * dtch36
      ch37 = chtab(37, itemp) + dtemp * dtch37
      ch38 = chtab(38, itemp) + dtemp * dtch38
      ch39 = chtab(39, itemp) + dtemp * dtch39
      ch40 = chtab(40, itemp) + dtemp * dtch40
      ch43 = chtab(43, itemp) + dtemp * dtch43  ! C
      ch44 = chtab(44, itemp) + dtemp * dtch44  ! Si (M)
      ch46 = chtab(46, itemp) + dtemp * dtch46
      ch47 = chtab(47, itemp) + dtemp * dtch47
      ch56 = chtab(56, itemp) + dtemp * dtch56
      ch92 = chtab(92, itemp) + dtemp * dtch92
      ch105 = chtab(105, itemp) + dtemp * dtch105
      ch106 = chtab(106, itemp) + dtemp * dtch106
      ch135 = chtab(135, itemp) + dtemp * dtch135
      ch136 = chtab(136, itemp) + dtemp * dtch136
      ch146 = chtab(146, itemp) + dtemp * dtch146
      ch147 = chtab(147, itemp) + dtemp * dtch147
      ch153 = chtab(153, itemp) + dtemp * dtch153
      ch154 = chtab(154, itemp) + dtemp * dtch154
      ch155 = chtab(155, itemp) + dtemp * dtch155
      ch156 = chtab(156, itemp) + dtemp * dtch156
      ch157 = chtab(157, itemp) + dtemp * dtch157
      ch158 = chtab(158, itemp) + dtemp * dtch158
      ch159 = chtab(159, itemp) + dtemp * dtch159
      ch162 = chtab(162, itemp) + dtemp * dtch162 ! recomb electron energy loss factor for HeI
c      print*, "[Chemistry rate_eq]: ch160 =", ch160
c

c H3+ dissociative recombination: branching fractions from UMIST06,
c based on MAC04
c
      ch106a = (2.34d0 / 6.7d0) * ch106   ! H2 + H
c
c C+ + H2 -> products: rate and branching ratio from GOW16, based on Wakelam et al 2010
c
      ch155a = (2.31d0 / 3.3d0) * ch155
#endif

c
c H2 collisional dissociation rates -- need special handling because of
c density dependence
c
      ncrinv   = (2d0 * abh2 * (ch6 - ch5) + ch5)
      h2var0   = 1d0 / ( 1d0 + yn * ncrinv)
      h2var1   = ch3**h2var0
      h2var2   = ch4**h2var0
      ykdh     = ch1 * h2var1
      ykdh2    = ch2 * h2var2

c
c Photoelectric effect:
c
      G_dust = G0 * fshield_dust
c
c Photodissociation -- diffuse field
c
      AV_mean = (-1.0/2.5) * log(fshield_dust)

      ykdis = phtab(1) * fshield_H2 * fshield_dust_H2 * G0
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
c      gamma_co = G0 * (1.235e-10) * fshield_dust * fshield_CO
      gamma_co  = G0*phtab(10)*exp(-1.0*phgam(10)*AV_mean)*fshield_CO
      gamma_chx =  2.94d-10 * G_dust
#endif
c
c Hydrogen ionization from point sources -- ***TODO***
c

      uv_ion  = 0e0
c      uv_heat = 0e0
      uv_ion_HeII = 0e0
      uv_ion_HeIII = 0e0
      uv_ion_cI_directly_to_c2p = 0e0
      Lyc_HeatingRate_HII = 0e0
      Lyc_HeatingRate_HeII = 0e0
      Lyc_HeatingRate_HeIII = 0e0
#ifdef TR_ONTHESPOT
      do I = 1, nrrec
        if (I==iHII) then
          uv_ion  = Lyc_uv_ion(iHII)
          Lyc_HeatingRate_HII = Lyc_HeatingRate(iHII)
        else if (I==iHeII) then
          uv_ion_HeII  = Lyc_uv_ion(iHeII)
          Lyc_HeatingRate_HeII = Lyc_HeatingRate(iHeII)
        else if (I==iHeIII) then
          uv_ion_HeIII  = Lyc_uv_ion(iHeIII)
          Lyc_HeatingRate_HeIII = Lyc_HeatingRate(iHeIII)
        endif
      enddo 
#endif
      uv_ion_cp = uv_ion_HeII * (4.015D-18 /5.729D-18) !umrechnen He->cp  sigma_eff for BB(40.000K)
      uv_ion_cI = uv_ion * (10.215D-18 /3.344D-18) ! effective sigmas for BB with T=40.000 from 13.6 to 100 eV
c      if (uv_ion_HeII > 1d-10) then
c        uv_ion_cI_directly_to_c2p = uv_ion_cI
c        uv_ion_cI = 0d0
c      endif
c      call calc_uv(Lyc_EnergyRate,Lyc_PhotonExceEnergy, uv_ion, uv_heat)

#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      ph2   = phtab(2) * exp(-1.0 * phgam(2) * AV_mean)  * G0
      ph10  = phtab(10) * fshield_CO*exp(-1.0*phgam(10)*AV_mean) * G0
      ph48  = phtab(48) * exp(-1.0 * phgam(48) * AV_mean) * G0
      ph49  = phtab(49) * exp(-1.0 * phgam(49) * AV_mean) * G0
      ph50  = phtab(50) * exp(-1.0 * phgam(50) * AV_mean) * G0
      ph51  = phtab(51) * exp(-1.0 * phgam(51)  * AV_mean) * G0  !! For Si
#endif
c
c X-ray ionization and heating
c
#ifdef XRAY_TR
c ** JMackey **
c Replace Simon's x-ray routine with JM's heating/ionization rates.
c
c Ionization+heating Rate Function:
c 1 - input flux (attenuated), from rpar(nrpar+2*nrrec-->nrpar+2*nrrec+ipar(3)).
c 2 - Hydrogen nuclei number density (mass_density*X_H/m_p), where X_H
c    is the mass fraction of hydrogen.  (this is yn)
c 3 - dR:  path length of ray through cell (cm); only needed if cell is
c    optically thick.
c (if tempxsec) - gas temperature
c 4 - input abundances of xab(H^0,H2,He^0,He^+)
c 5 - electron fraction = n(e)/(n(H)+2n(H2)+n(He))
c    This should(!) be abe/(abhp+abHI+2.0*abh2+abhe)
c
c 6 - output: ionization/diss. rates dy/dt (s^{-1}, already multiplied by abundance!)
c     rates xion[H->H+, H2->H2+, H2->2H, He->He+, He+->He2+]
c 7 - output: total heating rate (erg/cm^3/s)
c 8 - output: total absorbed X-ray power (erg/s/H-nucleus)
c
      xray_ye = abe/(abhp+abHI+2.0*abh2+abhe)
      xab(1) = abHI
      xab(2) = abh2
      xab(3) = abheI
      xab(4) = abhep
      call ionization_heating_rates(xray_Fin,yn,dl,
#ifdef XRAY_TEMPXSEC
     $                             temp,
#endif
     $                             xab,xray_ye,xion,GammaX,HX)
c ** JMackey **
#endif


c
c Rates for recombination on grain surfaces. These rates must be
c computed here because of their dependence on the electron number
c density.
c
      if (abe .eq. 0d0) then
c If the fractional ionization is zero, then there won't be any recombination,
c so the value we use for phi doesn't matter too much -- 1d20 is simply an
c arbitrary large number
c
        phi = 1d20
      else
        phi = G_dust * sqrt(temp) / (yn * abe)
      endif

      phi = max(1.0e-40, phi)
c
c HII:
c
      if (phi .lt. 1e-6) then
        h_gr = 1.225d-13 * dust_to_gas_ratio
      else
        hgrvar1  = 8.074d-6 * phi**1.378d0
        hgrvar2  = (1d0 + ch34 * phi**ch35)
        h_gr     = 1.225d-13 * dust_to_gas_ratio /
     $             (1d0 + hgrvar1 * hgrvar2)
      endif
c Exponentially reduce rate at high temperature to account for fact that we
c expect hot gas to be largely dust free. This also prevents us from getting
c an anomalously large grain recombination cooling rate when T is large
c and G_dust is very small or zero
c
      h_gr = h_gr * exp(-temp/3.0e4)


!#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
c ** JM **  GOW17 dust recombination (H+ already included above).
      khe_dust = 5.572e-14 * dust_to_gas_ratio
      kcp_dust = 45.58e-14 * dust_to_gas_ratio
      kmp_dust = 2.166e-14 * dust_to_gas_ratio
      if (phi .gt. 1e-6) then
        hgrvar1 = 3.185e-7*phi**1.512
        hgrvar2 = 1.0 + ch135 * phi**ch136
        khe_dust = khe_dust / (1.0 + hgrvar1*hgrvar2)

        hgrvar1 = 6.089e-3*phi**1.128
        hgrvar2 = 1.0 + ch36 * phi**ch37
        kcp_dust = kcp_dust / (1.0 + hgrvar1*hgrvar2)

        hgrvar1 = 5.678e-8*phi**1.874
        hgrvar2 = 1.0 + ch38 * phi**ch39
        kmp_dust = kmp_dust / (1.0 + hgrvar1*hgrvar2)
      endif

      khe_dust = khe_dust * exp(-temp/3.0e4)
      kcp_dust = kcp_dust * exp(-temp/3.0e4)
      kmp_dust = kmp_dust * exp(-temp/3.0e4)
!endif
c
c H2 formation rate on dust -- we calculate this here, as this allows us to
c vary tdust without needing to call cheminmo again
c
      iflag_contcool = ipar(7)

#ifdef XRAY_TR
      if ((yn .gt. 1d9 / dust_to_gas_ratio**1.5d0)
     $ .and.(iflag_contcool.eq.0)) then
        tdust = temp
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, HX,
     $                      G0RP, dl, 1)
      else
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, HX,
     $                      G0RP, dl, 0)
      endif
#else
      if ((yn .gt. 1d9 / dust_to_gas_ratio**1.5d0)
     $ .and.(iflag_contcool.eq.0)) then
        tdust = temp
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust,
     $                      G0RP, dl, 1)
      else
        call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust,
     $                      G0RP, dl, 0)
      endif
#endif
c
c Cosmic ray ionization rates: precomputed in cheminmo.
c
c HI:
      cr1  = crtab(1)*cr_rescale
      cr2  = crtab(2)*cr_rescale   ! H2 -> H2+ + e-
      cr10 = crtab(10)*cr_rescale  ! H2 -> H + H+ + e
      cr11 = crtab(11)*cr_rescale  ! H2 -> H + H

#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
c HeI:
      cr6 = crtab(6)*cr_rescale
c CI:
      cr3 = crtab(3)*cr_rescale
c CO: CO + CR -> C + O
      cr13 = crtab(13)*cr_rescale
      crco = cr1 * 6.52*cr_rescale
c
c O + H3+ -> products; rate and branching ratio from GOW16, rate based
c on de Ruette et al 2015. Tiny term in denominator is to avoid division
c by zero in the case where abh2 and abe are both zero
c
      rohx = cst(92)*abh2/(cst(92)*abh2 + ch154 * abe + 1e-40)

c CR-induced photochemistry
c C + gamma_CR -> C+ + e-
      cp1  = crphot(1)*cr_rescale
c CO + gamma_CR -> C + O (JM: use UMIST12 rate)
c do not apply cr_rescale scaling factor here, already included in cr1
      cp11 = cr1 *420.0*(temp/300.0)**1.17 * abco
      cpSi = cr1 * 8460.0 * abmI  !! UMIST12 [*2 (1-w)] [*2 (per H, not per H2)]
c
c O+ equilibrium abundance
c
c Formation rate due to charge transfer with H+ depends on O abundance,
c which is not yet known. However, if we add ch47 * abhp * abop to both
c sides of balance equation, can remove this dependence, as l.h.s. now
c depends only on abotot, which is known. The other advantage this has
c is that it allows us to avoid a division-by-zero error if abHI and abh2
c are both zero. Instead, we simply find that abop=abotot and aboI=0, which
c is correct for fully ionized gas
c
      num   = ch47 * abotot * abhp + ch158 * abhep * abohx
      denom = cst(93) * abH2 + ch46 * abHI + ch47 * abhp
c
      abop = num / denom
      aboI = abotot - abop
      abe  = abe + abop
      if (ISNAN(abe)) then
        print *, "abe nan after abop"
        print *, "abop must be NaN! Terms in num, denom:"
        print *, "Num: ", ch47, abotot, abhp, ch158, abhep, abohx
        print *, "Denom: ", cst(93), abH2, ch46, abHI, ch47, abhp
      endif
c
c H3+ equilibrium abundance
c
c ** JM **
c calculate fraction of H2+ that reacts in each of the three channels:
c fh3p = fraction that reacts with H2 to produce H3+
c fh2  = fraction that undergoes charge exchange with H --> H2 + H+
c f2h  = fraction that undergoes dissociative recombination --> 2H
      fh3p = ch105 * abh2  / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
      fh2  = cst(1) * abHI / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
      f2h  = ch30 * abe    / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
c
c ** JM ** Xray creation of H2+ in numerator
      num   = fh3p * ( (cr2 / yn) * abh2 + cst(75) * abhep * abh2
     $                 + ch24 * abhp * abh2
#ifdef XRAY_TR
     $                 + xion(irh2)/yn
#endif
     $ )

      denom = ch159 * abcI + ch153 * aboI + cst(59) * abco
     $      + ch106 * abe + cst(81) * abmI
c
      abh3p = num / denom
      if (abh3p .gt. 0.333) then
        print*, "ERROR! abh3p=",abh3p, "num,denom=",num,denom
        abh3p = 0.333
      endif
      abe   = abe + abh3p
      if (ISNAN(abe)) then
        print *, "abe nan after h3p"
        print *, "abh3p must be NaN! Terms in num, denom:"
!        print *, "Num terms: ", fh3p, cr2, yn, abh2, cst(75),
!     $    abhp, ch24, abhp, xion(irh2)
!        print *, "Ionization stuff: ", HX, xray_fin, temp, xion
      endif
#endif
c
c Protons
c
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb. (case B)
     $          -  h_gr * abhp  * 1d0
     $          + cst(75) * abhep * abh2 * (fh2)
     $          -  ch24 * abhp  * abh2 * (1.0 - fh2) ! CT: H+ + H2 -> H2+ + H
     $          +  ch9  * abhep * abh2) * yn    ! DCT: He+ + H2
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
     $          + cr10  * abh2
     $          + (fh2) * cr2 * abh2      ! H2+ + H -> H2 + H+
     $          + crco  * abco * abhI / (abhI + abh2 + 1.0e-40) ! GOW17 T2.10, mod SCOG
#ifdef XRAY_TR
     $          + xion(irh )              ! X-ray 2ndary ioniz. of H^0
     $          + xion(irh2) * (fh2)      ! X-ray 2ndary ioniz. of H_2
     $          + (xion(irh )+xion(irh2))*2.97*abohx/(abhi+abh2)  ! 2ndary ioniz. of OHx
     $          + (xion(irh )+xion(irh2))*3.92*abhcop/(abhi+abh2) ! 2ndary ioniz. of HCO+
#endif
     $          +  ph50 * abhcop                ! Photodissoc. of HCO+
     $          +  uv_ion * abHI
c
c Molecular hydrogen
c
      ydot(ih2) = (
     $  -   ch25  * abe   * abh2          ! Collisional dissoc.: H2 + e-
     $  -   ykdh  * abHI  * abh2          ! Collisional dissoc.: H2 + H
     $  -   ykdh2 * abh2  * abh2          ! Collisional dissoc.: H2 + H2
     $  +   Rdust * abHI  * 1d0
     $  -   ch24  * abhp  * abh2 * (1.0 - fh2)  ! CT: H2 + H+ -> H2+ + H
     $  - ch155   * abh2  * abcp          ! GOW17 T1.10
     $  + ch159   * abh3p * abcI          ! GOW17 T1.01
     $  + ch153   * abh3p * aboI * rohx   ! GOW17 T1.02
     $  - cst(93) * abop * abh2           ! GOW17 T1.04
     $  + ch157   * abchx * abhI          ! GOW17 T1.26
     $  - cst(54) * abcI  * abh2
     $  + cst(59) * abh3p * abco
     $  -     ch9 * abhep * abh2
     $  - cst(75) * abhep * abh2 * (fh3p +f2h)  !! including diss.rec. ** JM **
     $  +  ch106a * abh3p * abe
     $  + cst(81) * abh3p * abmI) * yn
     $  -     cr2 * abh2  * (2.0*fh3p +f2h)   ! Cosmic rays ! including diss.rec.
     $  -   ykdis * abh2                      ! Photodissoc:         H2 + PH
     $  -   cr10  * abh2
     $  -   cr11  * abh2
     $  -    crco * abco * abh2 / (abhI + abh2 + 1.0e-40)  ! GOW17 T2.10, mod SCOG
#ifdef XRAY_TR
     $  - xion(irh2) * (fh3p +f2h)      ! X-ray 2ndary ioniz. of H_2
     $  - xion(drh2)                    ! X-ray: 2ndary diss. by e-
#endif
     $  - uv_ion * abh2
c
c C+
c
      ydot(icp)  = (
     $      ch92 * abhep * abco
     $  -   ch40 * abcp  * abe      ! C+ recombination
#if CHEMISTRYNETWORK == 17
     $  +  ch160 * abc2p * abe      ! C++ recombination
#endif
     $  -kcp_dust * abcp
     $  - ch155  * abcp  * abh2
     $  - ch156  * abcp  * abohx    ! GOW17 T1.11
     $  +   ch43 * abcI  * abe
     $                          ) * yn
     $  +    ph2 * abcI
     $  +    cr3 * abcI
     $  +    cp1 * abcI
     $  + uv_ion_cI * abcI          ! photo-ionization C -> C+
#if CHEMISTRYNETWORK == 17
     $  - uv_ion_cp * abcp          ! photo-ionization C+ -> C++
#endif
#ifdef XRAY_TR
     $  +    (xion(irh )+xion(irh2))*3.92*abcI/(abhi+abh2)
     $  +    (xion(irh )+xion(irh2))*3.92*abco/(abhi+abh2)
     $  +    (xion(irh )+xion(irh2))*3.92*abchx/(abhi+abh2)
     $  + (xion(irh )+xion(irh2))*3.92*abhcop/(abhi+abh2)
     $  +    xion(irh )*1020.0*abcI
#endif
c
c CHx
c
      ydot(ichx) = (
     $           +  cst(54) * abcI * abh2
     $           -   ch157  * abchx * abhI  ! T1.26
     $           +   ch155a * abcp  * abh2  ! T1.09
     $           +   ch159  * abcI * abh3p  ! T1.01
     $           -  cst(85) * aboI * abchx) * yn
     $           -     ph48 * abchx
#ifdef XRAY_TR
     $           - (xion(irh2)*730.0 + xion(irh )*35.0)*2.0*abchx     !! XR-PHOT
     $           -    (xion(irh )+xion(irh2))*3.92*abchx/(abhi+abh2)  !! XR-2ndary
#endif

c
c OHx
c
      ydot(iohx) = (
     $           + ch153    * abh3p * aboI * rohx  ! GOW17 T1.02
     $           + cst(93)  * abop  * abh2 * rohx  ! GOW17 T1.04
     $           - ch156    * abcp  * abohx        ! GOW17 T1.11
     $           - cst(94)  * aboI  * abohx        ! GOW17 T1.27
     $           - ch158    * abhep * abohx        ! GOW17 T1.29
     $           -    ch147 * abcI * abohx
     $                                       ) * yn
     $           -     ph49 * abohx
#ifdef XRAY_TR
     $           - (xion(irh2)*508.0 + xion(irh )*87.6)*2.0*abohx     ! XR-PHOT
     $           -    (xion(irh )+xion(irh2))*2.97*abohx/(abhi+abh2)  ! XR-2ndary
#endif

c
c CO
c

      ydot(ico) = (cst(85) * abchx  * aboI
     $          +    ch147 * abohx  * abcI
     $          +     ch56 * abhcop * abe
     $          -  cst(59) * abh3p  * abco
     $          -     ch92 * abhep  * abco) * yn
     $          -     ph10 * abco
     $          +     ph50 * abhcop
     $          -     cr13 * abco
     $          -     cp11
     $          -     crco * abco * (1.-abhI / (abhI + abh2 + 1.0e-40))    ! GOW17 T2.10, mod SCOG
     $          -   uv_ion*3.861003861 * abco     ! photo-ionization
c      3.8.. factor from ratio of CO to H2 dissociation rate in 1 draine
#ifdef XRAY_TR
c     $          - xion(irh2)*2.7*sqrt(temp*abco*1.0e-3)          ! XR-PHOT  MHT96
     $          - xion(irh2) *420.0*(temp/300.0)**1.17 * abco     ! XR-PHOT  UMIST12
     $          - (xion(irh )+xion(irh2))*3.92*abco/(abhi+abh2)   ! XR-2ndary
#endif

c
c HCO+
c
      ydot(ihcop) = (cst(59) * abh3p  * abco
     $             + ch156   * abcp  * abohx    ! GOW17 T1.11
     $             -    ch56 * abhcop * abe) * yn
     $             -    ph50 * abhcop
     $             +   crco  * abco * abh2 / (abhI + abh2 + 1.0e-40)  ! GOW17 T2.10, mod SCOG
#ifdef XRAY_TR
     $             - (xion(irh )+xion(irh2))*3.92*abhcop/(abhi+abh2)  ! XR-2ndary
#endif

c
c M+
c
      ydot(imp) = (cst(81) * abh3p * abmI       !! H3+ charge exchange
     $          +     ch44 * abmI  * abe        !! collisional ionization
     $          - kmp_dust * abmp               !! grain recombination
     $          -    ch146 * abmp  * abe) * yn  !! radiative recombination
     $          +     ph51 * abmI               !! UV photoionization
     $          + cpSi                          !! Si + CRPHOT --> Si+ + e
#ifdef XRAY_TR
     $          + xion(irh )*6.67*abmI          !! XR-2ndary (Si)
     $    + (xion(irh2)*4320.0 + xion(irh)*1050.0) * 2.0 * abmI  !! XR-PHOT (Si)
#endif

c
c He+
c
      ydot(ihep) =  cr6  * abheI
     $      + (- ch9 * abhep * abh2
     $         - ch92 * abhep * abco
     $         - ch17 * abhep * abe
     $         + ch14 * abheI * abe
#if CHEMISTRYNETWORK == 17
     $         + ch20 * abhe2p * abe   ! HeIII recombination
#endif
     $         - ch158  * abhep * abohx   ! GOW17 T1.29
     $         - khe_dust * abhep
     $         - cst(75) * abhep * abh2 ) * yn
     $         + uv_ion_HeII * abheI    ! Photo-ionization of HeI ->HeII
     $         - uv_ion_HeIII * abhep   ! Photo-ionization HeII -> HeIII
#ifdef XRAY_TR
     $         + xion(irhe)                             ! XR-2ndary
#endif
#endif !chemnet15 or chemnet17

#if CHEMISTRYNETWORK == 17
c
c C++
c
      ydot(ic2p)  = (
     $  -   ch160 * abc2p * abe    ! C++ recombination
     $                         ) * yn
     $  + uv_ion_cp * abcp         ! photo-ionization C+ -> C++
     $  + uv_ion_cI_directly_to_c2p * abcI ! photo-ionization C -> C++
c
c HeIII
c
      ydot(ihe2p) = (
     $  -    ch20 * abhe2p * abe   ! HeIII recombination
     $                         ) * yn
     $  + uv_ion_HeIII * abhep     ! photo-ionization HeII -> HeIII
#endif !chemnet17
c
c
c
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      ydot(ihp) = (ch11 * abe   * abHI          ! Collisional ion: H  + e-
     $          -  ch13 * abhp  * abe           ! Gas-phase H recomb.(case B)
     $          -  h_gr * abhp  * 1e0) * yn     ! Grain H recomb.
     $          +  cr1  * abHI                  ! Cosmic ray ion.:  H + CR
c     $          + xr_ion* abHI                  ! XRay ionization
c     $          + 2.*uv_ion* abh2               ! Photo-ionization of H2
     $          + cr10  * abh2                  ! CR ionization of H2
     $          + uv_ion* abHI
c WHY is photoionization of HI missing? uv_ion *abHI?
c
c Molecular hydrogen
c
      ydot(ih2) = (-ch25  * abe  * abh2          ! Collisional dissoc.: H2 + e-
     $          -   ykdh  * abHI * abh2          ! Collisional dissoc.: H2 + H
     $          -   ykdh2 * abh2 * abh2          ! Collisional dissoc.: H2 + H2
     $          +   Rdust * abHI * 1e0 ) * yn    ! H2 formation on dust
     $          -   ykdis * abh2                 ! Photodissoc:         H2 + PH
     $      - 0.5e0 * cr2 * abh2   ! Fudge factor of two gives correct HI production
     $          -   cr10  * abh2
     $          -   cr11  * abh2
     $          -   uv_ion* abh2

#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
c
c CO -- from NL97
c
      k0 = 5e-16
      k1 = 5e-10
      abo = max(0.0, abundo - abco)
      if (abo .eq. 0.0) then
        beta = 0.0
      elseif (gamma_chx .lt. 1e-30) then
        beta = 1.0
      else
        beta = k1 * abo / (k1 * abo +
     $       gamma_chx / yn)
      endif

c the h2phih dissociation rate is a bit iffy it comes from ionising
c radiation and how that interacts with CO is not modeled
c the uvphih is fine though, the factor comes from the comparison
c the dissociation rates of H2 and CO at UV

      ydot(ico) = k0 * abcp * abh2 * beta * yn
     $          - gamma_co * abco
c Bug: already contained in gamma_co
c     $          - ykdis*3.861003861 * abco
c Here we only consider ionisation from the stellar radiation field but neglect dissociation. The photons below 13.6 eV are omitted. DS 02/19
     $          - uv_ion*3.861003861 * abco
c 3.8.. factor from ratio of CO to H2 dissociation rate in 1 draine
c field, 2.0e-10/5.18e-11, Roellig et al 2007 numbers
c this is an approximation as we do not know the actual dissocation rates at > 13.6 eV
#if CHEMISTRYNETWORK == 6
c
c C++
c
      ydot(ic2p)  = (
     $  -   ch160 * abc2p * abe    ! C++ recombination
     $                         ) * yn
     $  + uv_ion_cp * abcp         ! photo-ionization C+ -> C++
#endif
#endif
#endif

c Commented out, as very problem specific hack. DS 01/2019
cc ** JM **
cc If T>2e4 K, then send molecular abundances to zero on 1e5 year timescale.
cc if T>1e5 K, then send ionic species to full ionization on 1e5 year timescale.
c      if (temp .gt. 2.0e4) then
c        ydot(ih2)   = (1e-40-abh2  )/3.156e12
c        ydot(ico)   = (1e-40-abco  )/3.156e12
c        ydot(ichx)  = (1e-40-abchx )/3.156e12
c        ydot(iohx)  = (1e-40-abohx )/3.156e12
c        ydot(ihcop) = (1e-40-abhcop)/3.156e12
c      endif
c      if (temp .gt. 1.0e5) then
c        ydot(ihp)   = (1.0-abhp)/3.156e12
c        ydot(ihep)  = (abhe   -abhep)/3.156e12
c        ydot(imp)   = (abundm - abmp)/3.156e12
c        ydot(icp)   = (abundc - abcp)/3.156e12
c      endif
cc ** JM **


      if (no_cooling .eq. 1) then
        ydot(itmp) = 0.0
      else
c
c H2 collisional dissociation cooling:
c
      rates_chem(1) = 4.48d0 * eV * (ykdh  * abHI * abh2 +
     $                               ykdh2 * abh2 * abh2 +
     $                               ch25  * abe  * abh2 ) * yn**2
c
c H2 photodissociation heating -- 0.4eV per photodissociation:
c
      rates_chem(2) = - 4d-1 * eV * abH2 * yn * ykdis
c
#ifdef XRAY_TR
c ** JM ** Heating due to ionization and dissociative recomb of H2
c  H_2 + Xray --> H_2^+ +e^- --> 3 further pathways (MS05)
c  The three rates are in MS05, App.B.2.
      rates_chem(2) = rates_chem(2) -
     $ (10.9*f2h + 0.94*fh2 + 8.6*fh3p) * xion(irh2)*eV *yn
      rates_chem(2) = rates_chem(2) -  4d-1 * eV * xion(drh2) *yn
#endif

c
c Heating by UV pumping of H2 -- effective for n > n_crit
c
c We assume 2 eV per UV photon (Burton, Hollenbach & Tielens 1990),
c and adopt a pumping rate 8.5 times larger than the dissociation rate
c (Draine & Bertoldi 1996, table 2).
c
      cdex = (1.4 * exp(-1.81d4 / (temp + 1200)) * abH2
     $     +  1.0 * exp(-1d3 / temp) * abHI) * 1d-12 * sqrt(temp) * yn
      rdex = 2d-7
        fpheat = cdex / (cdex + rdex)

      rates_chem(3) = - 2d0 * eV * abH2 * yn * ykdis * 6.94e0 * fpheat
c
c H2 formation heating -- we assume a similar depedence on n / n_crit to
c the collisional dissociation rates, although ideally one would like to
c do a more careful, level-dependent treatment.
c
c We assume that a fraction f_ex of the energy goes into rotational &
c vibrational excitation of the H2, and a further fraction f_kin goes
c directly into kinetic energy (i.e. heat). The remaining energy goes
c directly into the dust
c
c Dust -- 4.48eV
c
        if (RH2 .gt. 0e0) then
c        call compute_h2_heating_fraction(temp, yn, abh2, f_ex)
c        f_ex = f_ex * h2_form_ex
        f_ex = 0d0
        else
        f_ex = 0d0
        endif
        f_tot = h2_form_kin + f_ex
c
      rates_chem(4) = - 4.48d0 * eV * f_tot * RH2
c
c HI collisional ionization cooling -- 13.6eV per ionization:
c
      rates_chem(5) = 13.6d0 * eV * ch11 * abe * abHI * yn**2
c
c HII recombination cooling -- we assume that kT ergs is radiated
c per recombination and that this is independent of whether the
c recombination occurs in the gas phase or on a grain surface
c
c SV 5.12.23: adding recombination energy loss factors to H and He
c recombinations (~2/3 kT are lost per gas phase recombination)
c
      rates_chem(6) = kboltz * temp * abhp
     $              * (ch13 * abe * ch161 + h_gr)*yn**2
#if CHEMISTRYNETWORK == 15 ||  CHEMISTRYNETWORK == 17
     $              + kboltz * temp * abhep
     $              * ch17 * abe * ch162 * yn**2
#endif
c
c Photoionization heating
c
c        rates_chem(7) = uv_heat * abHI * yn
c         photon_flux = uv_ion/ dl / (6.3e-18 *
c     $   (13.6/(Lyc_PhotonExceEnergy/1.602e-12 + 13.6)) ** 3.)
c         rates_chem(7) = -1.e0*min(ch13*yn**2.e0,photon_flux)
c     $                 *max(Lyc_PhotonExceEnergy - 3./2.*kboltz*temp,0.)
      rates_chem(7) = - Lyc_HeatingRate_HII * abHI * yn
#if CHEMISTRYNETWORK == 15  || CHEMISTRYNETWORK == 17
     $              - Lyc_HeatingRate_HeII * abheI * yn
     $              - Lyc_HeatingRate_HeIII * abhep * yn
#endif
c
c X-ray heating
c Add Coulomb heating rate from X-ray module.
#ifdef XRAY_TR
      rates_chem(8) = -GammaX
#endif
c
      ylam_chem = rates_chem(1)  + rates_chem(2)  + rates_chem(3)  +
     $            rates_chem(4)  + rates_chem(5)  + rates_chem(6)  +
     $            rates_chem(7)
#ifdef XRAY_TR
     $          + rates_chem(8)
#endif
      if (ylam_chem .ne. ylam_chem) then
        print *, "chem rates: ", rates_chem
      endif

c
c Initialize all abundances to zero, then just set the non-zero ones
c
        do i = 1, nabn
        abundances(i) = 0d0
        enddo
c
        abundances(1)  = abh2
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
        abundances(2)  = max(0e0, abundo - abco) * max(0e0, 1.0-abhp)
        abundances(5)  = abco
        abundances(7)  = abcp
#else
        abundances(2)  = abundo * max(0e0, 1.0-abhp)
        abundances(7)  = abundc
#endif
#if CHEMISTRYNETWORK == 6
        abundances(18) = abc2p
#endif
        abundances(9)  = abundsi
        abundances(10) = abe
        abundances(11) = abhp
        abundances(12) = abHI
c Approximate OII abundance
        abundances(17) = abundo * abhp
#endif
#if CHEMISTRYNETWORK == 15 || CHEMISTRYNETWORK == 17
        abundances(2)  = aboI
        abundances(5)  = abco
        abundances(6)  = abcI
        abundances(7)  = abcp
        abundances(8)  = abmI
        abundances(9)  = abmp
        abundances(10) = abe
        abundances(11) = abhp
        abundances(12) = abHI
        abundances(14) = abheI
        abundances(15) = abhep
        abundances(17) = abop
#endif
#if CHEMISTRYNETWORK == 17
        abundances(16) = abhe2p
        abundances(18) = abc2p
#endif
c
c Compute local Jeans length (for use in cool_func)
c
        rho = (1e0 + 4e0 * abhe) * mh * yn
        cs2 = (gamma * (gamma - 1e0) * energy / rho)
        L_jeans  = sqrt(pi * cs2 /
     $             (guniv * rho))
c

        if (uv_ion .gt. 1d-99) log_HII_region = 1

        !print *, "DUMP 2a", temp, yn, abundances
        call cool_func(temp, yn, L_jeans, divv, G_dust, abundances,
     $                  ylam, rates, 0, log_HII_region, cr_rescale,
     $                  iflag_contcool)
c
        ydot(itmp) = -(ylam + ylam_chem)
c        print *, temp, ydot(itmp), ylam, ylam_chem
c        if (temp < 2.0) then
c          print *, "SMALL TEMPERATURES! Here's Tdot: "
c          print *, "ydot = ", ydot
c          print *, "ylam = ", ylam
c          print *, "ylam_chem = ", ylam_chem
c          print *, "ipar(3) rpar = ", ipar(3), rpar
c          print *, "X-ray flux = ", rpar(14:)/(2.998E10)
c          print *, "rates_chem = ", rates_chem
c        endif
        if (ydot(itmp) .ne. ydot(itmp)) then
          print *, "NAN YDOT(ITMP) "
          print *, "ydot = ", ydot
          print *, "ylam = ", ylam
          print *, "ylam_chem = ", ylam_chem
          print *, "rates_chem = ", rates_chem
          print *, "ipar(3) rpar = ", ipar(3), rpar
#ifdef XRAY_TR
          print *, "X-ray flux = ", rpar(14:)/(2.998E10)
          print *, "Xray GammaX, HX: ", GammaX, HX
#endif
        endif
        !print *, "DUMP 2b", temp, yn, abundances, ydot(itmp)

c#ifdef XRAY_TR
c ** JM ** X-rays heat ionized gas a lot, so limit heating for T>1e6K
c      if (temp .gt. 1.0e5 .and. ydot(itmp) .gt. 0.0) then
c        ydot(itmp) = ydot(itmp) * exp(-(temp*1e-6)**3)
c      endif
c
c#endif


#ifdef THERMAL_INFO_DUMP
c
cc We'll leave the accretion heating holder here for the future.
c
      heat_accr   = 0.
      heat_cray    = rates(11)
      heat_photo_elec = rates(12)
      heat_h2_form = rates_chem(4)
      heat_h2_pdiss = rates_chem(2)
      cool_h2_line = rates(2)
      cool_cplus   = rates(16)
      cool_oxy     = rates(13)
      cool_co      = rates(7) + rates(21) + rates(22) +
     $               rates(8) + rates(9) + rates(24)
      cool_gas_dust = rates(1)
      cool_c        = rates(14)
      cool_Lyman_alpha = rates(25)
      cool_bremsstrahlung = rates(28)
#endif
#ifdef THERMAL_DUMP_ALL_RATES
      cool_rates(:) = rates(:)/yn
      cool_rates_chem(:) = rates_chem(:)/yn
#endif

      endif ! no cooling
      call validate_output(y,ydot,rpar)

      !print*, "rate_eq:"
      !print*, "rate_eq: y   ", y
      !print*, "rate_eq: ydot", ydot
c
c Now add the pdV term:  Actually, let's not! Gadget evolves entropy... :)
c
c      ydot(itmp) = ydot(itmp) - (gamma - 1d0) * energy * divv
c
! #endif /* CHEMCOOL */
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////              R A T E _ E Q                \\\\\\\\\\
c
c=======================================================================
! #endif /* not KOYAMA */
! #endif /* CHEMISTRYNETWORK == 15 */




#ifdef XRAY_TR
c=======================================================================
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////             GET_ABH3P                 \\\\\\\\\\
c
c=======================================================================
c
c - J.Mackey: 2016.10.05: turned this into a function.
c             2017.11.29: Include O+ abundance; made consistent with rate_eq
c
      subroutine get_abh3p(y, temp, yn, xion, fh3p, abe, abop)
      implicit NONE
#include "Flash.h"
#include "cool.h"
#include "non_eq.h"
#include "cons.h"
      real y(nspec), xion(5)
      real temp, dtemp, yn
      real fh3p, fh2, f2h, num, denom, denom2
      real abe, abhi, abci, aboi, abotot, abop, abmi
      real abco, abh2, abhcop, abhep, abhp, abohx
      real old_abe, new_abe, old_abh3p, abh3p
      real ch30, dtch30, ch105, dtch105, ch106, dtch106
      real ch153, dtch153, ch158, dtch158, ch159, dtch159
      REAL ch46, ch47, dtch46, dtch47, cr2
      integer ct, itemp
      INTEGER irh, irh2, drh2, irhe, irhep
      parameter (irh   = 1)
      parameter (irh2  = 2)
      parameter (drh2  = 3)
      parameter (irhe  = 4)
      parameter (irhep = 5)
c
      abh2   = max(1e-40, y(ih2))
      abhp   = max(1e-40, y(ihp))
      abco   = max(1e-40, y(ico))
      abhcop = max(1e-40, y(ihcop))
      abohx  = max(1e-40, y(iohx))
      abhep  = max(1e-40, y(ihep))

      abe   = max(y(ihp) + y(icp) + y(ihcop) + y(ihep) + y(imp)
     $            + 2. * y(ic2p), 1e-40)
      abHI  = max(1d0 - 2.*y(ih2) -y(ihp) -y(ichx) -y(iohx) -y(ihcop),
     $            1e-40)
      abcI  = max(abundc - y(icp) - y(ichx) - y(ico) - y(ihcop)
     $            - y(ic2p), 1e-40)
      abotot = max(abundo - abohx - abco - abhcop, 1e-40)
      abmI  = max(abundM - y(imp), 1e-40)
      cr2  = crtab(2)
c
c Compute temperature
c
      if (temp .le. 0d0) then
        itemp   = 1
        dtemp   = 0d0
        temp    = tmin
      elseif (temp .le. tmin) then
        itemp = 1
        dtemp = 0d0
      elseif (temp .ge. tmax) then
        itemp = nmd
        dtemp = 0d0
      else
        itemp = int(dlog10(temp) / dtlog) + 1
        dtemp = temp - temptab(itemp)
      endif
      dtch30 = dtchtab(30, itemp)
      dtch46 = dtchtab(46, itemp)
      dtch47 = dtchtab(47, itemp)
      dtch105 = dtchtab(105, itemp)
      dtch106 = dtchtab(106, itemp)
      dtch153 = dtchtab(153, itemp)
      dtch158 = dtchtab(158, itemp)
      dtch159 = dtchtab(159, itemp)
      ch30 = chtab(30, itemp) + dtemp * dtch30
      ch46 = chtab(46, itemp) + dtemp * dtch46
      ch47 = chtab(47, itemp) + dtemp * dtch47
      ch105 = chtab(105, itemp) + dtemp * dtch105
      ch106 = chtab(106, itemp) + dtemp * dtch106
      ch153 = chtab(153, itemp) + dtemp * dtch153
      ch158 = chtab(158, itemp) + dtemp * dtch158
      ch159 = chtab(159, itemp) + dtemp * dtch159
c
c O+ equilibrium abundance
c
c Formation rate due to charge transfer with H+ depends on O abundance,
c which is not yet known. However, if we add ch47 * abhp * abop to both
c sides of balance equation, can remove this dependence, as l.h.s. now
c depends only on abotot, which is known. The other advantage this has
c is that it allows us to avoid a division-by-zero error if abHI and abh2
c are both zero. Instead, we simply find that abop=abotot and aboI=0, which
c is correct for fully ionized gas
c
      num   = ch47 * abotot * abhp + ch158 * abhep * abohx
      denom = cst(93) * abH2 + ch46 * abHI + ch47 * abhp
c
      abop = num / denom
      aboI = abotot - abop
      old_abe = abe
      abe  = abe + abop
      new_abe = abe
c
c H3+ equilibrium abundance
c
c ** JM **
c calculate fraction of H2+ that reacts in each of the three channels:
c fh3p = fraction that reacts with H2 to produce H3+
c fh2  = fraction that undergoes charge exchange with H --> H2 + H
c f2h  = fraction that undergoes dissociative recombination --> 2H
      fh3p = ch105 * abh2  / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
      fh2  = cst(1) * abHI / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
      f2h  = ch30 * abe    / (ch105 * abh2 + cst(1) * abHI + ch30 * abe)
c
c ** JM ** Xray creation of H2+ in numerator
      num   = fh3p * ( (cr2 / yn) * abh2 + cst(75) * abhep * abh2
     $                 + xion(irh2)/yn )
c
#ifdef GONG
      denom = ch159 * abcI + ch153 * aboI + cst(59) * abco
     $      + ch106 * abe + cst(81) * abmI
c      print*, "H3+ denom", ch159, abcI, ch153, aboI, cst(59),
c     $        abco, ch106, abe, cst(81), abmI
#else
      denom2 = cst(84) * abcI + cst(86) * aboI + cst(59) * abco
     $      + ch106 * abe + cst(81) * abmI
      denom = denom2
#endif
c
      abh3p = num / denom
      fh3p = abh3p
      abe   = abe + abh3p
c
c      print*, "electrons: ", old_abe, new_abe, abe, abh3p, denom, denom2
c     $        , ch159
c
      return
      end
c=======================================================================
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E        //////////
c    //////////             GET_ABH3P                 \\\\\\\\\\
c
c=======================================================================
#endif

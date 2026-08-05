
#include "Flash.h"

c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
c     Written by S. Glover (AMNH, 2004-2005, AIP 2006-2007)
c
c ** JM **
c add extra parameters for the X-ray flux.
c - Fin: input flux array (attenuated) erg/cm2/s/bin
c - Nbins: Number of energy bins
c - HX: Xray energy absorbed in cell per H-nucleus (erg/s).
c ** JM **
c
      subroutine evolve_abundances(time, dl, yn, divv, energy,
     $                             current_redshift, non_eq_abundances,
     $                             fshield_H2, fshield_CO,
     $                             AV_mean, chi_mean,
     $                             dust_temp,
     $                             G0RP, cr_rescale
#ifdef TR_ONTHESPOT
     $                             ,Lyc_uv_ion,Lyc_HeatingRate
     $                             ,tr_osIHII, tr_osIHeII, tr_osIHeIII
#endif
#ifdef XRAY_TR
     !! ** JM ** Added extra variables for X-ray flux.
     $                            ,Fin, Nbins, HX
#endif
     $                            ,nstep_max, istatus, turbostep
     $                            ,iccool)
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
#include "Flash.h"
c#include "xray.h"



c
      REAL    time, dl, yn, divv, energy, current_redshift
      REAL    non_eq_abundances(num_non_eq_species)
      REAL    dust_temp
#ifdef TR_ONTHESPOT
      REAL, dimension(nrrec) :: Lyc_HeatingRate, Lyc_uv_ion
      INTEGER tr_osIHII, tr_osIHeII, tr_osIHeIII
#endif
      REAL    G0RP
      REAL    cr_rescale

#ifdef XRAY_TR
c ** JM **  X-ray params
      INTEGER Nbins
      !real*8, dimension(xdr_nbins) :: Fin
      real*8, dimension(Nbins) :: Fin
      REAL*8    HX      ! X-ray absorbed power per H-nucleus (erg/s)
c ** JM **  X-ray params
#endif

      INTEGER nstep_max, istatus
      INTEGER turbostep
      INTEGER iccool
      LOGICAL nothrow

      REAL    abe, abhep, abHI
      REAL    yntot, ekn, temp, gamma
      REAL    t_start, t_end, dt_est, dt_est2, deltat, t_in, t_out
      REAL    dt_correction
      REAL    deltat_old
      REAL    tdust_old
      REAL    fshield_dust, fshield_H2, fshield_CO, AV_mean,
     $        chi_mean, fshield_dust_H2
      REAL    RH2, Rdust
      integer new_step
c
      integer nstep
      REAL    BACKUP_DATA(20)
c
      REAL safety(nspec)
      integer counterMinDt
      real minDtBackup, currDtMin
      integer used_for_stepsize_control(nspec)
#if CHEMISTRYNETWORK == 4
      DATA safety /1e-1, 1e-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 5
      DATA safety /1e-1, 1e-1, 1e-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 6
      DATA safety /1e-1, 1e-1, 1e-1, 1e-1, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1, 1, 1/
#endif

#if CHEMISTRYNETWORK == 15
      DATA safety /ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT,
     $ ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT,
     $             DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1, 1, 1, 1, 1, 1, 1, 1, 1, 1/
#endif
#if CHEMISTRYNETWORK == 17
      DATA safety /ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT,
     $ ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT, ABUND_DT,
     $          ABUND_DT, ABUND_DT, DTCOOL_SCALE_FACTOR/
      DATA used_for_stepsize_control /1,1,1,1,1,1,1,1,1,1,1,1/
#endif
c
      integer I, J
c
c DVODE control variables -- NB atol declared in non_eq.h
c

c      REAL    rtol(nspec)
      integer itask, iopt, itol, mf, istate, iflag_dvi
      integer have_set_dvode
c
c Other DVODE variables
c
      integer lrw, liw
      parameter (lrw   = 22 + 9*nspec + 2*nspec*nspec)
c      parameter (lrw   = 22 + 9*nspec + nspec*nspec)
      parameter (liw   = 30 + nspec)
#ifdef XRAY_TR
c ** JM ** Added Nbins to the number of real parameters.
      !real, dimension(nrpar+xdr_nbins)  :: rpar
      real, dimension(nrpar+2*nrrec+Nbins)  :: rpar
      real maxdy  ! test if anything is changing.
c ** JM **
#else
      REAL rpar(nrpar+2*nrrec)
      real maxdy
#endif
      REAL    rwork(lrw)
      integer ipar(nipar), iwork(liw)
      REAL    y(nspec), ydot(nspec), y_old(nspec)
c
      external jac, rate_eq

      BACKUP_DATA(1) = time
      BACKUP_DATA(2) = dl
      BACKUP_DATA(3) = yn
      BACKUP_DATA(4) = divv
      BACKUP_DATA(5) = energy
      BACKUP_DATA(6) = current_redshift
      BACKUP_DATA(7) = non_eq_abundances(ih2)
      BACKUP_DATA(8) = non_eq_abundances(ihp)
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      BACKUP_DATA(9) = non_eq_abundances(ico)
#endif
#if CHEMISTRYNETWORK == 6
      BACKUP_DATA(10) = non_eq_abundances(ic2p)
#endif
      BACKUP_DATA(11) = fshield_H2
      BACKUP_DATA(12) = fshield_CO
      BACKUP_DATA(13) = AV_mean
      BACKUP_DATA(14) = chi_mean
      BACKUP_DATA(15) = dust_temp
#ifdef TR_ONTHESPOT
      BACKUP_DATA(16) = Lyc_uv_ion(1)
      BACKUP_DATA(17) = Lyc_HeatingRate(1)
#endif
      BACKUP_DATA(18) = G0RP
      BACKUP_DATA(19) = cr_rescale
!      BACKUP_DATA(20) = Lyc_HeatingRate

      ! print*, "time", time
      ! print*, "fshield_H2", fshield_H2
      ! print*, "fshield_CO", fshield_CO
      ! print*, "AV_mean", AV_mean
      ! print*, "chi_mean", chi_mean
      ! print*, "non_eq_abundances(ih2)", non_eq_abundances(ih2)
      ! print*, "non_eq_abundances(ihp)", non_eq_abundances(ihp)
      ! print*, "non_eq_abundances(ico)", non_eq_abundances(ico)
      ! print*, "dens", yn 
      ! print*, "energy", energy

c Check if the canary value istatus=101 for nothrow was provided.
c If yes, adhere to nstep_max parameter and avoid aborts at all cost
c Otherwise, assume the caller did use an old version of evolve_abundances,
c which does not set / care about those parameters
      if (istatus .eq. 4711) then
            nstep = nstep_max
            nothrow = .true.
      else
            nstep = 5000000
            nothrow = .false.
      endif

c
c Sanity check for our input time
c
      if (time .le. 0e0) then
#ifdef DEBUG_EVOLVE
        if (time .lt. -1e0) then
          print*, "Warning: negative timestep in evolve_abundances",
     $             time
        endif
#endif
        return
      endif
c
c Set redshift
c
      redshift = current_redshift
c
c Initialize abundance vector, treating the internal energy as
c an additional abundance.
c
c N.B. We assume that abundances are stored consecutively, with
c no gaps, and that itmp == nspec (i.e. that the energy is the last
c value stored)
c

      do I = 1, nspec-1
        y(I) = non_eq_abundances(I)
      enddo
      y(itmp)  = energy

c
c Setup initial parameters for DVODE integrator. These don't vary over
c the course of the timestep.
c
c Optional input supplied -- see iwork, rwork init. below:
c
      iopt = 1
c
c Evolve ODEs from t -> tout, allowing overshooting & interpolation
c [See DVODE documentation for more details].
c
      itask = 1
c
c Use same relative tolerance for each species, different absolute tolerances
c (since magnitudes may vary by many orders of magnitude).
c
      itol = 2
c
c Set relative tolerance -- for simplicity we keep this fixed at 0.1% for now
c
      rtol(1:nspec) = RTOL_FIX !!1e-4
c
c ODE integration method: 22 => BDF with an internally generated Jacobian
c
c      mf = -22
      mf = 22
c
c Start DVODE in start-up mode (which forces it to initialize all of its
c internal arrays). Subsequent calls within the loop below use istate = 2
c
      istate = 1
      have_set_dvode = 0
c
c Optional input -- increase maximum number of steps that DVODE can take
c internally without returning an to 20000, as 500 doesn't appear to be
c enough for portions of our parameter space
c
      do I = 5, 10
        rwork(I) = 0e0
        iwork(I) = 0
      enddo
c
      iwork(6) = 400000
c
      fshield_dust    = exp(-2.5e0 * AV_mean)
c      fshield_dust_H2 = exp(-3.5e0 * AV_mean)
c      fshield_dust_H2 = exp(-4.18e0 * AV_mean)  ! SCOG: from Heays+17
      fshield_dust_H2 = exp(-1.0*phgam(1)*AV_mean)

c
c Parameters to pass to rate_eq
c
      rpar(1) = yn
      rpar(2) = dl
      rpar(3) = divv
      rpar(4) = fshield_H2
      rpar(5) = fshield_CO
      rpar(6) = fshield_dust  ! Scaling used for CO, PE effect etc.
      rpar(7) = chi_mean
      rpar(8) = fshield_dust_H2
      rpar(9) = G0RP
      rpar(10) = cr_rescale
#ifdef TR_ONTHESPOT
c ** SV ** Add EUVE flux and energy flux to rpar
      do I = 1, nrrec
        rpar(nrpar+I) = Lyc_uv_ion(I)
        rpar(nrpar+nrrec+I) = Lyc_HeatingRate(I)
      enddo
      ipar(4) = tr_osIHII
      ipar(5) = tr_osIHeII
      ipar(6) = tr_osIHeIII
#endif
#ifdef XRAY_TR
c ** JM ** Add X-ray flux to list of parameters for DVODE
      do I = 1,Nbins !xdr_nbins
        rpar(nrpar+2*nrrec+I) = Fin(I)
      enddo
c *** JM *** Add Nbins to ipar so that rate_eq knows what it is.
      ipar(3) = Nbins
#else
      ipar(3) = 0
#endif
c      print*,"RPAR",rpar

!      print*,'evolve_abundances: rpar=',rpar,yn
c
      ipar(1) = 0
      ipar(7) = iccool
c
      t_start  = 0e0
      t_end    = time
      dt_correction = 1e0
      new_step = 1
c
c ===============================================================
      do I = 1, nstep
c
c If we did Euler integration on the last step, then set istate=1 so
c that DVODE will re-initialise arrays to the input values.
c
c        if (have_set_dvode .eq. 0) then
c          istate = 1
c        endif
c
c If we're starting a new step, then figure out a reasonable timestep to
c pass to DVODE -- we need to do this here to avoid DVODE producing -ve
c abundances or -ve temperatures, which can occur if it starts with too
c large a timestep.
c
c If we're not starting a new step, then dt has already been set.
c
        if (new_step .eq. 1) then
c
c First, save a copy of the current abundances
c
          do J = 1, nspec
            y_old(J) = y(J)
          enddo
          tdust_old = tdust
c
c Now compute ydot. On the first step, we must use rate_eq directly,
c but on subsequent steps, we can just call dvindy [if we're using
c the dvode integrator]
c
c          if (have_set_dvode .eq. 0) then !I .eq. 1) then
c            print*,'Calling rate_eq in step ',I
c            print*,nspec,t_start,y,ydot,rpar,ipar
            call rate_eq(nspec, t_start, y, ydot, rpar, ipar)
c          else
c            print*,'Calling dvindy in step ',I
c            print*,t_start,rwork(21),nspec,ydot,iflag_dvi
c            print*,'ydot ', ydot
c            call dvindy(t_start, 1, rwork(21), nspec, ydot, iflag_dvi)
c          endif
c
c Init. dt_est to an unrealistic value
c
          dt_est = 1e18

          do J = 1, nspec
            if (used_for_stepsize_control(J) .eq. 1) then
              if (ydot(J)   .ne. 0e0 .and.
     $            safety(J) .ne. 0e0 .and.
     $            y(J)      .gt. ATOL(J)) then
                dt_est = min(dt_est, abs(safety(J) * y(J) / ydot(J)))
              endif
            endif
          enddo
c
          ipar(2) = 0
c
c Don't evolve past the end of the global timestep
c
          deltat = min(dt_est * dt_correction, t_end - t_start)
          deltat_old=deltat
        endif
c
        t_in  = t_start
        t_out = t_start + deltat
        if (deltat .lt. 0e0) then
          print*, 'Negative timestep!'
          if (nothrow) then
            istatus = 1
            return
          else
            ABORT(1)
          endif !(nothrow)
        endif

        if (deltat .eq. 0e0) then
          if (t_end .eq. t_start) then
            return
          else
            print*, "times:", t_in, t_out, t_start, t_end, deltat
            print*, "params:", BACKUP_DATA
            print*, "dt_est:", dt_est
            print*, 'Zero timestep!'
            if (nothrow) then
              istatus = 2
              return
            else
              ABORT(2)
            endif !(nothrow)
          endif
        endif
c
ccc        if ( y(ico) .lt. 0 ) then
c          print*,'before dvode:'
c          print*,'nspec,y',nspec,y
ccc        endif
c         print*,'t_in,t_out=',t_in,t_out,t_start,deltat
c        print*,'dt_est=',dt_est, y,ydot
c        print*,'itol=',itol
c        print*,'atol=',ATOL
c        print*,'rtol=',rtol
c        print*,'itask,istate=',itask,istate
c ** JM **
c Check if anything is changing.  If not, then do a forward Euler
c integration.  If so, then call the big dvode function.
c If we do Euler integration, we set istate=2,ipar(1)=0 to indicate
c successful integration, so that new_state will be reset at the end
c of the loop.
c
c        ! Get max rate of change from y.
c        maxdy = 0.0
c        do J = 1, nspec
c          maxdy = max(maxdy, abs(ydot(J)/(y(J)+ATOL(J))))
c        enddo
c        ! If relative change over step is <1e-3, then just do Euler
c        if (maxdy*deltat .lt. 1.0e-3) then
c          do J = 1, nspec
c            y(J) = y(J) + ydot(J)*deltat
c          enddo
c          istate  = 2
c          ipar(1) = 0
c          have_set_dvode = 0
c        ! Otherwise call DVODE.
c        else
          call dvode(rate_eq, nspec, y, t_in, t_out, itol,
     $               itask, istate, iopt, rwork, lrw, iwork, liw, jac,
     $               mf, rpar, ipar)
c          have_set_dvode = 1
c        endif
c
c ** JM **

c        if (I .gt. 3000 .and. ipar(1) .eq. 0 .and. istate .eq. 2) then
c          print*, "i, nstep = ",i, nstep
c          print*, "dt", t_start, deltat, dt_est, t_end, t_out, t_in
c          !!print*, yn, dl, divv
c          print*, y_old
c          print*, ydot
c          print*, y
c          print*, " ******* continuing ....... ******** "
c        endif
c
c Check if the change rate calculated by dvode would have allowed for
c a larger timestep and adjust the next timestep accordingly.
c 
        if (turbostep .ne. 0) then
          dt_est2 = 1e18
c
          do J = 1, nspec
            if (used_for_stepsize_control(J) .eq. 1) then
              if (ydot(J)   .ne. 0e0 .and.
     $            safety(J) .ne. 0e0 .and.
     $            y_old(J)  .gt. ATOL(J)) then
                dt_est2 = min(dt_est2,
     $                        abs(safety(J) * y_old(J) *
     $                            deltat / (y(J) - y_old(J))))
              endif
            endif
          enddo
c
          if (turbostep .eq. 1) then
            dt_correction = sqrt(dt_est2 / dt_est) 
          else if (turbostep .eq. 2) then  
            dt_correction = dt_correction * sqrt(dt_est2 / deltat)
          else
            if (nothrow) then
              istatus = 1000
              return
            else
              ABORT(1000)
            endif
          endif  
        endif
c
c
c Ensure that dvode has produced sane output values -- if not, we
c need to shrink the timestep. Note that if ipar(1) is non-zero,
c or if istate != 2 (indicating an integrator error), we're going
c to shrink the timestep anyway, so we can skip the check.
c
        if (ipar(1) .eq. 0 .and. istate .eq. 2) then
          call check_abundances(y, ipar)
        endif
c
!******************
c        ipar(1)=0
c        istate =2
!******************

        if ((ipar(1) .ne. 0) .or. (istate .ne. 2)) then
#ifdef DEBUG_EVOLVE
         if (istate .ne. 2) then
           print*, "Error: istate = ", ipar,istate
         endif
         if (ipar(1) .ne. 0) then
           print*, "Error: ipar = ", ipar,istate
         endif
#endif
          do J = 1, nspec
            y(J) = y_old(J)
          enddo
          tdust    = tdust_old
          deltat       = deltat / 2e0
          dt_correction = 1e0
          if (deltat .lt. minDeltaT) then
            if (useAssistH2Decrease) then
c             We are in a high density case:
c             It might happen that H2 is recombining although it should
c             be ionised instead.
c             Check if this is the case
c             SOME DEBUG PRINTS
c              print*, "deltat is ", deltat
c              print*, y
c              print*, ydot

c             find species causing the minimum timestep
              counterMinDt = -1
              minDtBackup = 1e99
              currDtMin = 1e99
              do J = 1, nspec
              if (ydot(J)   .ne. 0e0 .and.
     $            safety(J) .ne. 0e0 .and.
     $            y(J)      .gt. ATOL(J)) then
                currDtMin = abs(safety(J) * y(J) / ydot(J))
                if (currDtMin .lt. minDtBackup) then
                  minDtBackup = currDtMin
                  if (ydot(J).lt.0) then
                    counterMinDt = -J
                  else
                    counterMinDt = J
                  endif
                endif
              endif
              enddo

              if (abs(counterMinDt).eq.ihp) then
c             ihp is causing trouble!!!
               if (counterMinDt.lt.0.0) then
c              ihp wants to recombine (ydot(ihp) is negative), force it
               y(ihp) = ATOL_HP
               else
c              we want to ionize (ydot(ihp) is positive)
               y(ih2) = ATOL_H2
               y(ihp) = 0.9
               goto 100 ! return
               endif
              endif
#ifdef TR_ONTHESPOT
              if (abs(counterMinDt).eq.ih2) then
               if ((ydot(ihp).gt.0.0) .and. (ydot(ih2).lt.0.0) .and.
     $            (y(ih2).lt.y(ihp)*1e-2) .and.
     $            (Lyc_HeatingRate(1).gt.0)) then
                 y(ih2) = ATOL_H2
                 y(ihp) = ATOL_HP
               endif
              endif
#endif
            else
              print*, "Excessively small timestep! ", deltat, dt_est,
     $                t_start, t_end, tdust, tdust_old
              print*, yn, dl, divv
              print*, y
              print*, ydot
              print*, ih2, ihp, ico
#ifdef TREE_RAD
              print*, "Tree rad totals: ", column_density_projection
#endif
#ifdef TREE_RAD_H2
              print*, "Tree rad h2: ", column_density_projection_h2
#endif
#ifdef TREE_RAD_CO
              print*, "Tree rad co: ", column_density_projection_co
#endif
              print*, "Try setting ch_useAssistH2Decrease = true."
              print*, "More comments in Chemistry_init."
              print*, "Good luck!"

              if (nothrow) then
                istatus = 100
                return
              else
                ABORT(100)
              endif !(nothrow)

            endif
          endif
          new_step = 0
          ipar(1)  = 0
          istate   = 1
        else
          t_start = t_start + deltat
          new_step = 1
        endif
c

c        print*,'t_start, dt, t_end', t_start,deltat,t_end

        if (t_start .ge. t_end) then
          goto 100
        endif

        if ((I+10).ge.nstep) then
          print*, "evolve_abundances: Pre-failure (nstep) data dump"
          print*, "step, nstep_max:", I, nstep
          print*, "deltat, dt_est, dt_est2:", deltat, dt_est, dt_est2
          print*, "dt_correction:", dt_correction
          print*, "t_start, t_end:", t_start, t_end
          print*, "yn, dl, divv:", yn, dl, divv
          print*, "y:", y
          print*, "ydot:", ydot
          print*, "tdust, tdust_old:", tdust, tdust_old
          print*, "ih2, ihp, ico:", ih2, ihp, ico
          print*, "evolve_abundances: End of data dump"
        endif

c         
      enddo

c End of nstep loop
c ===============================================================

      write(*,'(A)') 'EVOLVE_ABUNDANCES: nstep exceeded,'  
      print*, deltat
      print*, rpar
      print*, y
      print*, ydot
      print*, "Try setting ch_turbostep = 2"
      if (nothrow) then
        istatus = 10
        return
      else
        ABORT(0)
      endif !(nothrow)
      return

100   continue

      do I = 1, nspec-1
        non_eq_abundances(I) = y(I)
      enddo
      energy = y(itmp)

c Calculate final dust temp.
#if CHEMISTRYNETWORK == 4 || CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
      abHI  = max(0e0, 1e0 - 2e0 * y(ih2) - y(ihp))
      abe   = y(ihp) ! + abhep ! this appears to be a bug since abhep is never set DS + SG 2014
#endif
#if CHEMISTRYNETWORK == 15
      abe   = max(y(ihp) + y(icp) + y(ihcop) + y(ihep) + y(imp)
     $           , 1e-40)
      abHI  = max(1d0 - 2. * y(ih2) - y(ihp) - y(ichx) - y(iohx)
     $            - y(ihcop), 1e-40)
#endif
#if CHEMISTRYNETWORK == 17
      abe   = max(y(ihp) + y(icp) + y(ihcop) + y(ihep) + y(imp)
     $            + 2. * y(ic2p) + 2. * y(ihe2p) , 1e-40)
      abHI  = max(1d0 - 2. * y(ih2) - y(ihp) - y(ichx) - y(iohx)
     $            - y(ihcop), 1e-40)
#endif
      yntot = (1e0 + abhe - y(ih2) + abe) * yn
      temp  = 2e0 * energy / (3e0 * yntot * kboltz)
#ifdef XRAY_TR
      call calc_dust_temp(yn, chi_mean, temp, abHI, RH2,Rdust,HX,G0RP,
     $                    dl, 0)
#else
      call calc_dust_temp(yn, chi_mean, temp, abHI, RH2, Rdust, G0RP,
     $                    dl, 0)
#endif
      dust_temp = tdust
c

c     Everything went right: Set status code 999
      istatus = 999

      ! print*, "ih2", non_eq_abundances(ih2)
      ! print*, "ihp", non_eq_abundances(ihp)
      ! print*, "ico", non_eq_abundances(ico)
      ! print*, "energy", energy
      ! print*, "temp", temp
      ! print*, "G0RP", G0RP

      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////     E V O L V E _ A B U N D A N C E S     \\\\\\\\\\
c
c=======================================================================
c
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c
      subroutine check_abundances(y, ipar)
      implicit NONE
#include "param.h"
#include "cool.h"
#include "non_eq.h"
      REAL y(nspec), max_abundance(nspec)
      REAL current_carbon, current_oxygen
      REAL diff
      integer ipar(nipar)
      integer J
c
      max_abundance = 1.0e30
c
#if CHEMISTRYNETWORK == 4
        max_abundance(ih2)   = 0.5e0
        max_abundance(ihp)   = 1.0e0
#endif
c
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
        max_abundance(ih2)   = 0.5e0
        max_abundance(ihp)   = 1.0e0
        max_abundance(ico)   = min(abundc, abundo)
#endif
c
#if CHEMISTRYNETWORK == 6
        max_abundance(ic2p)  = abundc
#endif
c
#if (CHEMISTRYNETWORK == 15) || (CHEMISTRYNETWORK == 16) || (CHEMISTRYNETWORK == 17)
        max_abundance(ih2) = 0.5d0
        max_abundance(ihp) = 1.0d0
        max_abundance(icp)  = abundc
        max_abundance(ichx) = abundc
        max_abundance(iohx) = abundo
        max_abundance(ico)  = min(abundo, abundc)
        max_abundance(ihcop) = min(abundo, abundc)
        max_abundance(ihep) = abhe
        max_abundance(imp)  = abundM
#endif
#if CHEMISTRYNETWORK == 17
        max_abundance(ic2p) = abundc
        max_abundance(ihe2p) = abhe
#endif
c
c Check for NaNs
c
      do J = 1, nspec
        if (y(J) .ne. y(J)) then
          ipar(1) = 1
#ifdef DEBUG_EVOLVE
          print*, "NaN abundance", y
#endif
          return
        endif
      enddo
c
c Check for negative abundances. Small negative values (which have
c absolute values less than the corresponding absolute tolerance)
c are accepted, but are clipped to zero. For large negative values,
c we repeat the step with a smaller dt.
c
      do J = 1, nspec
        if (y(J) .lt. 0e0) then
          if (y(J) .lt. -2.0*ATOL(J)) then
            ipar(1) = 1
#ifdef DEBUG_EVOLVE
            print*, "check_abundances :: negative abundance", y
#endif
            return
          else
            y(J) = 1e-40
          endif
        endif
      enddo
c
c Don't allow ion abundances to exceed corresponding total metal
c abundances. Again, small errors are silently corrected, large
c errors lead to us repeating with a smaller timestep
c
      do J = 1, nspec
        if (J .ne. itmp .and. y(J) .gt. max_abundance(J)) then
          diff = y(J) - max_abundance(J)
          if (diff .lt. 1e-4 * y(J)) then
            y(J) = max_abundance(J)
          else
#ifdef DEBUG_EVOLVE
            print*, "Excess abundance", y ,J,y(J)
#endif
            ipar(1) = 1
            return
          endif
        endif
      enddo
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D   S U B R O U T I N E      //////////
c    //////////      C H E C K _ A B U N D A N C E S    \\\\\\\\\\
c
c=======================================================================
c=======================================================================
c
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================
c
      subroutine init_tolerances
      implicit NONE
#include "param.h"
#include "non_eq.h"
#include "cool.h"
c
c Set-up absolute tolerances
c
#if CHEMISTRYNETWORK == 4
        ATOL(ihp) = ATOL_HP
        ATOL(ih2) = ATOL_H2
#endif
c
#if CHEMISTRYNETWORK == 5 || CHEMISTRYNETWORK == 6
        ATOL(ihp) = ATOL_HP
        ATOL(ih2) = ATOL_H2
        ATOL(ico) = ATOL_CO
#endif
c
#if CHEMISTRYNETWORK == 6
        ATOL(ic2p)= ATOL_C2P
#endif
c
#if (CHEMISTRYNETWORK == 15) || (CHEMISTRYNETWORK == 16) || (CHEMISTRYNETWORK == 17)
        atol(ih2)   = ATOL_H2
        atol(ihp)   = ATOL_HP
        atol(icp)   = ATOL_CP
        atol(ichx)  = ATOL_CHX
        atol(iohx)  = ATOL_OHX
        atol(ico)   = ATOL_CO
        atol(ihcop) = ATOL_HCOP
        atol(ihep)  = ATOL_HEP
        atol(imp)   = ATOL_MP
#endif

#if CHEMISTRYNETWORK == 17
        atol(ic2p)  = ATOL_C2P
        atol(ihe2p)  = ATOL_HE2P
#endif
c
      ATOL(itmp)  = ATOL_TMP
c
      return
      end
c=======================================================================
c
c
c    \\\\\\\\\\        E N D     S U B R O U T I N E      //////////
c    //////////       I N I T _ T O L E R A N C E S       \\\\\\\\\\
c
c=======================================================================

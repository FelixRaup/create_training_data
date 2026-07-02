#include "Flash.h"


c=======================================================================
c
c    \\\\\\\\\\      B E G I N   S U B R O U T I N E      //////////
c    //////////            G E T _ K A P P A              \\\\\\\\\\
c
c=======================================================================
c
      subroutine tr_rpGetKappa(temp, kappa_R, kappa_P)
      
      implicit none
c
      REAL temp, kappa_R, kappa_P
      REAL tr_rpKappa
      LOGICAL tr_rpTsquare
c      
      tr_rpKappa = 1e3
      tr_rpTsquare = .true.
c
      if (tr_rpTsquare) then
        if (temp .lt. 150.0) then
          kappa_P = 0.001 * temp**2 
          kappa_R = 0.0003162277660 * temp**2
        elseif (temp .lt. 1500) then
          kappa_P = 0.001 * 150.0**2
          kappa_R = 0.0003162277660 * 150.0**2
        else
          kappa_P = 0.0
          kappa_R = 0.0
        endif
      else
          kappa_P = tr_rpKappa 
          kappa_R = tr_rpKappa 
      endif
c
      return
      end
      
c=======================================================================
c
c
c    \\\\\\\\\\        E N D    S U B R O U T I N E       //////////
c    //////////            G E T _ K A P P A              \\\\\\\\\\
c
c=======================================================================


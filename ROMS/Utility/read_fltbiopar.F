#include "cppdefs.h"
#if defined FLOATS && defined FLOAT_BIOLOGY
!!
!!svn $Id: read_fltbiopar.F 709 2014-01-23 20:09:38Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports biological floats model input       !
!!  parameters.                                                        !
!!                                                                     !
!!======================================================================
!!
# if defined FLOAT_OYSTER
#  include <oyster_floats_inp.h>
# endif
#else
      SUBROUTINE read_FltBioPar
      END SUBROUTINE read_FltBioPar
#endif

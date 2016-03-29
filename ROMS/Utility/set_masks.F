#include "cppdefs.h"
      MODULE set_masks_mod
#if defined MASKING || defined TIDAL_TURBINE 
!
!svn $Id: set_masks.F 709 2014-01-23 20:09:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!                                                                      !
!  With Portions Copyright (c) 2014 Thomas Roc and ITPower             !
!    Licensed under an Affero GPL style license                        !
!    See License_ROMS-Tidal.txt                                        !
!=======================================================================
!                                                                      !
!  These routines set internal Land/Sea masking arrays that are used   !
!  to process fields into output NetCDF files.  The Land grid points   !
!  are replaced by the _FillValue in the output files to  facilitate   !
!  post-processing with generic tools.                                 !
!                                                                      !
!  If point sources, insure that masks at point source locations are   !
!  set to water to avoid masking with _FillValue at those locations.   !
# ifdef WET_DRY
!
!  If wetting and drying, masks array are time dependent and changed   !
!  at every time-step in routine "wetdry".  Notice that time-average   !
!  masks are needed for the selected time window.
# endif
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC :: set_masks
# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
      PUBLIC :: set_avg_masks
# endif
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE set_masks (ng, tile, model)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
# if defined TIDAL_TURBINE
      USE mod_scalars
# endif
!========================================== End of the portion =========
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
# include "tile.h"
!
# ifdef PROFILE
      CALL wclock_on (ng, model, 2)
# endif
      CALL set_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
# if !defined MASKING
!================================================== Thomas Roc =========
#  if defined TIDAL_TURBINE
     &                     GRID(ng) % om_r,                             &
     &                     GRID(ng) % on_r,                             &
     &                     GRID(ng) % om_u,                             &
     &                     GRID(ng) % on_v,                             &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % h,                                &
#   ifdef TIDAL_TURBINE
     &                     GRID(ng) % imask_t,                          &
     &                     GRID(ng) % rmask_t,                          &
     &                     GRID(ng) % rmask1t,                          &
     &                     GRID(ng) % rmask2t,                          &
     &                     GRID(ng) % rmask3t,                          &
     &                     GRID(ng) % umask_t,                          &
     &                     GRID(ng) % vmask_t,                          &
#   endif
     &                     GRID(ng) % xr,                               &
     &                     GRID(ng) % yr,                               &
     &                     GRID(ng) % xu,                               &
     &                     GRID(ng) % yu,                               &
     &                     GRID(ng) % xv,                               &
     &                     GRID(ng) % yv,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w)
#  endif
# endif
!========================================== End of the portion =========
# ifdef MASKING
     &                     GRID(ng) % pmask,                            &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE
     &                     GRID(ng) % om_r,                             &
     &                     GRID(ng) % on_r,                             &
     &                     GRID(ng) % om_u,                             &
     &                     GRID(ng) % on_v,                             &
     &                     GRID(ng) % Hz,                               &
     &                     GRID(ng) % h,                                &
     &                     GRID(ng) % imask_t,                          &
     &                     GRID(ng) % rmask_t,                          &
     &                     GRID(ng) % rmask1t,                          &
     &                     GRID(ng) % rmask2t,                          &
     &                     GRID(ng) % rmask3t,                          &
     &                     GRID(ng) % umask_t,                          &
     &                     GRID(ng) % vmask_t,                          &
     &                     GRID(ng) % xr,                               &
     &                     GRID(ng) % yr,                               &
     &                     GRID(ng) % xu,                               &
     &                     GRID(ng) % yu,                               &
     &                     GRID(ng) % xv,                               &
     &                     GRID(ng) % yv,                               &
     &                     GRID(ng) % z_r,                              &
     &                     GRID(ng) % z_w,                              &
#  endif
!========================================== End of the portion =========
#  if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                     GRID(ng) % pmask_avg,                        &
     &                     GRID(ng) % rmask_avg,                        &
     &                     GRID(ng) % umask_avg,                        &
     &                     GRID(ng) % vmask_avg,                        &
#  endif
#  ifdef DIAGNOSTICS
     &                     GRID(ng) % pmask_dia,                        &
     &                     GRID(ng) % rmask_dia,                        &
     &                     GRID(ng) % umask_dia,                        &
     &                     GRID(ng) % vmask_dia,                        &
#  endif
     &                     GRID(ng) % pmask_io,                         &
     &                     GRID(ng) % rmask_io,                         &
     &                     GRID(ng) % umask_io,                         &
     &                     GRID(ng) % vmask_io)
# endif
# ifdef PROFILE
      CALL wclock_off (ng, model, 2)
# endif

      RETURN
      END SUBROUTINE set_masks
!
!***********************************************************************
      SUBROUTINE set_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
# if !defined MASKING
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE 
     &                           om_r, on_r, om_u, on_v, Hz, h,         &
     &                           imask_t, rmask_t, rmask1t, rmask2t,    &
     &                           rmask3t, umask_t, vmask_t,             &
     &                           xr, yr, xu, yu, xv, yv, z_r, z_w)             
#  endif
# endif
# ifdef MASKING 
     &                           pmask, rmask,                          &
     &                           umask, vmask,                          &
#  if defined TIDAL_TURBINE
     &                           om_r, on_r, om_u, on_v, Hz, h,         &
     &                           imask_t, rmask_t, rmask1t, rmask2t,    &
     &                           rmask3t, umask_t, vmask_t,             &
     &                           xr, yr, xu, yu, xv, yv, z_r, z_w,      &
#  endif
!========================================== End of the portion =========
#  if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
     &                           pmask_avg, rmask_avg,                  &
     &                           umask_avg, vmask_avg,                  &
#  endif
#  ifdef DIAGNOSTICS
     &                           pmask_dia, rmask_dia,                  &
     &                           umask_dia, vmask_dia,                  &
#  endif
     &                           pmask_io, rmask_io,                    &
     &                           umask_io, vmask_io)
# endif
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
# if defined TIDAL_TURBINE
      USE exchange_3d_mod
# endif
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
      USE distribute_mod, ONLY : mp_collect
# endif
!========================================== End of the portion =========
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
# ifdef ASSUMED_SHAPE
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#  endif
!========================================== End of the portion =========
#  ifdef MASKING
      real(r8), intent(in) :: pmask(LBi:,LBj:)
      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#   endif
#   ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: umask_dia(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_dia(LBi:,LBj:)
#   endif
      real(r8), intent(inout) :: pmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_io(LBi:,LBj:)
      real(r8), intent(inout) :: umask_io(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_io(LBi:,LBj:)
#  endif
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE
      real(r8), intent(out) :: imask_t(LBi:,LBj:,:)
      real(r8), intent(out) :: rmask_t(LBi:,LBj:,:)
      real(r8), intent(out) :: rmask1t(LBi:,LBj:,:)
      real(r8), intent(out) :: rmask2t(LBi:,LBj:,:)
      real(r8), intent(out) :: rmask3t(LBi:,LBj:,:)
      real(r8), intent(out) :: umask_t(LBi:,LBj:,:)
      real(r8), intent(out) :: vmask_t(LBi:,LBj:,:)       
      real(r8), intent(in) :: xr(LBi:,LBj:)       
      real(r8), intent(in) :: yr(LBi:,LBj:)       
      real(r8), intent(in) :: xu(LBi:,LBj:)       
      real(r8), intent(in) :: yu(LBi:,LBj:)       
      real(r8), intent(in) :: xv(LBi:,LBj:)       
      real(r8), intent(in) :: yv(LBi:,LBj:)       
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)       
      real(r8), intent(in) :: z_w(LBi:,LBj:,:)       
#  endif
!========================================== End of the portion =========
# else
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE
      real(r8), intent(in) :: om_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_r(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: om_u(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: on_v(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#  endif
!========================================== End of the portion =========
#  ifdef MASKING 
      real(r8), intent(in) :: pmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   if defined AVERAGES    || \
     (defined AD_AVERAGES && defined ADJOINT) || \
     (defined RP_AVERAGES && defined TL_IOMS) || \
     (defined TL_AVERAGES && defined TANGENT)
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#   endif
#   ifdef DIAGNOSTICS
      real(r8), intent(inout) :: pmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_dia(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_dia(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(inout) :: pmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_io(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_io(LBi:UBi,LBj:UBj)
#  endif
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
#  if defined TIDAL_TURBINE
      real(r8), intent(inout) :: imask_t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rmask_t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rmask1t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rmask2t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: rmask3t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: umask_t(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: vmask_t(LBi:UBi,LBj:UBj,N(ng))       
      real(r8), intent(in) :: xr(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: yr(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: xu(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: yu(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: xv(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: yv(LBi:UBi,LBj:UBj)       
      real(r8), intent(in) :: z_r(LBi:UBi,LBj:UBj,N(ng))       
      real(r8), intent(in) :: z_w(LBi:UBi,LBj:UBj,N(ng))       
#  endif
!========================================== End of the portion ========= 
# endif
!
!  Local variable declarations.
!
      integer :: i, is, j

!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
# if defined TIDAL_TURBINE
      integer :: ix, iy, iz, nx, ny, nz, k
      real(r8) :: xpos, ypos, zpos, rposx, rposy
      real(r8) :: xposmin, yposmin, zposmin, rposxmin, rposymin
      real(r8) :: xposmax, yposmax, zposmax, rposxmax, rposymax
      real(r8) :: test, ct, cd, ctke, cgls, lcpa, diam
      integer, dimension(Nturbines(ng)) :: Itrb, Jtrb, Ktrb
# endif
!========================================== End of the portion =========

# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize internal history files Land/Sea masks with its respective
!  application grid mask.
!-----------------------------------------------------------------------
!
# ifdef MASKING
      DO j=JstrP,JendP
        DO i=IstrP,IendP
          pmask_io(i,j)=pmask(i,j)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rmask_io(i,j)=rmask(i,j)
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          umask_io(i,j)=umask(i,j)
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vmask_io(i,j)=vmask(i,j)
        END DO
      END DO
!
!  Insure that masks at mass point source locations are set to water
!  to avoid masking with _FillValue at those locations.
!
      IF (LuvSrc(ng)) THEN
        DO is=1,Nsrc(ng)
          i=SOURCES(ng)%Isrc(is)
          j=SOURCES(ng)%Jsrc(is)
          IF (((IstrT.le.i).and.(i.le.IendT)).and.                        &
     &        ((JstrT.le.j).and.(j.le.JendT))) THEN
            IF (INT(SOURCES(ng)%Dsrc(is)).eq.0) THEN
              umask_io(i,j)=1.0_r8
            ELSE
              vmask_io(i,j)=1.0_r8
            END IF
          END IF
        END DO
      END IF
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_io)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_io)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_io)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_io)
      END IF

#  ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_io, rmask_io, umask_io, vmask_io)
#  endif
# endif

!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
# if defined TIDAL_TURBINE
        DO is=1,Nturbines(ng)         
          IF (SCALARS(ng)%Tc(is).gt.0) THEN
#  ifdef TIDAL_TURBINE
! Define turbines' positions and features
              xpos=SCALARS(ng)%Txpos(is)
              ypos=SCALARS(ng)%Typos(is)
              diam=SCALARS(ng)%Tdiam(is)
              ct=SCALARS(ng)%Tct(is)!Thrust coefficient
              ctke=SCALARS(ng)%Ttke(is)!Tke source coefficient
              cgls=SCALARS(ng)%Tgls(is)!Gls second-order correction-term coefficient
              lcpa=SCALARS(ng)%Tlc(is)*SIN(SCALARS(ng)%Tpa(is))!Blade apparent width
              cd=4.0_r8*((1.0_r8-((1.0_r8-ct)**0.5_r8))/                &
     &                   (1.0_r8+((1.0_r8-ct)**0.5_r8)))*lcpa!Drag coeff * lcpa
! Create turbulence term mask
#   ifdef I_ORIENTATION
            DO iz=1,N(ng) 
              DO iy=JstrR,Jend
                DO ix=IstrR,Iend
              zpos=-1.0_r8*SCALARS(ng)%Tzpos(is)*h(ix,iy)
!        WRITE(*,*) 'Turb.',is, 'zpos', zpos
              rposxmin=xpos-(0.5_r8*om_r(ix,iy))
              rposxmax=xpos+(0.5_r8*om_r(ix,iy))
            IF ((xr(ix,iy).gt.rposxmin).AND.(xr(ix,iy).lt.rposxmax)) THEN
               test=(((zpos-z_w(ix,iy,iz))**2.0_r8)+                    &
     &               ((ypos-yr(ix,iy))**2.0_r8))**0.5_r8 
             IF (test.lt.(diam/2.0_r8)) THEN
                           rmask_t(ix,iy,iz)=ct
                           rmask1t(ix,iy,iz)=ctke
                           rmask2t(ix,iy,iz)=cgls
                           rmask3t(ix,iy,iz)=lcpa
              END IF
            END IF
                END DO
              END DO
            END DO
#   endif
#   ifdef J_ORIENTATION
            DO iz=1,N(ng) 
              DO iy=JstrR,Jend
                DO ix=IstrR,Iend
              zpos=-1.0_r8*SCALARS(ng)%Tzpos(is)*h(ix,iy)
!        WRITE(*,*) 'Turb.',is, 'zpos', zpos
              rposymin=ypos-(0.5_r8*on_r(ix,iy))
              rposymax=ypos+(0.5_r8*on_r(ix,iy))
            IF ((yr(ix,iy).gt.rposymin).AND.(yr(ix,iy).lt.rposymax)) THEN
               test=(((zpos-z_w(ix,iy,iz))**2.0_r8)+                    &
     &               ((xpos-xr(ix,iy))**2.0_r8))**0.5_r8 
             IF (test.lt.(diam/2.0_r8)) THEN
                           rmask_t(ix,iy,iz)=ct
                           rmask1t(ix,iy,iz)=ctke
                           rmask2t(ix,iy,iz)=cgls
                           rmask3t(ix,iy,iz)=lcpa
              END IF
            END IF
                END DO
              END DO
            END DO
#   endif
!
! Create u-mask
!
#   ifdef I_ORIENTATION
            DO iz=1,N(ng) 
              DO iy=JstrR,Jend
                DO ix=IstrR,Iend
              zpos=-1.0_r8*SCALARS(ng)%Tzpos(is)*h(ix,iy)
!        WRITE(*,*) 'Turb.',is, 'zpos', zpos
              xposmin=xpos-(1.0_r8*om_u(ix,iy))
              xposmax=xpos+(1.0_r8*om_u(ix,iy))
            IF ((xu(ix,iy).gt.xposmin).AND.(xu(ix,iy).lt.xpos)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                &
     &               ((ypos-yu(ix,iy))**2.0_r8))**0.5_r8
             IF (test.lt.(diam/2.0_r8)) THEN
                           umask_t(ix,iy,iz)=cd
                           imask_t(ix,iy,iz)=1.0_r8
!        WRITE(*,*) 'Mask u, turb',is, 'i', ix, 'j', iy, 'k', iz
             END IF
            END IF
            IF ((xu(ix,iy).gt.xpos).AND.(xu(ix,iy).lt.xposmax)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                &
     &               ((ypos-yu(ix,iy))**2.0_r8))**0.5_r8
             IF (test.lt.(diam/2.0_r8)) THEN
                           umask_t(ix,iy,iz)=cd
                           imask_t(ix,iy,iz)=1.0_r8
!        WRITE(*,*) 'Mask u, turb',is, 'i', ix, 'j', iy, 'k', iz
             END IF
            END IF
                END DO
              END DO
            END DO
#   endif
#   ifdef J_ORIENTATION
            DO iz=1,N(ng) 
              DO iy=JstrR,Jend
                DO ix=IstrR,Iend
              zpos=-1.0_r8*SCALARS(ng)%Tzpos(is)*h(ix,iy)
!        WRITE(*,*) 'Turb.',is, 'zpos', zpos
              yposmin=ypos-(1.0_r8*on_v(ix,iy))
              yposmax=ypos+(1.0_r8*on_v(ix,iy))
            IF ((yv(ix,iy).gt.yposmin).AND.(yv(ix,iy).lt.ypos)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                &
     &               ((xpos-xv(ix,iy))**2.0_r8))**0.5_r8
             IF (test.lt.(diam/2.0_r8)) THEN
                           vmask_t(ix,iy,iz)=cd
                           imask_t(ix,iy,iz)=1.0_r8
!        WRITE(*,*) 'Mask v, turb',is, 'i', ix, 'j', iy, 'k', iz
             END IF
            END IF
            IF ((yv(ix,iy).gt.ypos).AND.(yv(ix,iy).lt.yposmax)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                &
     &               ((xpos-xv(ix,iy))**2.0_r8))**0.5_r8
             IF (test.lt.(diam/2.0_r8)) THEN
                           vmask_t(ix,iy,iz)=cd
                           imask_t(ix,iy,iz)=1.0_r8
!        WRITE(*,*) 'Mask v, turb',is, 'i', ix, 'j', iy, 'k', iz
             END IF
            END IF
                END DO
              END DO
            END DO
#   endif
#  endif
!
          END IF
        END DO

      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
#  ifdef TIDAL_TURBINE
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          imask_t)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rmask_t)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rmask1t)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rmask2t)
        CALL exchange_r3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          rmask3t)
        CALL exchange_u3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          umask_t)
        CALL exchange_v3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 1, N(ng),           &
     &                          vmask_t)
#  endif
      END IF

#  ifdef DISTRIBUTE
#   ifdef TIDAL_TURBINE
      CALL mp_exchange3d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    imask_t, rmask_t, umask_t, vmask_t)
      CALL mp_exchange3d (ng, tile, model, 3,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    rmask1t, rmask2t, rmask3t)
#   endif
#  endif
# endif
!========================================== End of the portion =========
# ifdef MASKING
#  if defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!-----------------------------------------------------------------------
!  Initialize average file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
#   ifdef WET_DRY
          pmask_avg(i,j)=0.0_r8
#   else
          pmask_avg(i,j)=pmask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
#   ifdef WET_DRY
          rmask_avg(i,j)=0.0_r8
#   else
          rmask_avg(i,j)=rmask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrP,IendT
#   ifdef WET_DRY
          umask_avg(i,j)=0.0_r8
#   else
          umask_avg(i,j)=umask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
#   ifdef WET_DRY
          vmask_avg(i,j)=0.0_r8
#   else
          vmask_avg(i,j)=vmask_io(i,j)
#   endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_avg)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_avg)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_avg)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_avg)
      END IF

#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_avg, rmask_avg, umask_avg, vmask_avg)
#   endif
#  endif

#  ifdef DIAGNOSTICS
!
!-----------------------------------------------------------------------
!  Initialize diagnostic file Land/Sea masks for time-averaged fields.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendP
        DO i=IstrP,IendP
#   ifdef WET_DRY
          pmask_dia(i,j)=0.0_r8
#   else
          pmask_dia(i,j)=pmask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrT,IendT
#   ifdef WET_DRY
          rmask_dia(i,j)=0.0_r8
#   else
          rmask_dia(i,j)=rmask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrT,JendT
        DO i=IstrP,IendT
#   ifdef WET_DRY
          umask_dia(i,j)=0.0_r8
#   else
          umask_dia(i,j)=umask_io(i,j)
#   endif
        END DO
      END DO

      DO j=JstrP,JendT
        DO i=IstrT,IendT
#   ifdef WET_DRY
          vmask_dia(i,j)=0.0_r8
#   else
          vmask_dia(i,j)=vmask_io(i,j)
#   endif
        END DO
      END DO
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_p2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pmask_dia)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          rmask_dia)
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          umask_dia)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vmask_dia)
      END IF

#   ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pmask_dia, rmask_dia, umask_dia, vmask_dia)
#   endif
#  endif

# endif

      RETURN
      END SUBROUTINE set_masks_tile

# if defined WET_DRY     && \
     defined AVERAGES    || \
    (defined AD_AVERAGES && defined ADJOINT) || \
    (defined RP_AVERAGES && defined TL_IOMS) || \
    (defined TL_AVERAGES && defined TANGENT)
!
!***********************************************************************
      SUBROUTINE set_avg_masks (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          pmask_avg, rmask_avg,                   &
     &                          umask_avg, vmask_avg)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
#  ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#  ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: pmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: rmask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: umask_avg(LBi:,LBj:)
      real(r8), intent(inout) :: vmask_avg(LBi:,LBj:)
#  else
      real(r8), intent(inout) :: pmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: rmask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: umask_avg(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: vmask_avg(LBi:UBi,LBj:UBj)
#  endif
!
!
!  Local variable declarations.
!
      integer :: i, j

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Return if time-averaging window is zero.
!-----------------------------------------------------------------------
!
      IF (nAVG(ng).eq.0) RETURN
!
!-----------------------------------------------------------------------
!  If last time-step of average window, convert time dependent counters
!  for wet points to time-averaged Land/Sea masks (dry=0, wet=1) for
!  the current average window period. Notice that a grid point is wet
!  if the count is greater than zero for the current time average
!  window.
!-----------------------------------------------------------------------
!
      IF ((iic(ng).gt.ntsAVG(ng)).and.                                  &
     &    (MOD(iic(ng)-1,nAVG(ng)).eq.0).and.                           &
     &    ((iic(ng).ne.ntstart(ng)).or.(nrrec(ng).eq.0))) THEN

        DO j=JstrP,JendP
          DO i=IstrP,IendP
            pmask_avg(i,j)=MIN(1.0_r8, pmask_avg(i,j))
          END DO
        END DO
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            rmask_avg(i,j)=MIN(1.0_r8, rmask_avg(i,j))
          END DO
        END DO
        DO j=JstrT,JendT
          DO i=IstrP,IendT
            umask_avg(i,j)=MIN(1.0_r8, umask_avg(i,j))
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            vmask_avg(i,j)=MIN(1.0_r8, vmask_avg(i,j))
          END DO
        END DO
!
        IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
          CALL exchange_p2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            pmask_avg)
          CALL exchange_r2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            rmask_avg)
          CALL exchange_u2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            umask_avg)
          CALL exchange_v2d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vmask_avg)
        END IF

#  ifdef DISTRIBUTE
        CALL mp_exchange2d (ng, tile, model, 4,                         &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      NghostPoints,                               &
     &                      EWperiodic(ng), NSperiodic(ng),             &
     &                      pmask_avg, rmask_avg, umask_avg, vmask_avg)
#  endif

      END IF

      RETURN

      END SUBROUTINE set_avg_masks
# endif
#endif
      END MODULE set_masks_mod

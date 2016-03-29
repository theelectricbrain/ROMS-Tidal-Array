      MODULE set_masks_mod
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
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC :: set_masks
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
      USE mod_scalars
!========================================== End of the portion =========
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL wclock_on (ng, model, 2)
      CALL set_masks_tile (ng, tile, model,                             &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
!================================================== Thomas Roc =========
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
     &                     GRID(ng) % z_w)
!========================================== End of the portion =========
      CALL wclock_off (ng, model, 2)
      RETURN
      END SUBROUTINE set_masks
!
!***********************************************************************
      SUBROUTINE set_masks_tile (ng, tile, model,                       &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                           om_r, on_r, om_u, on_v, Hz, h,         &
     &                           imask_t, rmask_t, rmask1t, rmask2t,    &
     &                           rmask3t, umask_t, vmask_t,             &
     &                           xr, yr, xu, yu, xv, yv, z_r, z_w)             
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
      USE mod_sources
!
      USE exchange_2d_mod
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      USE exchange_3d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d, mp_exchange3d
      USE distribute_mod, ONLY : mp_collect
!========================================== End of the portion =========
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: h(LBi:,LBj:)
!========================================== End of the portion =========
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
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
!========================================== End of the portion =========
!
!  Local variable declarations.
!
      integer :: i, is, j
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      integer :: ix, iy, iz, nx, ny, nz, k
      real(r8) :: xpos, ypos, zpos, rposx, rposy
      real(r8) :: xposmin, yposmin, zposmin, rposxmin, rposymin
      real(r8) :: xposmax, yposmax, zposmax, rposxmax, rposymax
      real(r8) :: test, ct, cd, ctke, cgls, lcpa, diam
      integer, dimension(Nturbines(ng)) :: Itrb, Jtrb, Ktrb
!========================================== End of the portion =========
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Initialize internal history files Land/Sea masks with its respective
!  application grid mask.
!-----------------------------------------------------------------------
!
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
        DO is=1,Nturbines(ng)         
          IF (SCALARS(ng)%Tc(is).gt.0) THEN
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
!
! Create u-mask
!
            DO iz=1,N(ng) 
              DO iy=JstrR,Jend
                DO ix=IstrR,Iend
              zpos=-1.0_r8*SCALARS(ng)%Tzpos(is)*h(ix,iy)
!        WRITE(*,*) 'Turb.',is, 'zpos', zpos
              xposmin=xpos-(1.0_r8*om_u(ix,iy))
              xposmax=xpos+(1.0_r8*om_u(ix,iy))
            IF ((xu(ix,iy).gt.xposmin).AND.(xu(ix,iy).lt.xpos)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                    &
     &               ((ypos-yu(ix,iy))**2.0_r8))**0.5_r8
             IF (test.lt.(diam/2.0_r8)) THEN
                           umask_t(ix,iy,iz)=cd
                           imask_t(ix,iy,iz)=1.0_r8
!        WRITE(*,*) 'Mask u, turb',is, 'i', ix, 'j', iy, 'k', iz
             END IF
            END IF
            IF ((xu(ix,iy).gt.xpos).AND.(xu(ix,iy).lt.xposmax)) THEN
               test=(((zpos-z_r(ix,iy,iz))**2.0_r8)+                    &
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
!
          END IF
        END DO
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
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
      END IF
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
!========================================== End of the portion =========
      RETURN
      END SUBROUTINE set_masks_tile
      END MODULE set_masks_mod

      MODULE gls_corstep_mod
!
!svn $Id: gls_corstep.F 709 2014-01-23 20:09:38Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!                                                                      !
!  With Portions Copyright (c) 2014 Thomas Roc and ITPower             !
!    Licensed under an Affero GPL style license                        !
!    See License_ROMS-Tidal.txt                                        !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine perfoms the corrector step for turbulent kinetic       !
!  energy and generic length scale prognostic variables, tke and       !
!  gls.                                                                !
!                                                                      !
!  References:                                                         !
!                                                                      !
!  Umlauf, L. and H. Burchard, 2001:  A generic length-scale           !
!    Equation for geophysical turbulence models.                       !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: gls_corstep
      CONTAINS
!
!***********************************************************************
      SUBROUTINE gls_corstep (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
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
      CALL wclock_on (ng, iNLM, 19)
      CALL gls_corstep_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
!================================================== Thomas Roc ========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                       GRID(ng) % om_r,                           &
     &                       GRID(ng) % on_r,                           &
!========================================== End of the portion =========
     &                       GRID(ng) % Huon,                           &
     &                       GRID(ng) % Hvom,                           &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       GRID(ng) % z_r,                            &
     &                       GRID(ng) % z_w,                            &
     &                       GRID(ng) % ZoBot,                          &
     &                       OCEAN(ng) % u,                             &
     &                       OCEAN(ng) % v,                             &
     &                       OCEAN(ng) % W,                             &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                       OCEAN(ng) % dragtke,                       &
     &                       OCEAN(ng) % draggls,                       &
!========================================== End of the portion =========
     &                       FORCES(ng) % bustr,                        &
     &                       FORCES(ng) % bvstr,                        &
     &                       FORCES(ng) % sustr,                        &
     &                       FORCES(ng) % svstr,                        &
     &                       MIXING(ng) % Akt,                          &
     &                       MIXING(ng) % Akv,                          &
     &                       MIXING(ng) % bvf,                          &
     &                       MIXING(ng) % Akk,                          &
     &                       MIXING(ng) % Akp,                          &
     &                       MIXING(ng) % Lscale,                       &
     &                       MIXING(ng) % gls,                          &
     &                       MIXING(ng) % tke)
      CALL wclock_off (ng, iNLM, 19)
      RETURN
      END SUBROUTINE gls_corstep
!
!***********************************************************************
      SUBROUTINE gls_corstep_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                             om_r, on_r,                          &
     &                             Huon, Hvom, Hz, pm, pn, z_r, z_w,    &
     &                             ZoBot,                               &
     &                             u, v, W,                             &
     &                             dragtke, draggls,                    &
!========================================== End of the portion =========
     &                             bustr, bvstr, sustr, svstr,          &
     &                             Akt, Akv, bvf,                       &
     &                             Akk, Akp, Lscale, gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      USE mod_grid
!========================================== End of the portion =========
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange3d, mp_exchange4d
      USE tkebc_mod, ONLY : tkebc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: on_r(LBi:,LBj:)
!========================================== End of the portion =========
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: ZoBot(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      real(r8), intent(out) :: dragtke(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: draggls(LBi:,LBj:,0:,:)
!========================================== End of the portion ========
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: bvf(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akt(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: Akv(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akk(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Akp(LBi:,LBj:,0:)
      real(r8), intent(inout) :: Lscale(LBi:,LBj:,0:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      logical :: Lmy25
      integer :: i, itrc, j, k
      real(r8), parameter :: Gadv = 1.0_r8/3.0_r8
      real(r8), parameter :: eps = 1.0E-10_r8
      real(r8) :: Zos_min
      real(r8) :: Gh, Gm, Kprod, Ls_unlmt, Ls_lmt, Pprod, Sh, Sm
      real(r8) :: cff, cff1, cff2, cff3
      real(r8) :: cmu_fac1, cmu_fac2, cmu_fac3, cmu_fac4
      real(r8) :: gls_c3, gls_exp1, gls_fac1, gls_fac2, gls_fac3
      real(r8) :: gls_fac4, gls_fac5, gls_fac6, ql, sqrt2, strat2
      real(r8) :: tke_exp1, tke_exp2, tke_exp3, tke_exp4, wall_fac
      real(r8) :: gls_d, gls_sigp_cb, ogls_sigp, sig_eff
      real(r8) :: L_sft
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      real(r8) :: c_t, cpx, c_p, c_d, c_gls, T_k, T_gls
      real(r8) :: coef1, coef2, coef3, coef4, mask, delta, kcff
!========================================== End of the portion =========
      real(r8), dimension(IminS:ImaxS) :: tke_fluxt
      real(r8), dimension(IminS:ImaxS) :: tke_fluxb
      real(r8), dimension(IminS:ImaxS) :: gls_fluxt
      real(r8), dimension(IminS:ImaxS) :: gls_fluxb
      real(r8), dimension(IminS:ImaxS) :: Zos_eff
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: BCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCK
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FCP
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dU
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: dV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: shear2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: buoy2
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Zob_min
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: curvP
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradK
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradP
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
!  Compute several constants.
!-----------------------------------------------------------------------
!
      Zos_min=MAX(Zos(ng),0.0001_r8)
      DO j=Jstr,Jend
        DO i=Istr,Iend
          Zob_min(i,j)=MAX(ZoBot(i,j),0.0001_r8)
        END DO
      END DO
!
      Lmy25=.FALSE.
      IF ((gls_p(ng).eq.0.0_r8).and.                                    &
     &    (gls_n(ng).eq.1.0_r8).and.                                    &
     &    (gls_m(ng).eq.1.0_r8)) THEN
        Lmy25=.TRUE.
      END IF
      L_sft=vonKar
      gls_sigp_cb=gls_sigp(ng)
      ogls_sigp=1.0_r8/gls_sigp_cb
!
      sqrt2=SQRT(2.0_r8)
      cmu_fac1=gls_cmu0(ng)**(-gls_p(ng)/gls_n(ng))
      cmu_fac2=gls_cmu0(ng)**(3.0_r8+gls_p(ng)/gls_n(ng))
      cmu_fac3=1.0_r8/gls_cmu0(ng)**2.0_r8
      cmu_fac4=(1.5_r8*gls_sigk(ng))**(1.0_r8/3.0_r8)/                  &
     &         (gls_cmu0(ng)**(4.0_r8/3.0_r8))
      gls_fac1=gls_n(ng)*gls_cmu0(ng)**(gls_p(ng)+1.0_r8)
      gls_fac2=gls_cmu0(ng)**(gls_p(ng))*gls_n(ng)*                     &
     &         vonKar**(gls_n(ng))
      gls_fac3=gls_cmu0(ng)**(gls_p(ng))*gls_n(ng)
      gls_fac4=gls_cmu0(ng)**(gls_p(ng))
      gls_fac5=0.56_r8**(0.5_r8*gls_n(ng))*gls_cmu0(ng)**gls_p(ng)
      gls_fac6=8.0_r8/gls_cmu0(ng)**6.0_r8
!
      gls_exp1=1.0_r8/gls_n(ng)
      tke_exp1=gls_m(ng)/gls_n(ng)
      tke_exp2=0.5_r8+gls_m(ng)/gls_n(ng)
      tke_exp3=0.5_r8+gls_m(ng)
      tke_exp4=gls_m(ng)+0.5_r8*gls_n(ng)
!
!-----------------------------------------------------------------------
!  Compute vertical velocity shear at W-points.
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)-1
        DO j=Jstrm1,Jendp1
          DO i=Istrm1,Iendp1
            cff=0.5_r8/(z_r(i,j,k+1)-z_r(i,j,k))
            shear2(i,j,k)=(cff*(u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)+      &
     &                          u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp)))**2+ &
     &                    (cff*(v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)+      &
     &                          v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp)))**2
          END DO
        END DO
      END DO
!
! Load Brunt-Vaisala frequency.
!
      DO k=1,N(ng)-1
        DO j=Jstr-1,Jend+1
          DO i=Istr-1,Iend+1
            buoy2(i,j,k)=bvf(i,j,k)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Time-step advective terms.
!-----------------------------------------------------------------------
!
!  At entry, it is assumed that the turbulent kinetic energy fields
!  "tke" and "gls", at time level "nnew", are set to its values at
!  time level "nstp" multiplied by the grid box thicknesses Hz
!  (from old time step and at W-points).
!
      DO k=1,N(ng)-1
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            gradK(i,j)=(tke(i,j,k,3)-tke(i-1,j,k,3))
            gradP(i,j)=(gls(i,j,k,3)-gls(i-1,j,k,3))
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Istr-1,j)=gradK(Istr,j)
              gradP(Istr-1,j)=gradP(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              gradK(Iend+2,j)=gradK(Iend+1,j)
              gradP(Iend+2,j)=gradP(Iend+1,j)
            END DO
          END IF
        END IF
!
!  Third-order, upstream bias advection with velocity dependent
!  hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=Istr-1,Iend+1
            curvK(i,j)=gradK(i+1,j)-gradK(i,j)
            curvP(i,j)=gradP(i+1,j)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            cff=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i-1,j)
              cff2=curvP(i-1,j)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FXK(i,j)=cff*0.5_r8*(tke(i-1,j,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FXP(i,j)=cff*0.5_r8*(gls(i-1,j,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            gradK(i,j)=(tke(i,j,k,3)-tke(i,j-1,k,3))
            gradP(i,j)=(gls(i,j,k,3)-gls(i,j-1,k,3))
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jstr-1)=gradK(i,Jstr)
              gradP(i,Jstr-1)=gradP(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              gradK(i,Jend+2)=gradK(i,Jend+1)
              gradP(i,Jend+2)=gradP(i,Jend+1)
            END DO
          END IF
        END IF
        DO j=Jstr-1,Jend+1
          DO i=Istr,Iend
            curvK(i,j)=gradK(i,j+1)-gradK(i,j)
            curvP(i,j)=gradP(i,j+1)-gradP(i,j)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            cff=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            IF (cff.gt.0.0_r8) THEN
              cff1=curvK(i,j-1)
              cff2=curvP(i,j-1)
            ELSE
              cff1=curvK(i,j)
              cff2=curvP(i,j)
            END IF
            FEK(i,j)=cff*0.5_r8*(tke(i,j-1,k,3)+tke(i,j,k,3)-           &
     &                           Gadv*cff1)
            FEP(i,j)=cff*0.5_r8*(gls(i,j-1,k,3)+gls(i,j,k,3)-           &
     &                           Gadv*cff2)
          END DO
        END DO
!
!  Time-step horizontal advection.
!
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FXK(i+1,j)-FXK(i,j)+                   &
     &                           FEK(i,j+1)-FEK(i,j))
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FXP(i+1,j)-FXP(i,j)+                   &
     &                           FEP(i,j+1)-FEP(i,j))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
          END DO
        END DO
      END DO
!
! Compute vertical advection.
!
      DO j=Jstr,Jend
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            cff=0.5*(W(i,j,k)+W(i,j,k-1))
            FCK(i,k)=cff*(cff1*(tke(i,j,k-1,3)+                         &
     &                          tke(i,j,k  ,3))-                        &
     &                    cff2*(tke(i,j,k-2,3)+                         &
     &                          tke(i,j,k+1,3)))
            FCP(i,k)=cff*(cff1*(gls(i,j,k-1,3)+                         &
     &                          gls(i,j,k  ,3))-                        &
     &                    cff2*(gls(i,j,k-2,3)+                         &
     &                          gls(i,j,k+1,3)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
         DO i=Istr,Iend
          cff=0.5_r8*(W(i,j,0)+W(i,j,1))
          FCK(i,1)=cff*(cff1*tke(i,j,0,3)+                              &
     &                  cff2*tke(i,j,1,3)-                              &
     &                  cff3*tke(i,j,2,3))
          FCP(i,1)=cff*(cff1*gls(i,j,0,3)+                              &
     &                  cff2*gls(i,j,1,3)-                              &
     &                  cff3*gls(i,j,2,3))
          cff=0.5_r8*(W(i,j,N(ng))+W(i,j,N(ng)-1))
          FCK(i,N(ng))=cff*(cff1*tke(i,j,N(ng)  ,3)+                    &
     &                      cff2*tke(i,j,N(ng)-1,3)-                    &
     &                      cff3*tke(i,j,N(ng)-2,3))
          FCP(i,N(ng))=cff*(cff1*gls(i,j,N(ng)  ,3)+                    &
     &                      cff2*gls(i,j,N(ng)-1,3)-                    &
     &                      cff3*gls(i,j,N(ng)-2,3))
        END DO
!
!  Time-step vertical advection term.
!
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff=dt(ng)*pm(i,j)*pn(i,j)
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-                            &
     &                      cff*(FCK(i,k+1)-FCK(i,k))
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-                            &
     &                      cff*(FCP(i,k+1)-FCP(i,k))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
          END DO
        END DO
!
!----------------------------------------------------------------------
!  Compute vertical mixing, turbulent production and turbulent
!  dissipation terms.
!----------------------------------------------------------------------
!
!  Set term for vertical mixing of turbulent fields.
!
        cff=-0.5_r8*dt(ng)
        DO i=Istr,Iend
          DO k=2,N(ng)-1
            FCK(i,k)=cff*(Akk(i,j,k)+Akk(i,j,k-1))/Hz(i,j,k)
            FCP(i,k)=cff*(Akp(i,j,k)+Akp(i,j,k-1))/Hz(i,j,k)
            CF(i,k)=0.0_r8
          END DO
          FCP(i,1)=0.0_r8
          FCP(i,N(ng))=0.0_r8
          FCK(i,1)=0.0_r8
          FCK(i,N(ng))=0.0_r8
        END DO
!
!  Compute production and dissipation terms.
!
        DO i=Istr,Iend
          DO k=1,N(ng)-1
!
!  Compute shear and bouyant production of turbulent energy (m3/s3)
!  at W-points (ignore small negative values of buoyancy).
!
            strat2=buoy2(i,j,k)
            IF (strat2.gt.0.0_r8) THEN
              gls_c3=gls_c3m(ng)
            ELSE
              gls_c3=gls_c3p(ng)
            END IF
            Kprod=shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-               &
     &            strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
            Pprod=gls_c1(ng)*shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng))-    &
     &            gls_c3*strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
!
!  If negative production terms, then add buoyancy to dissipation terms
!  (BCK and BCP) below, using "cff1" and "cff2" as the on/off switch.
!
            cff1=1.0_r8
            IF (Kprod.lt.0.0_r8) THEN
              Kprod=Kprod+strat2*(Akt(i,j,k,itemp)-Akt_bak(itemp,ng))
              cff1=0.0_r8
            END IF
            cff2=1.0_r8
            IF (Pprod.lt.0.0_r8) THEN
              Pprod=Pprod+gls_c3*strat2*(Akt(i,j,k,itemp)-              &
     &                                   Akt_bak(itemp,ng))
              cff2=0.0_r8
            END IF
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
!------------------------Turbine features-------------------------------
          delta=om_r(i,j)
          kcff=1000.0_r8!density anomaly correction
          c_t=GRID(ng)%rmask_t(i,j,k)
          cpx=c_t*((1-c_t)**0.5_r8)!=Ct*root(1-Ct)*lcpa
          coef4=kcff*GRID(ng)%rmask3t(i,j,k)!Rhie-Show Smearing correction
          c_p=GRID(ng)%rmask1t(i,j,k)
          c_d=4.0_r8*((1.0_r8-((1.0_r8-c_t)**0.5_r8))/                  &
     &               (1.0_r8+((1.0_r8-c_t)**0.5_r8))) 
          c_gls=GRID(ng)%rmask2t(i,j,k)          
!--------------Turbulence length scale correction term------------------
          coef1=(gls_cmu0(ng)**3.0_r8)*((tke(i,j,k,nstp)**1.5_r8)       &
     &           /Lscale(i,j,k))
          coef2=(shear2(i,j,k)*(Akv(i,j,k)-Akv_bak(ng)))**2.0_r8
          T_gls=((cpx/delta)**2.0_r8)*(c_gls*coef2)/coef1 
!--------------Turbulent kinetic energy correction term------------ 
          coef3=ABS(0.25_r8*(u(i,j,k,nnew)+u(i,j,k+1,nnew)              &
     &             + u(i,j+1,k,nnew)+u(i,j+1,k+1,nnew)))
          coef3=coef3+(ABS(0.25_r8*(u(i-1,j,k,nstp)+u(i-1,j,k+1,nstp)   &
     &             + u(i-1,j+1,k,nstp)+u(i-1,j+1,k+1,nstp))))/2.0_r8
          T_k=coef4*(cpx/(delta**2.0_r8))*((c_p*(coef3**3.0_r8))        &
     &        -(c_d*coef3*tke(i,j,k,nstp)))
!
!  Time-step shear and buoyancy production terms.
!
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+                            &
     &                      dt(ng)*cff*Kprod
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+                            &
     &                      dt(ng)*cff*Pprod*gls(i,j,k,nstp)/           &
     &                      MAX(tke(i,j,k,nstp),gls_Kmin(ng))
            tke(i,j,k,nnew)=tke(i,j,k,nnew)+cff*dt(ng)*T_k
            gls(i,j,k,nnew)=gls(i,j,k,nnew)+cff*dt(ng)*T_gls
            dragtke(i,j,k,nnew)=cff*dt(ng)*T_k
            draggls(i,j,k,nnew)=cff*dt(ng)*T_gls
!========================================== End of the portion =========
!
!  Compute dissipation of turbulent energy (m3/s3).
!
            wall_fac=1.0_r8
            IF (Lmy25) THEN
!
!  Parabolic wall function,  L = ds db / (ds + db).
!
!
            wall_fac=1.0_r8+gls_E2/(vonKar*vonKar)*                     &
     &                (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*           &
     &                 tke(i,j,k,nstp)**(-tke_exp1)*                    &
     &                 (1.0_r8/ (z_w(i,j,k)-z_w(i,j,0))))**2+           &
     &                0.25_r8/(vonKar*vonKar)*                          &
     &                (gls(i,j,k,nstp)**( gls_exp1)*cmu_fac1*           &
     &                 tke(i,j,k,nstp)**(-tke_exp1)*                    &
     &                 (1.0_r8/ (z_w(i,j,N(ng))-z_w(i,j,k))))**2
            END IF
!
            BCK(i,k)=cff*(1.0_r8+dt(ng)*                                &
     &                    gls(i,j,k,nstp)**(-gls_exp1)*cmu_fac2*        &
     &                    tke(i,j,k,nstp)**( tke_exp2)+                 &
     &                    dt(ng)*(1.0_r8-cff1)*strat2*                  &
     &                    (Akt(i,j,k,itemp)-Akt_bak(itemp,ng))/         &
     &                    tke(i,j,k,nstp))-                             &
     &                    FCK(i,k)-FCK(i,k+1)
            BCP(i,k)=cff*(1.0_r8+dt(ng)*gls_c2(ng)*wall_fac*            &
     &                    gls(i,j,k,nstp)**(-gls_exp1)*cmu_fac2*        &
     &                    tke(i,j,k,nstp)**( tke_exp2)+                 &
     &                    dt(ng)*(1.0_r8-cff2)*gls_c3*strat2*           &
     &                    (Akt(i,j,k,itemp)-Akt_bak(itemp,ng))/         &
     &                    tke(i,j,k,nstp))-                             &
     &                    FCP(i,k)-FCP(i,k+1)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Time-step dissipation and vertical diffusion terms implicitly.
!-----------------------------------------------------------------------
!
!  Set Dirichlet surface and bottom boundary conditions. Compute
!  surface roughness from wind stress (Charnok) and set Craig and
!  Banner wave breaking surface flux, if appropriate.
!
        DO i=Istr,Iend
          tke(i,j,N(ng),nnew)=MAX(cmu_fac3*0.5_r8*                      &
     &                            SQRT((sustr(i,j)+sustr(i+1,j))**2+    &
     &                                 (svstr(i,j)+svstr(i,j+1))**2),   &
     &                            gls_Kmin(ng))
          tke(i,j,0,nnew)=MAX(cmu_fac3*0.5_r8*                          &
     &                        SQRT((bustr(i,j)+bustr(i+1,j))**2+        &
     &                             (bvstr(i,j)+bvstr(i,j+1))**2),       &
     &                        gls_Kmin(ng))
          Zos_eff(i)=Zos_min
          gls(i,j,N(ng),nnew)=MAX(gls_cmu0(ng)**gls_p(ng)*              &
     &                            tke(i,j,N(ng),nnew)**gls_m(ng)*       &
     &                            (L_sft*Zos_eff(i))**gls_n(ng),        &
     &                            gls_Pmin(ng))
          cff=gls_fac4*(vonKar*Zob_min(i,j))**(gls_n(ng))
          gls(i,j,0,nnew)=MAX(cff*tke(i,j,0,nnew)**(gls_m(ng)),         &
     &                        gls_Pmin(ng))
        END DO
!
!  Solve tri-diagonal system for turbulent kinetic energy.
!
        DO i=Istr,Iend
          tke_fluxt(i)=0.0_r8
          tke_fluxb(i)=0.0_r8
!
          cff=1.0_r8/BCK(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCK(i,N(ng)-1)
          tke(i,j,N(ng)-1,nnew)=cff*(tke(i,j,N(ng)-1,nnew)+tke_fluxt(i))
        END DO
        DO i=Istr,Iend
          DO k=N(ng)-2,1,-1
            cff=1.0_r8/(BCK(i,k)-CF(i,k+1)*FCK(i,k+1))
            CF(i,k)=cff*FCK(i,k)
            tke(i,j,k,nnew)=cff*(tke(i,j,k,nnew)-                       &
     &                           FCK(i,k+1)*tke(i,j,k+1,nnew))
          END DO
          tke(i,j,1,nnew)=tke(i,j,1,nnew)-cff*tke_fluxb(i)
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            tke(i,j,k,nnew)=tke(i,j,k,nnew)-CF(i,k)*tke(i,j,k-1,nnew)
          END DO
        END DO
!
!  Solve tri-diagonal system for generic statistical field.
!
        DO i=Istr,Iend
          cff=0.5_r8*(tke(i,j,N(ng),nnew)+tke(i,j,N(ng)-1,nnew))
          gls_fluxt(i)=dt(ng)*gls_fac3*cff**gls_m(ng)*                  &
     &                 L_sft**(gls_n(ng))*                              &
     &                 (Zos_eff(i)+0.5_r8*Hz(i,j,N(ng)))**              &
     &                 (gls_n(ng)-1.0_r8)*                              &
     &                 0.5_r8*(Akp(i,j,N(ng))+Akp(i,j,N(ng)-1))
          cff=0.5_r8*(tke(i,j,0,nnew)+tke(i,j,1,nnew))
          gls_fluxb(i)=dt(ng)*gls_fac2*(cff**gls_m(ng))*                &
     &                 (0.5_r8*Hz(i,j,1)+Zob_min(i,j))**                &
     &                 (gls_n(ng)-1.0_r8)*                              &
     &                 0.5_r8*(Akp(i,j,0)+Akp(i,j,1))
!
          cff=1.0_r8/BCP(i,N(ng)-1)
          CF(i,N(ng)-1)=cff*FCP(i,N(ng)-1)
          gls(i,j,N(ng)-1,nnew)=cff*(gls(i,j,N(ng)-1,nnew)-gls_fluxt(i))
        END DO
        DO i=Istr,Iend
          DO k=N(ng)-2,1,-1
            cff=1.0_r8/(BCP(i,k)-CF(i,k+1)*FCP(i,k+1))
            CF(i,k)=cff*FCP(i,k)
            gls(i,j,k,nnew)=cff*(gls(i,j,k,nnew)-                       &
     &                           FCP(i,k+1)*gls(i,j,k+1,nnew))
          END DO
          gls(i,j,1,nnew)=gls(i,j,1,nnew)-cff*gls_fluxb(i)
        END DO
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            gls(i,j,k,nnew)=gls(i,j,k,nnew)-CF(i,k)*gls(i,j,k-1,nnew)
          END DO
        END DO
!
!---------------------------------------------------------------------
!  Compute vertical mixing coefficients (m2/s).
!---------------------------------------------------------------------
!
        DO i=Istr,Iend
          DO k=1,N(ng)-1
!
!  Compute turbulent length scale (m).
!
            tke(i,j,k,nnew)=MAX(tke(i,j,k,nnew),gls_Kmin(ng))
            gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_Pmin(ng))
            IF (gls_n(ng).ge.0.0_r8) THEN
              gls(i,j,k,nnew)=MIN(gls(i,j,k,nnew),gls_fac5*             &
     &                            tke(i,j,k,nnew)**(tke_exp4)*          &
     &                            (SQRT(MAX(0.0_r8,                     &
     &                                  buoy2(i,j,k)))+eps)**           &
     &                            (-gls_n(ng)))
            ELSE
              gls(i,j,k,nnew)=MAX(gls(i,j,k,nnew),gls_fac5*             &
     &                            tke(i,j,k,nnew)**(tke_exp4)*          &
     &                            (SQRT(MAX(0.0_r8,                     &
     &                                  buoy2(i,j,k)))+eps)**           &
     &                            (-gls_n(ng)))
            END IF
            Ls_unlmt=MAX(eps,                                           &
     &                   gls(i,j,k,nnew)**( gls_exp1)*cmu_fac1*         &
     &                   tke(i,j,k,nnew)**(-tke_exp1))
            IF (buoy2(i,j,k).gt.0.0_r8) THEN
              Ls_lmt=MIN(Ls_unlmt,                                      &
     &                 SQRT(0.56_r8*tke(i,j,k,nnew)/                    &
     &                      (MAX(0.0_r8,buoy2(i,j,k))+eps)))
            ELSE
              Ls_lmt=Ls_unlmt
            END IF
!
! Recompute gls based on limited length scale
!
            gls(i,j,k,nnew)=MAX(gls_cmu0(ng)**gls_p(ng)*                &
     &                          tke(i,j,k,nnew)**gls_m(ng)*             &
     &                          Ls_lmt**gls_n(ng), gls_Pmin(ng))
!
!  Compute nondimensional stability functions for tracers (Sh) and
!  momentum (Sm).
!
            Gh=MIN(gls_Gh0,-buoy2(i,j,k)*Ls_lmt*Ls_lmt/                 &
     &                    (2.0_r8*tke(i,j,k,nnew)))
            Gh=MIN(Gh,Gh-(Gh-gls_Ghcri)**2/                             &
     &                    (Gh+gls_Gh0-2.0_r8*gls_Ghcri))
            Gh=MAX(Gh,gls_Ghmin)
!
!  Compute shear number.
!
            Gm=(gls_b0/gls_fac6-gls_b1*Gh+gls_b3*gls_fac6*(Gh**2))/     &
     &         (gls_b2-gls_b4*gls_fac6*Gh)
            Gm=MIN(Gm,shear2(i,j,k)*Ls_lmt*Ls_lmt/                      &
     &                    (2.0_r8*tke(i,j,k,nnew)))
!
!  Compute stability functions
!
            cff=gls_b0-gls_b1*gls_fac6*Gh+gls_b2*gls_fac6*Gm+           &
     &          gls_b3*gls_fac6**2*Gh**2-gls_b4*gls_fac6**2*Gh*Gm+      &
     &          gls_b5*gls_fac6**2*Gm*Gm
            Sm=(gls_s0-gls_s1*gls_fac6*Gh+gls_s2*gls_fac6*Gm)/cff
            Sh=(gls_s4-gls_s5*gls_fac6*Gh+gls_s6*gls_fac6*Gm)/cff
            Sm=MAX(Sm,0.0_r8)
            Sh=MAX(Sh,0.0_r8)
!
!  Relate Canuto stability to ROMS notation
!
            Sm=Sm*sqrt2/gls_cmu0(ng)**3
            Sh=Sh*sqrt2/gls_cmu0(ng)**3
!
!  Compute vertical mixing (m2/s) coefficients of momentum and
!  tracers.  Average ql over the two timesteps rather than using
!  the new Lscale and just averaging tke.
!
            ql=sqrt2*0.5_r8*(Ls_lmt*SQRT(tke(i,j,k,nnew))+              &
     &                       Lscale(i,j,k)*SQRT(tke(i,j,k,nstp)))
            Akv(i,j,k)=Akv_bak(ng)+Sm*ql
            DO itrc=1,NAT
              Akt(i,j,k,itrc)=Akt_bak(itrc,ng)+Sh*ql
            END DO
!
!  Compute vertical mixing (m2/s) coefficents of turbulent kinetic
!  energy and generic statistical field.
!
            Akk(i,j,k)=Akk_bak(ng)+                                     &
     &                 Sm*ql/gls_sigk(ng)
            Akp(i,j,k)=Akp_bak(ng)+Sm*ql*ogls_sigp
!
!  Save limited length scale.
!
            Lscale(i,j,k)=Ls_lmt
          END DO
!
!  Compute vertical mixing coefficients at the surface and bottom.
!
!
          Akv(i,j,N(ng))=Akv_bak(ng)+L_sft*Zos_eff(i)*gls_cmu0(ng)*     &
     &                   SQRT(tke(i,j,N(ng),nnew))
          Akv(i,j,0)=Akv_bak(ng)+vonKar*Zob_min(i,j)*gls_cmu0(ng)*      &
     &               SQRT(tke(i,j,0,nnew))
!
          Akk(i,j,N(ng))=Akk_bak(ng)+Akv(i,j,N(ng))/gls_sigk(ng)
          Akk(i,j,0)=Akk_bak(ng)+Akv(i,j,0)/gls_sigk(ng)
          Akp(i,j,N(ng))=Akp_bak(ng)+Akv(i,j,N(ng))*ogls_sigp
          Akp(i,j,0)=Akp_bak(ng)+Akv(i,j,0)/gls_sigp(ng)
!
          DO itrc=1,NAT
            Akt(i,j,N(ng),itrc)=Akt_bak(itrc,ng)
            Akt(i,j,0,itrc)=Akt_bak(itrc,ng)
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nnew, nstp,                                      &
     &                 gls, tke)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,nnew))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          Akv)
        DO itrc=1,NAT
          CALL exchange_w3d_tile (ng, tile,                             &
     &                            LBi, UBi, LBj, UBj, 0, N(ng),         &
     &                            Akt(:,:,:,itrc))
        END DO
      END IF
      CALL mp_exchange3d (ng, tile, iNLM, 3,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,nnew),                              &
     &                    gls(:,:,:,nnew),                              &
     &                    Akv)
      CALL mp_exchange4d (ng, tile, iNLM, 1,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng), 1, NAT,         &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    Akt)
      RETURN
      END SUBROUTINE gls_corstep_tile
      END MODULE gls_corstep_mod

      MODULE rhs3d_mod
!
!svn $Id: rhs3d.F 719 2014-03-13 22:25:13Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!                                                                      !
!  With Portions Copyright (c) 2014 Thomas Roc and ITPower             !
!    Licensed under an Affero GPL style license                        !
!    See License_ROMS-Tidal.txt                                        !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  This subroutine evaluates right-hand-side terms for 3D momentum     !
!  and tracers equations.                                              !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC  :: rhs3d
!
      CONTAINS
!
!***********************************************************************
      SUBROUTINE rhs3d (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_coupling
      USE mod_forces
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
!
      USE pre_step3d_mod, ONLY : pre_step3d
      USE prsgrd_mod, ONLY : prsgrd
      USE uv3dmix_mod, ONLY : uv3dmix2
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
!-----------------------------------------------------------------------
!  Initialize computations for new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
!
      CALL pre_step3d (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
!-----------------------------------------------------------------------
!
      CALL prsgrd (ng, tile)
!
!-----------------------------------------------------------------------
!  Compute right-hand-side terms for the 3D momentum equations.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 21)
      CALL rhs3d_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj,                              &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nrhs(ng),                                        &
     &                 GRID(ng) % Hz,                                   &
     &                 GRID(ng) % Huon,                                 &
     &                 GRID(ng) % Hvom,                                 &
     &                 GRID(ng) % fomn,                                 &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                 GRID(ng) % om_r,                                 &
     &                 GRID(ng) % om_u,                                 &
     &                 GRID(ng) % om_v,                                 &
     &                 GRID(ng) % on_u,                                 &
     &                 GRID(ng) % on_v,                                 &
     &                 GRID(ng) % pm,                                   &
     &                 GRID(ng) % pn,                                   &
     &                 GRID(ng) % imask_t,                              &
     &                 GRID(ng) % umask_t,                              &
     &                 GRID(ng) % vmask_t,                              &
     &                 OCEAN(ng) % dragfrc,                             &
     &                 OCEAN(ng) % dragpwr,                             &
! End of the portion                                                   !
!=======================================================================
     &                 FORCES(ng) % bustr,                              &
     &                 FORCES(ng) % bvstr,                              &
     &                 FORCES(ng) % sustr,                              &
     &                 FORCES(ng) % svstr,                              &
     &                 OCEAN(ng) % u,                                   &
     &                 OCEAN(ng) % v,                                   &
     &                 OCEAN(ng) % W,                                   &
     &                 COUPLING(ng) % rufrc,                            &
     &                 COUPLING(ng) % rvfrc,                            &
     &                 OCEAN(ng) % ru,                                  &
     &                 OCEAN(ng) % rv)
      CALL wclock_off (ng, iNLM, 21)
!
!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
!
      CALL uv3dmix2 (ng, tile)
      RETURN
      END SUBROUTINE rhs3d
!
!***********************************************************************
      SUBROUTINE rhs3d_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nrhs,                                      &
     &                       Hz, Huon, Hvom,                            &
     &                       fomn,                                      &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                       om_r,                                      &
! End of the portion                                                   !
!=======================================================================
     &                       om_u, om_v, on_u, on_v, pm, pn,            &
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
     &                       imask_t, umask_t, vmask_t,                 &
     &                       dragfrc, dragpwr,                          &
! End of the portion                                                   !
!=======================================================================
     &                       bustr, bvstr,                              &
     &                       sustr, svstr,                              &
     &                       u, v, W,                                   &
     &                       rufrc, rvfrc,                              &
     &                       ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_clima
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
!
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: fomn(LBi:,LBj:)
      real(r8), intent(in) :: om_u(LBi:,LBj:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: on_v(LBi:,LBj:)
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      real(r8), intent(in) :: om_r(LBi:,LBj:)
      real(r8), intent(in) :: imask_t(LBi:,LBj:,:)
      real(r8), intent(in) :: umask_t(LBi:,LBj:,:)
      real(r8), intent(in) :: vmask_t(LBi:,LBj:,:)
      real(r8), intent(out) :: dragfrc(LBi:,LBj:,:,:)
      real(r8), intent(out) :: dragpwr(LBi:,LBj:,:,:)
! End of the portion                                                   !
!=======================================================================
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: bustr(LBi:,LBj:)
      real(r8), intent(in) :: bvstr(LBi:,LBj:)
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
      real(r8), intent(in) :: u(LBi:,LBj:,:,:)
      real(r8), intent(in) :: v(LBi:,LBj:,:,:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
      real(r8), intent(out) :: rufrc(LBi:,LBj:)
      real(r8), intent(out) :: rvfrc(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: Gadv = -0.25_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8) :: fac, fac1, fac2
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
      integer :: indice
      real(r8) :: coef1, coef2, coef3
      real(r8) :: cd, kcff
      real(r8) :: drag_f, drag_t, dx, dy
      real(r8) :: u_mid, v_mid, we
      real(r8) :: Hlay, surf_u, surf_v, surf_cell, surf_cyl
! End of the portion                                                   !
!=======================================================================
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: DC
      real(r8), dimension(IminS:ImaxS,0:N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Huxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Hvxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: UFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Uwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: VFe
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: Vwrk
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: uxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vee
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: vxx
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
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
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
! End of the portion                                                   !
!=======================================================================
!
!-----------------------------------------------------------------------
!  Apply surface stress as a bodyforce: determine the thickness (m)
!  of the surface layer; then add in surface stress as a bodyfoce.
!-----------------------------------------------------------------------
!
!================================================== Thomas Roc =========
! The following portion is Copyright (c) 2014 Thomas Roc and ITPower   !
!
!  Apply mid-depth stress as a bodyforce: determine the thickness (m)
!  of the mid-depth layer; then add in mid-depth stress as a bodyfoce.
!
      DO j=JstrV-1,Jend
        DO i=IstrU-1,Iend
          wrk(i,j)=0.0_r8
        END DO
      END DO
!Start K-LOOP
      DO k=1,N(ng)
        DO j=Jstr,Jend
          DO i=IstrU,Iend
!
!Pressure drop in term of stress in u-direction, BODYFORCE formulation 
!imposes such a formulation
!
!Turbine is located at the rho point therefore a 
!split implementation between the surrounding u points is needed.
!
!turbine features
          cd=umask_t(i,j,k)
          kcff=1000.0_r8 !Density anomaly correction
          indice=imask_t(i,j,k)
          dy=on_u(i,j)
          Hlay=0.5_r8*(Hz(i-1,j,k)+Hz(i,j,k)) !Vertical layer thickness
          surf_u=dy*Hlay
!
!--------------------Rhie/Show implementation scheme---------------------
          cff1=0.5_r8*(u((i+(indice-1)),j,k,nrhs)+u((i+indice),j,k,nrhs))
          cff2=0.25_r8*(v(i  ,j  ,k,nrhs)+                              &
     &                  v(i  ,j+1,k,nrhs)+                              &
     &                  v(i-1,j  ,k,nrhs)+                              &
     &                  v(i-1,j+1,k,nrhs))
          u_mid=SQRT((cff1*cff1)+(cff2*cff2))
          drag_t=(kcff/om_u(i+(indice-1),j))                            &
     &         *0.25_r8*cd*u_mid*cff1*surf_u
          dragfrc(i,j,k,nrhs)=drag_t
          dragpwr(i,j,k,nrhs)=drag_t*u_mid
          ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-drag_t 
          END DO
        END DO
       END DO
! End of the portion                                                   !
!=======================================================================
!
      K_LOOP : DO k=1,N(ng)
!
!-----------------------------------------------------------------------
!  Add in nudging of 3D momentum climatology.
!-----------------------------------------------------------------------
!
        IF (LnudgeM3CLM(ng)) THEN
          DO j=Jstr,Jend
            DO i=IstrU,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i-1,j,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i  ,j,k))*               &
     &            om_u(i,j)*on_u(i,j)
              ru(i,j,k,nrhs)=ru(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i-1,j,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%uclm(i,j,k)-                    &
     &                        u(i,j,k,nrhs))
            END DO
          END DO
          DO j=JstrV,Jend
            DO i=Istr,Iend
              cff=0.25_r8*(CLIMA(ng)%M3nudgcof(i,j-1,k)+                &
     &                     CLIMA(ng)%M3nudgcof(i,j  ,k))*               &
     &            om_v(i,j)*on_v(i,j)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)+                            &
     &                       cff*(Hz(i,j-1,k)+Hz(i,j,k))*               &
     &                       (CLIMA(ng)%vclm(i,j,k)-                    &
     &                        v(i,j,k,nrhs))
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Add in horizontal advection of momentum.
!-----------------------------------------------------------------------
!
!  Compute diagonal [UFx,VFe] and off-diagonal [UFe,VFx] components
!  of tensor of momentum flux due to horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrUm1,Iendp1
            uxx(i,j)=u(i-1,j,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
     &               u(i+1,j,k,nrhs)
            Huxx(i,j)=Huon(i-1,j,k)-2.0_r8*Huon(i,j,k)+Huon(i+1,j,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Istr,j)=uxx (Istr+1,j)
              Huxx(Istr,j)=Huxx(Istr+1,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              uxx (Iend+1,j)=uxx (Iend,j)
              Huxx(Iend+1,j)=Huxx(Iend,j)
            END DO
          END IF
        END IF
!
!  Third-order, upstream bias u-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=Jstr,Jend
          DO i=IstrU-1,Iend
            cff1=u(i  ,j,k,nrhs)+                                       &
     &           u(i+1,j,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=uxx(i,j)
            ELSE
              cff=uxx(i+1,j)
            END IF
            UFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Huon(i  ,j,k)+                                    &
     &                Huon(i+1,j,k)+                                    &
     &                Gadv*0.5_r8*(Huxx(i  ,j)+                         &
     &                             Huxx(i+1,j)))
          END DO
        END DO
        DO j=Jstrm1,Jendp1
          DO i=IstrU,Iend
            uee(i,j)=u(i,j-1,k,nrhs)-2.0_r8*u(i,j,k,nrhs)+              &
     &               u(i,j+1,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jstr-1)=uee(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=IstrU,Iend
              uee(i,Jend+1)=uee(i,Jend)
            END DO
          END IF
        END IF
        DO j=Jstr,Jend+1
          DO i=IstrU-1,Iend
           Hvxx(i,j)=Hvom(i-1,j,k)-2.0_r8*Hvom(i,j,k)+Hvom(i+1,j,k)
          END DO
        END DO
        DO j=Jstr,Jend+1
          DO i=IstrU,Iend
            cff1=u(i,j  ,k,nrhs)+                                       &
     &           u(i,j-1,k,nrhs)
            cff2=Hvom(i,j,k)+Hvom(i-1,j,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=uee(i,j-1)
            ELSE
              cff=uee(i,j)
            END IF
            UFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Hvxx(i  ,j)+                    &
     &                                  Hvxx(i-1,j)))
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istrm1,Iendp1
            vxx(i,j)=v(i-1,j,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i+1,j,k,nrhs)
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Istr-1,j)=vxx(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=JstrV,Jend
              vxx(Iend+1,j)=vxx(Iend,j)
            END DO
          END IF
        END IF
        DO j=JstrV-1,Jend
          DO i=Istr,Iend+1
           Huee(i,j)=Huon(i,j-1,k)-2.0_r8*Huon(i,j,k)+Huon(i,j+1,k)
          END DO
        END DO
!
!  Third-order, upstream bias v-momentum advection with velocity
!  dependent hyperdiffusion.
!
        DO j=JstrV,Jend
          DO i=Istr,Iend+1
            cff1=v(i  ,j,k,nrhs)+                                       &
     &           v(i-1,j,k,nrhs)
            cff2=Huon(i,j,k)+Huon(i,j-1,k)
            IF (cff2.gt.0.0_r8) THEN
              cff=vxx(i-1,j)
            ELSE
              cff=vxx(i,j)
            END IF
            VFx(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (cff2+Gadv*0.5_r8*(Huee(i,j  )+                    &
     &                                  Huee(i,j-1)))
          END DO
        END DO
        DO j=JstrVm1,Jendp1
          DO i=Istr,Iend
            vee(i,j)=v(i,j-1,k,nrhs)-2.0_r8*v(i,j,k,nrhs)+              &
     &               v(i,j+1,k,nrhs)
            Hvee(i,j)=Hvom(i,j-1,k)-2.0_r8*Hvom(i,j,k)+Hvom(i,j+1,k)
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jstr)=vee (i,Jstr+1)
              Hvee(i,Jstr)=Hvee(i,Jstr+1)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              vee (i,Jend+1)=vee (i,Jend)
              Hvee(i,Jend+1)=Hvee(i,Jend)
            END DO
          END IF
        END IF
        DO j=JstrV-1,Jend
          DO i=Istr,Iend
            cff1=v(i,j  ,k,nrhs)+                                       &
     &           v(i,j+1,k,nrhs)
            IF (cff1.gt.0.0_r8) THEN
              cff=vee(i,j)
            ELSE
              cff=vee(i,j+1)
            END IF
            VFe(i,j)=0.25_r8*(cff1+Gadv*cff)*                           &
     &               (Hvom(i,j  ,k)+                                    &
     &                Hvom(i,j+1,k)+                                    &
     &                Gadv*0.5_r8*(Hvee(i,j  )+                         &
     &                             Hvee(i,j+1)))
          END DO
        END DO
!
!  Add in horizontal advection.
!
        DO j=Jstr,Jend
          DO i=IstrU,Iend
            cff1=UFx(i,j)-UFx(i-1,j)
            cff2=UFe(i,j+1)-UFe(i,j)
            cff=cff1+cff2
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
          END DO
        END DO
        DO j=JstrV,Jend
          DO i=Istr,Iend
            cff1=VFx(i+1,j)-VFx(i,j)
            cff2=VFe(i,j)-VFe(i,j-1)
            cff=cff1+cff2
            rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
          END DO
        END DO
      END DO K_LOOP
!
      J_LOOP : DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Add in vertical advection.
!-----------------------------------------------------------------------
!
        cff1=9.0_r8/16.0_r8
        cff2=1.0_r8/16.0_r8
        DO k=2,N(ng)-2
          DO i=IstrU,Iend
            FC(i,k)=(cff1*(u(i,j,k  ,nrhs)+                             &
     &                     u(i,j,k+1,nrhs))-                            &
     &               cff2*(u(i,j,k-1,nrhs)+                             &
     &                     u(i,j,k+2,nrhs)))*                           &
     &              (cff1*(W(i  ,j,k)+                                  &
     &                     W(i-1,j,k))-                                 &
     &               cff2*(W(i+1,j,k)+                                  &
     &                     W(i-2,j,k)))
          END DO
        END DO
        DO i=IstrU,Iend
          FC(i,N(ng))=0.0_r8
          FC(i,N(ng)-1)=(cff1*(u(i,j,N(ng)-1,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs))-                    &
     &                   cff2*(u(i,j,N(ng)-2,nrhs)+                     &
     &                         u(i,j,N(ng)  ,nrhs)))*                   &
     &                  (cff1*(W(i  ,j,N(ng)-1)+                        &
     &                         W(i-1,j,N(ng)-1))-                       &
     &                   cff2*(W(i+1,j,N(ng)-1)+                        &
     &                         W(i-2,j,N(ng)-1)))
          FC(i,1)=(cff1*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,2,nrhs))-                                &
     &             cff2*(u(i,j,1,nrhs)+                                 &
     &                   u(i,j,3,nrhs)))*                               &
     &            (cff1*(W(i  ,j,1)+                                    &
     &                   W(i-1,j,1))-                                   &
     &             cff2*(W(i+1,j,1)+                                    &
     &                   W(i-2,j,1)))
          FC(i,0)=0.0_r8
        END DO
        DO k=1,N(ng)
          DO i=IstrU,Iend
            cff=FC(i,k)-FC(i,k-1)
            ru(i,j,k,nrhs)=ru(i,j,k,nrhs)-cff
          END DO
        END DO
        IF (j.ge.JstrV) THEN
          cff1=9.0_r8/16.0_r8
          cff2=1.0_r8/16.0_r8
          DO k=2,N(ng)-2
            DO i=Istr,Iend
              FC(i,k)=(cff1*(v(i,j,k  ,nrhs)+                           &
     &                       v(i,j,k+1,nrhs))-                          &
     &                 cff2*(v(i,j,k-1,nrhs)+                           &
     &                       v(i,j,k+2,nrhs)))*                         &
     &                (cff1*(W(i,j  ,k)+                                &
     &                       W(i,j-1,k))-                               &
     &                 cff2*(W(i,j+1,k)+                                &
     &                       W(i,j-2,k)))
            END DO
          END DO
          DO i=Istr,Iend
            FC(i,N(ng))=0.0_r8
            FC(i,N(ng)-1)=(cff1*(v(i,j,N(ng)-1,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs))-                  &
     &                     cff2*(v(i,j,N(ng)-2,nrhs)+                   &
     &                           v(i,j,N(ng)  ,nrhs)))*                 &
     &                    (cff1*(W(i,j  ,N(ng)-1)+                      &
     &                           W(i,j-1,N(ng)-1))-                     &
     &                     cff2*(W(i,j+1,N(ng)-1)+                      &
     &                           W(i,j-2,N(ng)-1)))
            FC(i,1)=(cff1*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,2,nrhs))-                              &
     &               cff2*(v(i,j,1,nrhs)+                               &
     &                     v(i,j,3,nrhs)))*                             &
     &              (cff1*(W(i,j  ,1)+                                  &
     &                     W(i,j-1,1))-                                 &
     &               cff2*(W(i,j+1,1)+                                  &
     &                     W(i,j-2,1)))
            FC(i,0)=0.0_r8
          END DO
          DO k=1,N(ng)
            DO i=Istr,Iend
              cff=FC(i,k)-FC(i,k-1)
              rv(i,j,k,nrhs)=rv(i,j,k,nrhs)-cff
            END DO
          END DO
        END IF
!
!-----------------------------------------------------------------------
!  Compute forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. If not
!  body force stresses, add in the difference between surface and
!  bottom stresses.
!
        DO i=IstrU,Iend
          rufrc(i,j)=ru(i,j,1,nrhs)
        END DO
        DO k=2,N(ng)
          DO i=IstrU,Iend
            rufrc(i,j)=rufrc(i,j)+ru(i,j,k,nrhs)
          END DO
        END DO
        DO i=IstrU,Iend
          cff=om_u(i,j)*on_u(i,j)
          cff1= sustr(i,j)*cff
          cff2=-bustr(i,j)*cff
          rufrc(i,j)=rufrc(i,j)+cff1+cff2
        END DO
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            rvfrc(i,j)=rv(i,j,1,nrhs)
          END DO
          DO k=2,N(ng)
            DO i=Istr,Iend
              rvfrc(i,j)=rvfrc(i,j)+rv(i,j,k,nrhs)
            END DO
          END DO
          DO i=Istr,Iend
            cff=om_v(i,j)*on_v(i,j)
            cff1= svstr(i,j)*cff
            cff2=-bvstr(i,j)*cff
            rvfrc(i,j)=rvfrc(i,j)+cff1+cff2
          END DO
        END IF
      END DO J_LOOP
      RETURN
      END SUBROUTINE rhs3d_tile
      END MODULE rhs3d_mod

      MODULE gls_prestep_mod
!
!svn $Id: gls_prestep.F 709 2014-01-23 20:09:38Z arango $
!=======================================================================
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.txt                   Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This routine perfoms the predictor step for turbulent kinetic       !
!  energy prognostic variables, tke and gls. A NON-conservative,       !
!  but constancy preserving, auxiliary advective substep for tke       !
!  gls equations is carried out. The result of this substep will       !
!  be used to compute advective terms in the corrector substep.        !
!  No dissipation terms are included here.                             !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: gls_prestep
      CONTAINS
!
!***********************************************************************
      SUBROUTINE gls_prestep (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_mixing
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
      CALL gls_prestep_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nstp(ng), nnew(ng),                        &
     &                       GRID(ng) % Huon,                           &
     &                       GRID(ng) % Hvom,                           &
     &                       GRID(ng) % Hz,                             &
     &                       GRID(ng) % pm,                             &
     &                       GRID(ng) % pn,                             &
     &                       OCEAN(ng) % W,                             &
     &                       MIXING(ng) % gls,                          &
     &                       MIXING(ng) % tke)
      CALL wclock_off (ng, iNLM, 19)
      RETURN
      END SUBROUTINE gls_prestep
!
!***********************************************************************
      SUBROUTINE gls_prestep_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             nstp, nnew,                          &
     &                             Huon, Hvom, Hz, pm, pn, W,           &
     &                             gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_ncparam
      USE mod_scalars
!
      USE exchange_3d_mod, ONLY : exchange_w3d_tile
      USE mp_exchange_mod, ONLY : mp_exchange3d
      USE tkebc_mod, ONLY : tkebc_tile
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew
!
      real(r8), intent(in) :: Huon(LBi:,LBj:,:)
      real(r8), intent(in) :: Hvom(LBi:,LBj:,:)
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: pm(LBi:,LBj:)
      real(r8), intent(in) :: pn(LBi:,LBj:)
      real(r8), intent(in) :: W(LBi:,LBj:,0:)
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      integer :: i, indx, j, k
      real(r8), parameter :: Gamma = 1.0_r8/6.0_r8
      real(r8) :: cff, cff1, cff2, cff3, cff4
      real(r8), dimension(IminS:ImaxS,N(ng)) :: CF
      real(r8), dimension(IminS:ImaxS,N(ng)) :: FC
      real(r8), dimension(IminS:ImaxS,N(ng)) :: FCL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng)) :: Hz_half
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: EF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FE
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FEL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FX
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: FXL
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: XF
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: grad
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: gradL
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
!  Predictor step for advection of turbulent kinetic energy variables.
!-----------------------------------------------------------------------
!
! Start computation of auxiliary time step fields tke(:,:,:,n+1/2) and
! gls(:,:,:,n+1/2) with computation of horizontal advection terms and
! auxiliary grid-box height field Hz_new()=Hz(:,:,k+1/2,n+1/2);
! This is effectivey an LF step with subsequent interpolation of the
! result half step back, using AM3 weights. The LF step and
! interpolation are perfomed as a single operation, which results in
! weights cff1,cff2,cff3 below.
!
! Either centered fourth-order accurate or standard second order
! accurate versions are supported.
!
! At the same time prepare for corrector step for tke,gls: set tke,
! gls(:,:,:,nnew) to  tke,gls(:,:,:,nstp) multiplied by the
! corresponding grid-box height. This needs done at this time because
! array Hz(:,:,:) will overwritten after 2D time stepping with the
! values computed from zeta(:,:,n+1) rather than zeta(:,:,n), so that
! the old-time-step Hz will be no longer awailable.
!
      DO k=1,N(ng)-1
!
!  Fourth-order, centered differences advection.
!
        DO j=Jstr,Jend
          DO i=Istrm1,Iendp2
            grad (i,j)=(tke(i,j,k,nstp)-tke(i-1,j,k,nstp))
            gradL(i,j)=(gls(i,j,k,nstp)-gls(i-1,j,k,nstp))
          END DO
        END DO
        IF (.not.(CompositeGrid(iwest,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Western_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Istr-1,j)=grad (Istr,j)
              gradL(Istr-1,j)=gradL(Istr,j)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(ieast,ng).or.EWperiodic(ng))) THEN
          IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
            DO j=Jstr,Jend
              grad (Iend+2,j)=grad (Iend+1,j)
              gradL(Iend+2,j)=gradL(Iend+1,j)
            END DO
          END IF
        END IF
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend
          DO i=Istr,Iend+1
            XF(i,j)=0.5_r8*(Huon(i,j,k)+Huon(i,j,k+1))
            FX (i,j)=XF(i,j)*                                           &
     &               0.5_r8*(tke(i-1,j,k,nstp)+tke(i,j,k,nstp)-         &
     &                       cff*(grad (i+1,j)-grad (i-1,j)))
            FXL(i,j)=XF(i,j)*                                           &
     &               0.5_r8*(gls(i-1,j,k,nstp)+gls(i,j,k,nstp)-         &
     &                       cff*(gradL(i+1,j)-gradL(i-1,j)))
          END DO
        END DO
!
        DO j=Jstrm1,Jendp2
          DO i=Istr,Iend
            grad (i,j)=(tke(i,j,k,nstp)-tke(i,j-1,k,nstp))
            gradL(i,j)=(gls(i,j,k,nstp)-gls(i,j-1,k,nstp))
          END DO
        END DO
        IF (.not.(CompositeGrid(isouth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jstr-1)=grad (i,Jstr)
              gradL(i,Jstr-1)=gradL(i,Jstr)
            END DO
          END IF
        END IF
        IF (.not.(CompositeGrid(inorth,ng).or.NSperiodic(ng))) THEN
          IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
            DO i=Istr,Iend
              grad (i,Jend+2)=grad (i,Jend+1)
              gradL(i,Jend+2)=gradL(i,Jend+1)
            END DO
          END IF
        END IF
        cff=1.0_r8/6.0_r8
        DO j=Jstr,Jend+1
          DO i=Istr,Iend
            EF(i,j)=0.5_r8*(Hvom(i,j,k)+Hvom(i,j,k+1))
            FE (i,j)=EF(i,j)*                                           &
     &               0.5_r8*(tke(i,j-1,k,nstp)+tke(i,j,k,nstp)-         &
     &                       cff*(grad (i,j+1)-grad (i,j-1)))
            FEL(i,j)=EF(i,j)*                                           &
     &               0.5_r8*(gls(i,j-1,k,nstp)+gls(i,j,k,nstp)-         &
     &                       cff*(gradL(i,j+1)-gradL(i,j-1)))
          END DO
        END DO
!
!  Time-step horizontal advection.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff1=1.0_r8
          cff2=0.0_r8
          cff3=0.5_r8*dt(ng)
          indx=nstp
        ELSE
          cff1=0.5_r8+Gamma
          cff2=0.5_r8-Gamma
          cff3=(1.0_r8-Gamma)*dt(ng)
          indx=3-nstp
        END IF
        DO j=Jstr,Jend
          DO i=Istr,Iend
            cff=0.5_r8*(Hz(i,j,k)+Hz(i,j,k+1))
            cff4=cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k)=cff-cff4*(XF(i+1,j)-XF(i,j)+                 &
     &                               EF(i,j+1)-EF(i,j))
            tke(i,j,k,3)=cff*(cff1*tke(i,j,k,nstp)+                     &
     &                        cff2*tke(i,j,k,indx))-                    &
     &                   cff4*(FX (i+1,j)-FX (i,j)+                     &
     &                         FE (i,j+1)-FE (i,j))
            gls(i,j,k,3)=cff*(cff1*gls(i,j,k,nstp)+                     &
     &                        cff2*gls(i,j,k,indx))-                    &
     &                   cff4*(FXL(i+1,j)-FXL(i,j)+                     &
     &                         FEL(i,j+1)-FEL(i,j))
            tke(i,j,k,nnew)=cff*tke(i,j,k,nstp)
            gls(i,j,k,nnew)=cff*gls(i,j,k,nstp)
          END DO
        END DO
      END DO
!
! Compute vertical advection term.
!
      DO j=Jstr,Jend
        cff1=7.0_r8/12.0_r8
        cff2=1.0_r8/12.0_r8
        DO k=2,N(ng)-1
          DO i=Istr,Iend
            CF(i,k)=0.5_r8*(W(i,j,k)+W(i,j,k-1))
            FC (i,k)=CF(i,k)*(cff1*(tke(i,j,k-1,nstp)+                  &
     &                              tke(i,j,k  ,nstp))-                 &
     &                        cff2*(tke(i,j,k-2,nstp)+                  &
     &                              tke(i,j,k+1,nstp)))
            FCL(i,k)=CF(i,k)*(cff1*(gls(i,j,k-1,nstp)+                  &
     &                              gls(i,j,k  ,nstp))-                 &
     &                        cff2*(gls(i,j,k-2,nstp)+                  &
     &                              gls(i,j,k+1,nstp)))
          END DO
        END DO
        cff1=1.0_r8/3.0_r8
        cff2=5.0_r8/6.0_r8
        cff3=1.0_r8/6.0_r8
        DO i=Istr,Iend
          CF(i,1)=0.5*(W(i,j,0)+W(i,j,1))
          FC (i,1)=CF(i,1)*(cff1*tke(i,j,0,nstp)+                       &
     &                      cff2*tke(i,j,1,nstp)-                       &
     &                      cff3*tke(i,j,2,nstp))
          FCL(i,1)=CF(i,1)*(cff1*gls(i,j,0,nstp)+                       &
     &                      cff2*gls(i,j,1,nstp)-                       &
     &                      cff3*gls(i,j,2,nstp))
          CF(i,N(ng))=0.5*(W(i,j,N(ng))+W(i,j,N(ng)-1))
          FC (i,N(ng))=CF(i,N(ng))*(cff1*tke(i,j,N(ng)  ,nstp)+         &
     &                              cff2*tke(i,j,N(ng)-1,nstp)-         &
     &                              cff3*tke(i,j,N(ng)-2,nstp))
          FCL(i,N(ng))=CF(i,N(ng))*(cff1*gls(i,j,N(ng)  ,nstp)+         &
     &                              cff2*gls(i,j,N(ng)-1,nstp)-         &
     &                              cff3*gls(i,j,N(ng)-2,nstp))
        END DO
!
!  Time-step vertical advection term.
!
        IF (iic(ng).eq.ntfirst(ng)) THEN
          cff3=0.5_r8*dt(ng)
        ELSE
          cff3=(1.0_r8-Gamma)*dt(ng)
        END IF
        DO k=1,N(ng)-1
          DO i=Istr,Iend
            cff4=cff3*pm(i,j)*pn(i,j)
            Hz_half(i,j,k)=Hz_half(i,j,k)-cff4*(CF(i,k+1)-CF(i,k))
            cff1=1.0_r8/Hz_half(i,j,k)
            tke(i,j,k,3)=cff1*(tke(i,j,k,3)-                            &
     &                         cff4*(FC (i,k+1)-FC (i,k)))
            gls(i,j,k,3)=cff1*(gls(i,j,k,3)-                            &
     &                         cff4*(FCL(i,k+1)-FCL(i,k)))
          END DO
        END DO
      END DO
!
!  Apply lateral boundary conditions.
!
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 3, nstp,                                         &
     &                 gls, tke)
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          tke(:,:,:,3))
        CALL exchange_w3d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj, 0, N(ng),           &
     &                          gls(:,:,:,3))
      END IF
      CALL mp_exchange3d (ng, tile, iNLM, 2,                            &
     &                    LBi, UBi, LBj, UBj, 0, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tke(:,:,:,3),                                 &
     &                    gls(:,:,:,3))
      RETURN
      END SUBROUTINE gls_prestep_tile
      END MODULE gls_prestep_mod

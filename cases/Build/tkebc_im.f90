      MODULE tkebc_mod
!
!svn $Id: tkebc_im.F 709 2014-01-23 20:09:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine sets lateral boundary conditions for turbulent      !
!  kinetic energy and turbulent length scale variables associated      !
!  with the Mellor and Yamada or GOTM closures.                        !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: tkebc_tile
      CONTAINS
!
!***********************************************************************
      SUBROUTINE tkebc (ng, tile, nout)
!***********************************************************************
!
      USE mod_param
      USE mod_mixing
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, nout
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
      CALL tkebc_tile (ng, tile,                                        &
     &                 LBi, UBi, LBj, UBj, N(ng),                       &
     &                 IminS, ImaxS, JminS, JmaxS,                      &
     &                 nout, nstp(ng),                                  &
     &                 MIXING(ng)% gls,                                 &
     &                 MIXING(ng)% tke)
      RETURN
      END SUBROUTINE tkebc
!
!***********************************************************************
      SUBROUTINE tkebc_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj, UBk,                   &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       nout, nstp,                                &
     &                       gls, tke)
!***********************************************************************
!
      USE mod_param
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nout, nstp
!
      real(r8), intent(inout) :: gls(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: tke(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8), parameter :: eps = 1.0e-20_r8
      real(r8) :: Ce, Cx, cff, dKde, dKdt, dKdx
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
!  Lateral boundary conditions at the western edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!
!  Western edge, implicit upstream radiation condition.
!
        IF (LBC(iwest,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend+1
              grad(Istr-1,j)=tke(Istr-1,j  ,k,nstp)-                    &
     &                       tke(Istr-1,j-1,k,nstp)
              grad(Istr  ,j)=tke(Istr  ,j  ,k,nstp)-                    &
     &                       tke(Istr  ,j-1,k,nstp)
              gradL(Istr-1,j)=gls(Istr-1,j  ,k,nstp)-                   &
     &                        gls(Istr-1,j-1,k,nstp)
              gradL(Istr  ,j)=gls(Istr  ,j  ,k,nstp)-                   &
     &                        gls(Istr  ,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                dKdt=tke(Istr,j,k,nstp)-tke(Istr  ,j,k,nout)
                dKdx=tke(Istr,j,k,nout)-tke(Istr+1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(Istr,j  )+                              &
     &                     grad(Istr,j+1))).gt.0.0_r8) THEN
                  dKde=grad(Istr,j  )
                ELSE
                  dKde=grad(Istr,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
                tke(Istr-1,j,k,nout)=(cff*tke(Istr-1,j,k,nstp)+         &
     &                                Cx *tke(Istr  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   grad(Istr-1,j  )-              &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   grad(Istr-1,j+1))/             &
     &                               (cff+Cx)
                dKdt=gls(Istr,j,k,nstp)-gls(Istr  ,j,k,nout)
                dKdx=gls(Istr,j,k,nout)-gls(Istr+1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(Istr,j  )+                             &
     &                     gradL(Istr,j+1))).gt.0.0_r8) THEN
                  dKde=gradL(Istr,j  )
                ELSE
                  dKde=gradL(Istr,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
                gls(Istr-1,j,k,nout)=(cff*gls(Istr-1,j,k,nstp)+         &
     &                                Cx *gls(Istr  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   gradL(Istr-1,j  )-             &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   gradL(Istr-1,j+1))/            &
     &                               (cff+Cx)
              END IF
            END DO
          END DO
!
!  Western edge, gradient boundary condition.
!
        ELSE IF (LBC(iwest,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
                gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
              END IF
            END DO
          END DO
!
!  Western edge, closed boundary condition.
!
        ELSE IF (LBC(iwest,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%west(j)) THEN
                tke(Istr-1,j,k,nout)=tke(Istr,j,k,nout)
                gls(Istr-1,j,k,nout)=gls(Istr,j,k,nout)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the eastern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!
!  Eastern edge, implicit upstream radiation condition.
!
        IF (LBC(ieast,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend+1
              grad(Iend  ,j)=tke(Iend  ,j  ,k,nstp)-                    &
     &                       tke(Iend  ,j-1,k,nstp)
              grad(Iend+1,j)=tke(Iend+1,j  ,k,nstp)-                    &
     &                       tke(Iend+1,j-1,k,nstp)
              gradL(Iend  ,j)=gls(Iend  ,j  ,k,nstp)-                   &
     &                        gls(Iend  ,j-1,k,nstp)
              gradL(Iend+1,j)=gls(Iend+1,j  ,k,nstp)-                   &
     &                        gls(Iend+1,j-1,k,nstp)
            END DO
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                dKdt=tke(Iend,j,k,nstp)-tke(Iend  ,j,k,nout)
                dKdx=tke(Iend,j,k,nout)-tke(Iend-1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(Iend,j  )+                              &
     &                     grad(Iend,j+1))).gt.0.0_r8) THEN
                  dKde=grad(Iend,j  )
                ELSE
                  dKde=grad(Iend,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
                tke(Iend+1,j,k,nout)=(cff*tke(Iend+1,j,k,nstp)+         &
     &                                Cx *tke(Iend  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   grad(Iend+1,j  )-              &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   grad(Iend+1,j+1))/             &
     &                               (cff+Cx)
                dKdt=gls(Iend,j,k,nstp)-gls(Iend  ,j,k,nout)
                dKdx=gls(Iend,j,k,nout)-gls(Iend-1,j,k,nout)
                IF ((dKdt*dKdx).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(Iend,j  )+                             &
     &                     gradL(Iend,j+1))).gt.0.0_r8) THEN
                  dKde=gradL(Iend,j  )
                ELSE
                  dKde=gradL(Iend,j+1)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=dKdt*dKdx
                Ce=MIN(cff,MAX(dKdt*dKde,-cff))
                gls(Iend+1,j,k,nout)=(cff*gls(Iend+1,j,k,nstp)+         &
     &                                Cx *gls(Iend  ,j,k,nout)-         &
     &                                MAX(Ce,0.0_r8)*                   &
     &                                   gradL(Iend+1,j  )-             &
     &                                MIN(Ce,0.0_r8)*                   &
     &                                   gradL(Iend+1,j+1))/            &
     &                               (cff+Cx)
              END IF
            END DO
          END DO
!
!  Eastern edge, gradient boundary condition.
!
        ELSE IF (LBC(ieast,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
                gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
              END IF
            END DO
          END DO
!
!  Eastern edge, closed boundary condition.
!
        ELSE IF (LBC(ieast,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO j=Jstr,Jend
              IF (LBC_apply(ng)%east(j)) THEN
                tke(Iend+1,j,k,nout)=tke(Iend,j,k,nout)
                gls(Iend+1,j,k,nout)=gls(Iend,j,k,nout)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the southern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!
!  Southern edge, implicit upstream radiation condition.
!
        IF (LBC(isouth,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jstr  )=tke(i  ,Jstr  ,k,nstp)-                    &
     &                       tke(i-1,Jstr  ,k,nstp)
              grad(i,Jstr-1)=tke(i  ,Jstr-1,k,nstp)-                    &
     &                     tke(i-1,Jstr-1,k,nstp)
              gradL(i,Jstr  )=gls(i  ,Jstr  ,k,nstp)-                   &
     &                      gls(i-1,Jstr  ,k,nstp)
              gradL(i,Jstr-1)=gls(i  ,Jstr-1,k,nstp)-                   &
     &                      gls(i-1,Jstr-1,k,nstp)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                dKdt=tke(i,Jstr,k,nstp)-tke(i,Jstr  ,k,nout)
                dKde=tke(i,Jstr,k,nout)-tke(i,Jstr+1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(i  ,Jstr)+                              &
     &                     grad(i+1,Jstr))).gt.0.0_r8) THEN
                  dKdx=grad(i  ,Jstr)
                ELSE
                  dKdx=grad(i+1,Jstr)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde, eps)
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
                Ce=dKdt*dKde
                tke(i,Jstr-1,k,nout)=(cff*tke(i,Jstr-1,k,nstp)+         &
     &                                Ce *tke(i,Jstr  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   grad(i  ,Jstr-1)-              &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   grad(i+1,Jstr-1))/             &
     &                               (cff+Ce)
                dKdt=gls(i,Jstr,k,nstp)-gls(i,Jstr  ,k,nout)
                dKde=gls(i,Jstr,k,nout)-gls(i,Jstr+1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(i  ,Jstr)+                             &
     &                     gradL(i+1,Jstr))).gt.0.0_r8) THEN
                  dKdx=gradL(i  ,Jstr)
                ELSE
                  dKdx=gradL(i+1,Jstr)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
                Ce=dKdt*dKde
                gls(i,Jstr-1,k,nout)=(cff*gls(i,Jstr-1,k,nstp)+         &
     &                                Ce *gls(i,Jstr  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   gradL(i  ,Jstr-1)-             &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   gradL(i+1,Jstr-1))/            &
     &                               (cff+Ce)
              END IF
            END DO
          END DO
!
!  Southern edge, gradient boundary condition.
!
        ELSE IF (LBC(isouth,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
                gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
              END IF
            END DO
          END DO
!
!  Southern edge, closed boundary condition.
!
        ELSE IF (LBC(isouth,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%south(i)) THEN
                tke(i,Jstr-1,k,nout)=tke(i,Jstr,k,nout)
                gls(i,Jstr-1,k,nout)=gls(i,Jstr,k,nout)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Lateral boundary conditions at the northern edge.
!-----------------------------------------------------------------------
!
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!
!  Northern edge, implicit upstream radiation condition.
!
        IF (LBC(inorth,isMtke,ng)%radiation) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend+1
              grad(i,Jend  )=tke(i  ,Jend  ,k,nstp)-                    &
     &                       tke(i-1,Jend  ,k,nstp)
              grad(i,Jend+1)=tke(i  ,Jend+1,k,nstp)-                    &
     &                       tke(i-1,Jend+1,k,nstp)
              gradL(i,Jend  )=gls(i  ,Jend  ,k,nstp)-                   &
     &                        gls(i-1,Jend  ,k,nstp)
              gradL(i,Jend+1)=gls(i  ,Jend+1,k,nstp)-                   &
     &                        gls(i-1,Jend+1,k,nstp)
            END DO
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                dKdt=tke(i,Jend,k,nstp)-tke(i,Jend  ,k,nout)
                dKde=tke(i,Jend,k,nout)-tke(i,Jend-1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(grad(i  ,Jend)+                              &
     &                     grad(i+1,Jend))).gt.0.0_r8) THEN
                  dKdx=grad(i  ,Jend)
                ELSE
                  dKdx=grad(i+1,Jend)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
                Ce=dKdt*dKde
                tke(i,Jend+1,k,nout)=(cff*tke(i,Jend+1,k,nstp)+         &
     &                                Ce *tke(i,Jend  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   grad(i  ,Jend+1)-              &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   grad(i+1,Jend+1))/             &
     &                               (cff+Ce)
                dKdt=gls(i,Jend,k,nstp)-gls(i,Jend  ,k,nout)
                dKde=gls(i,Jend,k,nout)-gls(i,Jend-1,k,nout)
                IF ((dKdt*dKde).lt.0.0_r8) dKdt=0.0_r8
                IF ((dKdt*(gradL(i  ,Jend)+                             &
     &                     gradL(i+1,Jend))).gt.0.0_r8) THEN
                  dKdx=gradL(i  ,Jend)
                ELSE
                  dKdx=gradL(i+1,Jend)
                END IF
                cff=MAX(dKdx*dKdx+dKde*dKde,eps)
                Cx=MIN(cff,MAX(dKdt*dKdx,-cff))
                Ce=dKdt*dKde
                gls(i,Jend+1,k,nout)=(cff*gls(i,Jend+1,k,nstp)+         &
     &                                Ce *gls(i,Jend  ,k,nout)-         &
     &                                MAX(Cx,0.0_r8)*                   &
     &                                   gradL(i  ,Jend+1)-             &
     &                                MIN(Cx,0.0_r8)*                   &
     &                                   gradL(i+1,Jend+1))/            &
     &                               (cff+Ce)
              END IF
            END DO
          END DO
!
!  Northern edge, gradient boundary condition.
!
        ELSE IF (LBC(inorth,isMtke,ng)%gradient) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
                gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
              END IF
            END DO
          END DO
!
!  Northern edge, closed boundary condition.
!
        ELSE IF (LBC(inorth,isMtke,ng)%closed) THEN
          DO k=0,N(ng)
            DO i=Istr,Iend
              IF (LBC_apply(ng)%north(i)) THEN
                tke(i,Jend+1,k,nout)=tke(i,Jend,k,nout)
                gls(i,Jend+1,k,nout)=gls(i,Jend,k,nout)
              END IF
            END DO
          END DO
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Boundary corners.
!-----------------------------------------------------------------------
!
      IF (.not.(EWperiodic(ng).or.NSperiodic(ng))) THEN
        IF (DOMAIN(ng)%SouthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jstr-1)) THEN
            DO k=0,N(ng)
              tke(Istr-1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (tke(Istr  ,Jstr-1,k,nout)+     &
     &                                   tke(Istr-1,Jstr  ,k,nout))
              gls(Istr-1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (gls(Istr  ,Jstr-1,k,nout)+     &
     &                                   gls(Istr-1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%SouthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%south(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jstr-1)) THEN
            DO k=0,N(ng)
              tke(Iend+1,Jstr-1,k,nout)=0.5_r8*                         &
     &                                  (tke(Iend  ,Jstr-1,k,nout)+     &
     &                                   tke(Iend+1,Jstr  ,k,nout))
               gls(Iend+1,Jstr-1,k,nout)=0.5_r8*                        &
     &                                   (gls(Iend  ,Jstr-1,k,nout)+    &
     &                                    gls(Iend+1,Jstr  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthWest_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Istr-1).and.                          &
     &        LBC_apply(ng)%west (Jend+1)) THEN
            DO k=0,N(ng)
              tke(Istr-1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (tke(Istr  ,Jend+1,k,nout)+     &
     &                                   tke(Istr-1,Jend  ,k,nout))
              gls(Istr-1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (gls(Istr  ,Jend+1,k,nout)+     &
     &                                   gls(Istr-1,Jend  ,k,nout))
            END DO
          END IF
        END IF
        IF (DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          IF (LBC_apply(ng)%north(Iend+1).and.                          &
     &        LBC_apply(ng)%east (Jend+1)) THEN
            DO k=0,N(ng)
              tke(Iend+1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (tke(Iend  ,Jend+1,k,nout)+     &
     &                                   tke(Iend+1,Jend  ,k,nout))
              gls(Iend+1,Jend+1,k,nout)=0.5_r8*                         &
     &                                  (gls(Iend  ,Jend+1,k,nout)+     &
     &                                   gls(Iend+1,Jend  ,k,nout))
            END DO
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE tkebc_tile
      END MODULE tkebc_mod

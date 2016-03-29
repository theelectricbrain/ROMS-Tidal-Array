      MODULE prsgrd_mod
!
!svn $Id: prsgrd.F 709 2014-01-23 20:09:38Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the baroclinic hydrostatic pressure gradient  !
!  term.                                                               !
!                                                                      !
!=======================================================================
!
      implicit none
      PRIVATE
      PUBLIC  :: prsgrd
      CONTAINS
      SUBROUTINE prsgrd (ng, tile)
!
!svn $Id: prsgrd31.h 709 2014-01-23 20:09:38Z arango $
!***********************************************************************
!  Copyright (c) 2002-2014 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                           Hernan G. Arango   !
!****************************************** Alexander F. Shchepetkin ***
!                                                                      !
!  This subroutine evaluates the  baroclinic  hydrostatic  pressure    !
!  gradient term using  the STANDARD density Jacobian  or  WEIGHTED    !
!  density Jacobian scheme of Song (1998). Both of these approaches    !
!  compute horizontal differences of density before of the vertical    !
!  integration.                                                        !
!                                                                      !
!  The pressure gradient terms (m4/s2) are loaded into right-hand-     !
!  side arrays "ru" and "rv".                                          !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Song, Y.T., 1998:  A general pressure gradient formulation for    !
!      numerical ocean models. Part I: Scheme design and diagnostic    !
!      analysis, Monthly Weather Rev., 126, 3213-3230.                 !
!                                                                      !
!***********************************************************************
!
      USE mod_param
      USE mod_grid
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
      CALL wclock_on (ng, iNLM, 23)
      CALL prsgrd_tile (ng, tile,                                       &
     &                  LBi, UBi, LBj, UBj,                             &
     &                  IminS, ImaxS, JminS, JmaxS,                     &
     &                  nrhs(ng),                                       &
     &                  GRID(ng) % Hz,                                  &
     &                  GRID(ng) % om_v,                                &
     &                  GRID(ng) % on_u,                                &
     &                  GRID(ng) % z_r,                                 &
     &                  GRID(ng) % z_w,                                 &
     &                  OCEAN(ng) % rho,                                &
     &                  OCEAN(ng) % ru,                                 &
     &                  OCEAN(ng) % rv)
      CALL wclock_off (ng, iNLM, 23)
      RETURN
      END SUBROUTINE prsgrd
!
!***********************************************************************
      SUBROUTINE prsgrd_tile (ng, tile,                                 &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS,               &
     &                        nrhs,                                     &
     &                        Hz, om_v, on_u, z_r, z_w,                 &
     &                        rho,                                      &
     &                        ru, rv)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nrhs
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
      real(r8), intent(in) :: om_v(LBi:,LBj:)
      real(r8), intent(in) :: on_u(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(in) :: z_w(LBi:,LBj:,0:)
      real(r8), intent(in) :: rho(LBi:,LBj:,:)
      real(r8), intent(inout) :: ru(LBi:,LBj:,0:,:)
      real(r8), intent(inout) :: rv(LBi:,LBj:,0:,:)
!
!  Local variable declarations.
!
      integer :: i, j, k
      real(r8) :: fac, fac1, fac2, fac3
      real(r8) :: cff1, cff2, cff3, cff4
      real(r8), dimension(IminS:ImaxS) :: phie
      real(r8), dimension(IminS:ImaxS) :: phix
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
!  Calculate pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
      fac1=0.5_r8*g/rho0
      fac2=1000.0_r8*g/rho0
      fac3=0.25_r8*g/rho0
      DO j=Jstr,Jend
        DO i=IstrU,Iend
          cff1=z_w(i  ,j,N(ng))-z_r(i  ,j,N(ng))+                       &
     &         z_w(i-1,j,N(ng))-z_r(i-1,j,N(ng))
          phix(i)=fac1*(rho(i,j,N(ng))-rho(i-1,j,N(ng)))*cff1
          phix(i)=phix(i)+                                              &
     &            (fac2+fac1*(rho(i,j,N(ng))+rho(i-1,j,N(ng))))*        &
     &            (z_w(i,j,N(ng))-z_w(i-1,j,N(ng)))
          ru(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i-1,j,N(ng)))*   &
     &                       phix(i)*on_u(i,j)
        END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
        DO k=N(ng)-1,1,-1
          DO i=IstrU,Iend
            cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-                           &
     &           rho(i,j,k  )-rho(i-1,j,k  )
            cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+                           &
     &           z_r(i,j,k  )-z_r(i-1,j,k  )
            phix(i)=phix(i)+                                            &
     &              fac3*(cff1*cff3-cff2*cff4)
            ru(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i-1,j,k))*             &
     &                     phix(i)*on_u(i,j)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient.
!
        IF (j.ge.JstrV) THEN
          DO i=Istr,Iend
            cff1=z_w(i,j  ,N(ng))-z_r(i,j  ,N(ng))+                     &
     &           z_w(i,j-1,N(ng))-z_r(i,j-1,N(ng))
            phie(i)=fac1*(rho(i,j,N(ng))-rho(i,j-1,N(ng)))*cff1
            phie(i)=phie(i)+                                            &
     &              (fac2+fac1*(rho(i,j,N(ng))+rho(i,j-1,N(ng))))*      &
     &              (z_w(i,j,N(ng))-z_w(i,j-1,N(ng)))
            rv(i,j,N(ng),nrhs)=-0.5_r8*(Hz(i,j,N(ng))+Hz(i,j-1,N(ng)))* &
     &                         phie(i)*om_v(i,j)
          END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate.
!
          DO k=N(ng)-1,1,-1
            DO i=Istr,Iend
              cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-                         &
     &             rho(i,j,k  )-rho(i,j-1,k  )
              cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+                         &
     &             z_r(i,j,k  )-z_r(i,j-1,k  )
              phie(i)=phie(i)+                                          &
     &                fac3*(cff1*cff3-cff2*cff4)
              rv(i,j,k,nrhs)=-0.5_r8*(Hz(i,j,k)+Hz(i,j-1,k))*           &
     &                       phie(i)*om_v(i,j)
            END DO
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE prsgrd_tile
      END MODULE prsgrd_mod

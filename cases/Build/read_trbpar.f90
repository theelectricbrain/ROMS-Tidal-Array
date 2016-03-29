      SUBROUTINE read_TrbPar (model, inp, out, Lwrite)
!
!  read_trbpar.F 2013-01-17 13:40:00
!================================================== Thomas Roc==========
!                                                                      !
!  Copyright (c) 2014 Thomas Roc and ITPower                           !
!    Licensed under an Affero GPL style license                        !
!    See License_ROMS-Tidal.txt                                        !
!=======================================================================
!                                                                      !
!  This routine reads and reports turbine input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      logical :: find_file
      integer :: Npts, Nval
      integer :: C, i, j, igrid, ng, status
      integer :: decode_line, load_i, load_l, load_r
      real(r8) :: xpos, ypos, zpos
      real(r8) :: Ct, Ctke, Cgls
      real(r8) :: Lc, Pa, Diam
      integer, dimension(Ngrids) :: is
      real(r8), dimension(100) :: Rval
      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval
!
!-----------------------------------------------------------------------
!  Read turbine parameters.
!-----------------------------------------------------------------------
!
!  Allocate turbines parameters 
!
      DO WHILE (.TRUE.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          SELECT CASE (TRIM(KeyWord))
            CASE ('Lturbines')
              Npts=load_l(Nval, Cval, Ngrids, Lturbines)
            CASE ('NTURBINES')
              Npts=load_i(Nval, Rval, Ngrids, Nturbines)
! For allocation alternative using allocate_turbines
            CASE ('POS')
              DO ng=1,Ngrids
                allocate ( SCALARS(ng) % Tc(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Txpos(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Typos(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tzpos(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tct(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Ttke(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tgls(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tlc(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tpa(Nturbines(ng)) )
                allocate ( SCALARS(ng) % Tdiam(Nturbines(ng)) )
              END DO
              is(1:Ngrids)=0
              DO WHILE (.TRUE.)
                READ (inp,*,ERR=10,END=10) igrid, C, xpos, ypos, zpos,  &
     &                                     Ct, Ctke, Cgls, Lc, Pa, Diam
                ng=MAX(1,ABS(igrid))
                is(ng)=is(ng)+1
                SCALARS(ng)%Tc(is(ng))=C
                SCALARS(ng)%Txpos(is(ng))=xpos
                SCALARS(ng)%Typos(is(ng))=ypos
                SCALARS(ng)%Tzpos(is(ng))=zpos
                SCALARS(ng)%Tct(is(ng))=Ct
                SCALARS(ng)%Ttke(is(ng))=Ctke
                SCALARS(ng)%Tgls(is(ng))=Cgls
                SCALARS(ng)%Tlc(is(ng))=Lc
                SCALARS(ng)%Tpa(is(ng))=Pa
                SCALARS(ng)%Tdiam(is(ng))=Diam
              END DO
 10           DO ng=1,Ngrids
                IF (Nturbines(ng).ne.is(ng)) THEN
                  IF (Master) WRITE (*,40) Nturbines(ng), is(ng)
                  exit_flag=4
                  RETURN
                END IF
              END DO
          END SELECT
        END IF
      END DO
  20  IF (Master) WRITE (*,50) line
      exit_flag=4
      RETURN
  30  CONTINUE
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lturbines(ng)) THEN
            WRITE (out,70) Nturbines(ng), 'Number of turbines deployed.'
            WRITE (out,60) ng
            WRITE (out,*)
            DO i=1,Nturbines(ng)
              WRITE (out,100) 'No', i,                                  &
     &                           SCALARS(ng)%Txpos(i),                  &
     &                           SCALARS(ng)%Typos(i),                  &
     &                           SCALARS(ng)%Tzpos(i),                  &
     &                           SCALARS(ng)%Tct(i),                    &
     &                           SCALARS(ng)%Ttke(i),                   &
     &                           SCALARS(ng)%Tgls(i),                   &
     &                           SCALARS(ng)%Tlc(i),                    &
     &                           SCALARS(ng)%Tpa(i),                    &
     &                           SCALARS(ng)%Tdiam(i)
            END DO
          END IF
        END DO
      END IF
  40  FORMAT (/,' READ_TrbPar - Inconsistent number of turbines, ',     &
     &        'Nturbines = ',2i8,/,15x,'change input script values.')
  50  FORMAT (/,' READ_TrbPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Turbines Parameters, Grid: ',i2.2,                  &
     &        /,  ' =============================',/,/,                 &
     &        14x,'xpos',5x,'ypos',4x,'zpos',3x,'Ct',                   &
     &        3x,'Ctke',3x,'Cgls',3x,'Lc',3x,'Pa',3x,'Diam',/)
  70  FORMAT (/,1x,i3,2x,a,t30,a)
 100  FORMAT (a,i3," :",1x,d9.3,1x,d9.3,1x,f5.2,1x,f5.2,1x,f5.2,1x,f5.2,&
     &        1x,f5.2,1x,f5.2,1x,f5.2)
 110  FORMAT (/,' READ_TrbPAR - variable info not yet loaded, ', a)
      RETURN
      END SUBROUTINE read_TrbPar

MODULE usrdef_istate
   !!==============================================================================
   !!                       ***  MODULE usrdef_istate  ***
   !!
   !!                      ===  OVERFLOW configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!==============================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce , ONLY : glamt
   USE usrdef_nam , ONLY : rn_T0, rn_T1
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10074 2018-08-28 16:15:49Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here OVERFLOW configuration
      !!
      !! ** Method  : - set temprature field
      !!              - set salinity   field
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s]
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  ::   jk,jj,ji     ! dummy loop indices
      REAL(wp) ::   zdam   ! location of dam [Km]
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : OVERFLOW configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   Ocean at rest, with a constant salinity (not used as rho=F(T) '
      IF(lwp) WRITE(numout,*) '                 and a vertical density front with a 2 kg/m3 difference located at glam=20km'
      IF(lwp) WRITE(numout,*) '                 (i.e. a temperature difference of 10 degrees with rn_a0 = 0.2'
      !
      !  rn_a0 =  0.2   !  thermal expension coefficient (nn_eos= 1)
      !  rho = rau0 - rn_a0 * (T-10)
      !  delta_T = 10 degrees  ==>>  delta_rho = 10 * rn_a0 = 2 kg/m3
      !
      pu  (:,:,:) = 0._wp        ! ocean at rest
      pv  (:,:,:) = 0._wp
      pssh(:,:)   = 0._wp
      !
      !                          ! T & S profiles
      zdam = 20._wp                 ! density front position in kilometers
      pts(:,:,:,jp_tem) = rn_T1 * ptmask(:,:,:)
      DO jk = 1, jpk
        DO jj = 1, jpj
          DO ji = 1, jpi
            IF ( glamt(ji,jj) <= zdam ) pts(ji,jj,jk,jp_tem) = rn_T0 * ptmask(ji,jj,jk)
          END DO
        END DO
      END DO
      !
      !!an =1. highlight spurious oscillations (centered schemes) but
      ! pts(:,:,:,jp_sal) = 1._wp * ptmask(:,:,:) !
      !!an =T, really close to T
      pts(:,:,:,jp_sal) = pts(:,:,:,jp_tem)
      !
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate

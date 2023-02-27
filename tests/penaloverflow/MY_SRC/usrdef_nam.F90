MODULE usrdef_nam
   !!======================================================================
   !!                       ***  MODULE  usrdef_nam  ***
   !!
   !!                     ===  OVERFLOW configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce  , ONLY: nimpp , njmpp            ! i- & j-indices of the local domain
   USE dom_oce  , ONLY: ln_zco, ln_zps, ln_sco   ! flag of type of coordinate
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called by nemogcm.F90

   !                              !!* namusr_def namelist *!!
   INTEGER, PUBLIC ::   nn_ovf     ! 
   LOGICAL, PUBLIC ::   ln_ovf     ! 

   REAL(wp), PUBLIC ::   rn_dx     ! resolution in meters defining the horizontal domain size
   REAL(wp), PUBLIC ::   rn_dz     ! resolution in meters defining the vertical   domain size
   REAL(wp), PUBLIC ::   rn_T1     ! surrounding temp
   REAL(wp), PUBLIC ::   rn_T0     ! dense temp
   !                              !!* penalisation *!!
   REAL(wp), PUBLIC ::   rn_abp             ! alpha boundary parameter                                       [-]
   INTEGER , PUBLIC ::   nn_abp             ! if T rpo deduced from bathy else indicator function
   INTEGER , PUBLIC ::   nn_cnp             ! number of cell on which is smoothed the porosity (phi)         [-]
   REAL(wp), PUBLIC ::   rn_fsp             ! friction parameter 1/epsilon of the permeability               [1/s]
   INTEGER, PUBLIC  ::   nn_fsp             ! friction parameter 1/epsilon of the permeability               [1/s]
   INTEGER, PUBLIC  ::   nn_wef             ! where friction is applied
   INTEGER , PUBLIC ::   nn_smo          ! number of cell on which is smoothed the porosity (phi)         [-]

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_nam.F90 11536 2019-09-11 13:54:18Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, kperio )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here OVERFLOW configuration
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*)              , INTENT(out) ::   cd_cfg          ! configuration name
      INTEGER                       , INTENT(out) ::   kk_cfg          ! configuration resolution
      INTEGER                       , INTENT(out) ::   kpi, kpj, kpk   ! global domain sizes
      INTEGER                       , INTENT(out) ::   kperio          ! lateral global domain b.c.
      !
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ ln_zco, ln_zps, ln_sco, rn_dx, rn_dz, rn_T1, rn_T0,      &
         &                 rn_abp, nn_abp,  nn_cnp, rn_fsp, nn_fsp, nn_wef, nn_smo, & ! penalisation parameters
         &                 ln_ovf, nn_ovf                                             ! larger bathy

      !!----------------------------------------------------------------------
      !
      REWIND( numnam_cfg )          ! Namelist namusr_def (exist in namelist_cfg only)
      READ  ( numnam_cfg, namusr_def, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namusr_def in configuration namelist' )
      !
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'OVERFLOW'           ! name & resolution (not used)
      kk_cfg = INT( rn_dx )
      !
      ! Global Domain size:  OVERFLOW domain is  200 km x 3 grid-points x 2000 m
      IF ( ln_ovf ) THEN
            IF(nn_ovf<1) CALL ctl_stop('nn_ovf must be >= 1') 
            kpi = INT( 200.e3 / ( REAL(nn_ovf,wp) * rn_dx ) ) + 2               ! [m] gridspacing BIGGER cells
            IF(lwp) WRITE(numout,*) 'before kpi=',kpi
            kpi = kpi * REAL(nn_ovf,wp)
            IF(lwp) WRITE(numout,*) 'after  kpi=',kpi
      ELSE
            kpi = INT( 200.e3 / rn_dx ) + 2
      END IF
      kpj = 3
      kpk = INT(  2000. / rn_dz ) + 1
      !                             ! control print
      WRITE(numout,*) '   '
      WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
      WRITE(numout,*) '~~~~~~~~~~~ '
      WRITE(numout,*) '   Namelist namusr_def : OVERFLOW test case'
      WRITE(numout,*) '      type of vertical coordinate : '
      WRITE(numout,*) '         z-coordinate flag                     ln_zco = ', ln_zco
      WRITE(numout,*) '         z-partial-step coordinate flag        ln_zps = ', ln_zps
      WRITE(numout,*) '         s-coordinate flag                     ln_sco = ', ln_sco
      WRITE(numout,*) '      horizontal resolution                    rn_dx  = ', rn_dx, ' meters'
      WRITE(numout,*) '      vertical   resolution                    rn_dz  = ', rn_dz, ' meters'
      WRITE(numout,*) '      OVERFLOW domain = 200 km x 3 grid-points x 2000 m'
      WRITE(numout,*) '         resulting global domain size :        jpiglo = ', kpi
      WRITE(numout,*) '                                               jpjglo = ', kpj
      WRITE(numout,*) '                                               jpkglo = ', kpk
      !
      !                             ! Set the lateral boundary condition of the global domain
      kperio = 0                    ! OVERFLOW configuration : close basin
      !
      WRITE(numout,*) '   '
      WRITE(numout,*) '   Lateral boundary condition of the global domain'
      WRITE(numout,*) '      OVERFLOW : closed basin                  jperio = ', kperio
      !
#if defined key_ubs
      WRITE(numout,*) ' key_ubs active'
#else
      WRITE(numout,*) ' key_ubs inactive'
#endif
#if defined key_bvp
      WRITE(numout,*) ' key_bvp active'
#else
      WRITE(numout,*) ' key_bvp inactive'
#endif
      !
#if defined key_w_bvp
      WRITE(numout,*) ' key_w_bvp active'
#else
      WRITE(numout,*) ' key_w_bvp inactive'
#endif
      WRITE(numout,*) '   Namelist namusr_def : Brinkman Penalisation Parameters'
      WRITE(numout,*) '                                              rn_abp = ', rn_abp
      WRITE(numout,*) '                                              nn_abp = ', nn_abp
      WRITE(numout,*) '                                              nn_cnp = ', nn_cnp
      WRITE(numout,*) '                                              rn_fsp = ', rn_fsp
      WRITE(numout,*) '                                              nn_fsp = ', nn_fsp
      WRITE(numout,*) '                                              nn_wef = ', nn_wef
      WRITE(numout,*) '                                              nn_smo = ', nn_smo
      WRITE(numout,*) '                                                       '
      WRITE(numout,*) '                                              ln_ovf = ', ln_ovf
      WRITE(numout,*) '                                              nn_ovf = ', nn_ovf
      !
      IF( rn_fsp<=0 .AND. (nn_fsp == 11 .OR. nn_fsp == 21 .OR. nn_fsp == 31 ) ) CALL ctl_stop( 'usr_def_nam: friction cannot be negative or nil with this choice of nn_fsp' )
      !
      IF( .NOT. ln_zco .AND. nn_abp <  1 ) CALL ctl_stop( 'usr_def_nam: choose nn_abp accordingly to ln_zco' )
      IF( .NOT. ln_zps .AND. nn_abp >= 1 ) CALL ctl_stop( 'usr_def_nam: choose nn_abp accordingly to ln_zps' )
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam

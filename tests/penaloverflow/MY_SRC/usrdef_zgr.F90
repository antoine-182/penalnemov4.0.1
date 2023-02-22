MODULE usrdef_zgr
   !!======================================================================
   !!                     ***  MODULE usrdef_zgr  ***
   !!
   !!                       ===  OVERFLOW case  ===
   !!
   !! user defined :  vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-08  (G. Madec, S. Flavoni)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system (required)
   !!       zgr_z1d   : reference 1D z-coordinate
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce ,  ONLY: mi0, mi1, nimpp, njmpp   ! ocean space and time domain
!!an
   USE dom_oce ,  ONLY: glamt, glamu, gphit, gphiu                    ! ocean space and time domain
!!an
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE
   !
!!an
   INTEGER, PARAMETER, PRIVATE  ::   nT  = 1
   INTEGER, PARAMETER, PRIVATE  ::   nU  = 2
   INTEGER, PARAMETER, PRIVATE  ::   nW  = 3
!!an
   !
   PUBLIC   usr_def_zgr   ! called by domzgr.F90

  !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 11077 2019-06-05 14:13:02Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v , pe3f ,               &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw,                      &   !     -      -      -
      &                    k_top  , k_bot    )                             ! top & bottom ocean level
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(in   ) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags ( read in namusr_def )
      LOGICAL                   , INTENT(  out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(  out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors
      INTEGER , DIMENSION(:,:)  , INTENT(  out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   ji, jj, jk,jx,jz         ! dummy indices
      INTEGER  ::   ik                ! local integers
      REAL(wp) ::   zfact, z1_jpkm1   ! local scalar
      REAL(wp) ::   ze3min            ! local scalar
      REAL(wp) ::   z1d, zrpostar,zr1,zr2     ! local scalar
      REAL(wp), DIMENSION(jpi,jpj)      ::   zht, zhu, z2d   ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk)  ::   z3d, z3du,z3d3         ! 3D workspace

      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : OVERFLOW configuration (z(ps)- or s-coordinate closed box ocean without cavities)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ! already set in usrdef_nam.F90 by reading the namusr_def namelist except for ISF
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      !
      !                       !==  UNmasked meter bathymetry  ==!
      !
      ! western continental shelf (500m deep) and eastern deep ocean (2000m deep)
      ! (set through the jpk and jpi (see userdef_nam.F90))
      ! with a hyperbolic tangent transition centered at 40km
      ! NB: here glamt is in kilometers
      !
      ! zht(:,:) = + (  500. + 0.5 * 1500. * ( 1.0 + tanh( (glamt(:,:) - 40.) / 7. ) )  )
      DO ji=1,jpi
        DO jj=1,jpj
          zht(ji,jj) = profilz(glamt(ji,jj))
        END DO
      END DO
      !
      ! at u-point: averaging zht
      DO ji = 1, jpim1
         zhu(ji,:) = 0.5_wp * ( zht(ji,:) + zht(ji+1,:) )
      END DO
      CALL lbc_lnk( 'usrdef_zgr', zhu, 'U', 1. )     ! boundary condition: this mask the surrouding grid-points
!!an
! don't know how but alright
      !                                ! ==>>>  set by hand non-zero value on first/last columns & rows
      DO ji = mi0(1), mi1(1)              ! first row of global domain only
         zhu(ji,2) = zht(ji,2)
      END DO
       DO ji = mi0(jpiglo), mi1(jpiglo)   ! last  row of global domain only
         zhu(ji,2) = zht(ji,2)
      END DO
      zhu(:,1) = zhu(:,2)
      zhu(:,3) = zhu(:,2)
!!an
      !
      !!an no vertical dimension yet
      CALL zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !!an pdet dimensioned
#if defined key_bvp
      ! 1) Definition of the porosity field
      IF ( nn_abp >= 1 ) THEN 
         rpot(:,:,:) = rn_abp
         DO ji = 2, jpim1
            DO jk = 1, jpkm1
               CALL zgr_pse (ji,2,jk,glamu,pdepw_1d,rpot, nT)
            END DO
         END DO
         WHERE ( pdept_1d(:) >= profilz(glamt(1,2)) )   rpot(1  ,2,:) = rn_abp
         WHERE ( pdept_1d(:) >= profilz(glamt(jpi,2)) ) rpot(jpi,2,:) = rn_abp
      ELSE 
         rpot(:,:,:) = 1._wp
         DO ji = 1, jpi
         WHERE ( pdept_1d(:) >= profilz(glamt(ji,2)) ) rpot(ji,2,:) = rn_abp
         ! WHERE ( pdept_1d(:) >= profilz(glamt(ji,2)) ) rpot(ji,2,:) = 1.e-10
         END DO
      ENDIF
      !
      !
      !! Shapiro filter (S=1/2) or linear
      !------------------------ smoothing along x ---------------------!
      !------------------------ ----------------- ---------------------!
      z3d (:,:,:) = rpot(:,:,:)
      DO jx = 1, nn_smo
        IF(lwp) WRITE(numout,*) ' smooth  zer ',jx
        DO jk = 1,jpk
            DO jj = 1,jpj
              DO ji = 2,jpim1
                z3d (ji,jj,jk) = 0.25_wp * rpot(ji-1,jj,jk)+ 0.5_wp * rpot(ji,jj,jk) + 0.25_wp* rpot(ji+1,jj,jk )
                ! z3d (ji,jj,jk) = ( rpot(ji-1,jj,jk) + rpot(ji,jj,jk) + rpot(ji+1,jj,jk ) ) / 3._wp
            END DO
          END DO
        END DO
        CALL lbc_lnk( 'usrdef_zgr', z3d, 'T', 1._wp, kfillmode=jpfillcopy)
        rpot(:,:,:) = z3d(:,:,:)
        !------------------------ smoothing along z ---------------------!
        DO jk = 2,jpkm1
          DO jj = 1,jpj
            DO ji = 1,jpi
                z3d(ji,jj,jk)  = 0.25_wp * rpot(ji,jj,jk-1)+ 0.5_wp * rpot(ji,jj,jk) + 0.25_wp* rpot(ji,jj,jk+1)
                ! z3d(ji,jj,jk)  = ( rpot(ji,jj,jk-1) + rpot(ji,jj,jk) + rpot(ji,jj,jk+1) ) / 3._wp
            END DO
          END DO
        END DO
        rpot(:,:,:) = z3d(:,:,:)
      END DO
      !
      ! CALL lbc_lnk( 'usrdef_zgr', rpot, 'T', 1._wp)
      !
      !!------------------------ broadcast to U and W ---------------------!
      !------------------------ --------------------- ---------------------!
      rpou(:,:,:) = rpot(:,:,:) ; rpow(:,:,:) = rpot(:,:,:)
      IF (nn_cnp == 0) THEN   ! demisomme
        DO jk = 2, jpk
           rpow(:,:,jk) = 0.5_wp * ( rpot(:,:,jk) + rpot(:,:,jk-1) )   ! sens z > 0
        END DO
        DO ji = 1, jpim1
           rpou(ji,:,:) = 0.5_wp * ( rpot(ji,:,:) + rpot(ji+1,:,:) )
        END DO
      ELSE IF (nn_cnp  == 1) THEN ! minimum
        DO jk = 2, jpk
          DO jj= 1, jpj
            DO ji = 1, jpim1
              rpow(ji,jj,jk) = MIN( rpot(ji,jj,jk), rpot(ji,jj,jk-1) )   ! sens z > 0
              rpou(ji,jj,jk) = MIN( rpot(ji,jj,jk), rpot(ji+1,jj,jk) )
            END DO
          END DO
        END DO
      ENDIF
      !!------------------------ inverse ---------------------!
      DO jk = 1, jpk
        DO jj = 1, jpj
           DO ji = 1, jpi
             r1_rpot(ji,jj,jk) = 1._wp / rpot(ji,jj,jk)
             r1_rpou(ji,jj,jk) = 1._wp / rpou(ji,jj,jk)
             r1_rpow(ji,jj,jk) = 1._wp / rpow(ji,jj,jk)
           END DO
        END DO
      END DO
    !
    !
    !!------------------------ impermeability ---------------------!
    !------------------------- -------------- ---------------------!
    !
   bmpu(:,:,:) = 1e-10 ; z3d3(:,:,:) = 1._wp
   ! indicator
   IF (nn_wef == 1) THEN                               ! scale on rux
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 2,jpim1
               z3d3(ji,jj,jk) =  0.5_wp * MAX( ( rpou(ji,jj,jk) + rpou(ji+1,jj,jk) )  , & 
               &                               ( rpou(ji,jj,jk) + rpou(ji-1,jj,jk) )    ) / rpou(ji,jj,jk)  ! U-point
            END DO
         END DO
      END DO
   ELSE IF (nn_wef == 2) THEN                       ! scale on rx ! eq here at max(rpou/rpoti, rpou/rpoti+1 )
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 2,jpi
               !z3d3(ji,jj,jk) =  MAX( rpou(ji,jj,jk), rpou(ji-1,jj,jk) ) / rpot(ji,jj,jk) ! T-point
               z3d3(ji,jj,jk) =  rpou(ji,jj,jk) / MIN(rpot(ji,jj,jk),rpot(ji+1,jj,jk) )    ! U-point (bof car pas contrainte stabolite)
            END DO
         END DO
      END DO
   ENDIF
   !
   IF      ( nn_fsp == 1  .OR. nn_fsp == 2 ) THEN                         ! implicit
      !bmpu(:,:,:) = ABS( 1._wp - z3d3(:,:,:) ) * rn_fsp                    ! T node (r<1, pas instable)
      bmpu(:,:,:) = MAX( z3d3(:,:,:) - 1._wp , 0._wp ) * rn_fsp             ! (T point)
      !bmpu(:,:,:) = ( MAX(z3d3(ji,:,:),z3d3(ji+1,:,:)) - 1._wp ) * rn_fsp ! tombera
   ELSE IF ( nn_fsp == 3 )                    THEN                         ! op. splitting
      !bmpu(:,:,:) = ABS( 1._wp - z3d3(:,:,:) ) * rn_fsp / z3d3(:,:,:)     
      bmpu(:,:,:) = MAX( z3d3(:,:,:) - 1._wp , 0._wp ) * rn_fsp  / z3d3(:,:,:)     
   ! ELSE IF ( nn_fsp == 11 .OR. nn_fsp == 21 ) THEN
   !    bmpu(:,:,:) = 1._wp ; WHERE( z3d3(:,:,:) > 1._wp) bmpu(:,:,:) = rn_fsp                  ! 
   ENDIF
    !
#endif
      !
      !
      !                       !==  top masked level bathymetry  ==!  (all coordinates)
      !
      ! no ocean cavities : top ocean level is ONE, except over land
      ! the ocean basin surrounded by land (1 grid-point) set through lbc_lnk call as jperio=0
      z2d(:,:) = 1._wp                    ! surface ocean is the 1st level
      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )        ! closed basin since jperio = 0 (see userdef_nam.F90)
      k_top(:,:) = NINT( z2d(:,:) )
      !
      !
      !
      IF ( ld_sco ) THEN      !==  s-coordinate  ==!   (terrain-following coordinate)
         !
         k_bot(:,:) = jpkm1 * k_top(:,:)  !* bottom ocean = jpk-1 (here use k_top as a land mask)
         !
         !                                !* terrain-following coordinate with e3.(k)=cst)
         !                                !  OVERFLOW case : identical with j-index (T=V, U=F)
         z1_jpkm1 = 1._wp / REAL( jpkm1 , wp)
         DO jk = 1, jpk
         ! z1_jpkm1 nombre de niveau dans l'ocean
            pdept(:,:,jk) = zht(:,:) * z1_jpkm1 * ( REAL( jk   , wp ) - 0.5_wp )
            pdepw(:,:,jk) = zht(:,:) * z1_jpkm1 * ( REAL( jk-1 , wp )          )
            pe3t (:,:,jk) = zht(:,:) * z1_jpkm1
            pe3u (:,:,jk) = zhu(:,:) * z1_jpkm1
            pe3v (:,:,jk) = zht(:,:) * z1_jpkm1
            pe3f (:,:,jk) = zhu(:,:) * z1_jpkm1
            pe3w (:,:,jk) = zht(:,:) * z1_jpkm1
            pe3uw(:,:,jk) = zhu(:,:) * z1_jpkm1
            pe3vw(:,:,jk) = zht(:,:) * z1_jpkm1
         END DO
      ENDIF
      !
      !
      IF ( ld_zco ) THEN      !==  z-coordinate  ==!   (step-like topography)
         !
#if defined key_bvp
         IF (      nn_abp == -1) THEN 
            k_bot(:,:) = jpkm1 * k_top(:,:)     ! here use k_top as a land mask
         ELSE IF ( nn_abp == 0 ) THEN
            k_bot(:,:) = jpkm1 ! last wet cell
            DO jk = jpkm1, 1, -1
               WHERE( rpot(:,:,jk) <= rn_abp )   k_bot(:,:) = MIN(jk,jpkm1)
            END DO
         ENDIF
#else
          !                                !* bottom ocean compute from the depth of grid-points
          k_bot(:,:) = jpkm1 * k_top(:,:)     ! here use k_top as a land mask
          DO jk = 1, jpkm1
             WHERE( pdept_1d(jk) < zht(:,:) .AND. zht(:,:) <= pdept_1d(jk+1) )   k_bot(:,:) = jk * k_top(:,:)
          END DO
#endif
         !                                !* horizontally uniform coordinate (reference z-co everywhere)
         DO jk = 1, jpk
            pdept(:,:,jk) = pdept_1d(jk)
            pdepw(:,:,jk) = pdepw_1d(jk)
            pe3t (:,:,jk) = pe3t_1d (jk)
            pe3u (:,:,jk) = pe3t_1d (jk)
            pe3v (:,:,jk) = pe3t_1d (jk)
            pe3f (:,:,jk) = pe3t_1d (jk)
            pe3w (:,:,jk) = pe3w_1d (jk)
            pe3uw(:,:,jk) = pe3w_1d (jk)
            pe3vw(:,:,jk) = pe3w_1d (jk)
         END DO
      ENDIF
      !
      !
      IF ( ld_zps ) THEN      !==  zps-coordinate  ==!   (partial bottom-steps)
         !
#if defined key_bvp
      !
      ! 2) definition mask
      IF ( nn_abp >= 1) THEN 
         k_bot(:,:) = jpkm1 * k_top(:,:) ! last wet cell
         DO jk = jpkm1, 1, -1
            WHERE( rpot(:,:,jk) <= rn_abp )   k_bot(:,:) = MIN(jk,jpkm1)
         END DO
      ENDIF
      ! 3) penalisation of the vertical scale factors (done in domain.F90)
      DO jk = 1, jpk                      ! initialization to the reference z-coordinate
         pdept(:,:,jk) = pdept_1d(jk)
         pdepw(:,:,jk) = pdepw_1d(jk)
         pe3t (:,:,jk) = pe3t_1d (jk)
         pe3u (:,:,jk) = pe3t_1d (jk)
         pe3v (:,:,jk) = pe3t_1d (jk)
         pe3f (:,:,jk) = pe3t_1d (jk)
         pe3w (:,:,jk) = pe3w_1d (jk)
         pe3uw(:,:,jk) = pe3w_1d (jk)
         pe3vw(:,:,jk) = pe3w_1d (jk)
      END DO
#else
        ! classic cells
         ze3min = 0.1_wp * rn_dz   ! 10% of the width
         IF(lwp) WRITE(numout,*) '   minimum thickness of the partial cells = 10 % of e3 = ', ze3min
         !
         !                                !* bottom ocean compute from the depth of grid-points
         k_bot(:,:) = jpkm1 ! avant dernier niveau (definit au point T)
         DO jk = jpkm1, 1, -1
            WHERE( zht(:,:) < pdepw_1d(jk) + ze3min )   k_bot(:,:) = jk-1
         END DO
         !
!!an pdep = integral of e3
         !                                !* vertical coordinate system
         DO jk = 1, jpk                      ! initialization to the reference z-coordinate
            pdept(:,:,jk) = pdept_1d(jk)
            pdepw(:,:,jk) = pdepw_1d(jk)
            pe3t (:,:,jk) = pe3t_1d (jk)
            pe3u (:,:,jk) = pe3t_1d (jk)
            pe3v (:,:,jk) = pe3t_1d (jk)
            pe3f (:,:,jk) = pe3t_1d (jk)
            pe3w (:,:,jk) = pe3w_1d (jk)
            pe3uw(:,:,jk) = pe3w_1d (jk)
            pe3vw(:,:,jk) = pe3w_1d (jk)
         END DO
         DO jj = 1, jpj                      ! bottom scale factors and depth at T- and W-points
            DO ji = 1, jpi
              ik = k_bot(ji,jj)
              ! ik    last water cell
              ! ik+1  first land cell (from 0)
              ! derniere cellule mouille - choisit la profondeur la plus fidèle
              pdepw(ji,jj,ik+1) = MIN( zht(ji,jj) , pdepw_1d(ik+1) )
              ! adapt the interface (thickness and column depth)
              ! e3t in the middle of W-point
              pe3t (ji,jj,ik  ) = pdepw(ji,jj,ik+1) - pdepw(ji,jj,ik)

              ! why ? unecessary or couldnt be : pe3t(ji,jj,ik+1) = 2*pe3t(ji,jj,ik+1)-pe3t(ji,jj,ik) ?
              ! careful, with zps, vertical sum of e3t is different from 2000m
              pe3t (ji,jj,ik+1) = pe3t (ji,jj,ik  )
              !
              pdept(ji,jj,ik  ) = pdepw(ji,jj,ik  ) + pe3t (ji,jj,ik  ) * 0.5_wp
              pdept(ji,jj,ik+1) = pdepw(ji,jj,ik+1) + pe3t (ji,jj,ik+1) * 0.5_wp

              ! centering W point on water cell center                                                        |
              ! pe3w (ji,jj,ik+1) = pdept(ji,jj,ik+1) - pdept(ji,jj,ik)              ! = pe3t (ji,jj,ik  )    |
              ! pe3w (ji,jj,ik  ) = pdept(ji,jj,ik  ) - pdept(ji,jj,ik-1)            ! st caution ik > 1      |

              ! change U height to zhu                                                                        ||
              ! z1d = MIN( zhu(ji,jj) , pdepw_1d(ik+1) )                             !                        ||
              ! pe3u (ji,jj,ik)   = z1d - pdepw(ji,jj,ik)                            !                        ||
              ! pe3u (ji,jj,ik+1) = pe3t (ji,jj,ik  )                                !                        ||
            END DO
         END DO
         !                                   ! bottom scale factors and depth at  U-, V-, UW and VW-points
         !                                   ! usually Computed as the minimum of neighbooring scale factors
         pe3u (:,:,:) = pe3t(:,:,:)          ! HERE OVERFLOW configuration :
         pe3v (:,:,:) = pe3t(:,:,:)          !    e3 increases with i-index and identical with j-index
         pe3f (:,:,:) = pe3t(:,:,:)          !    so e3 minimum of (i,i+1) points is (i) point
         pe3uw(:,:,:) = pe3w(:,:,:)          !    in j-direction e3v=e3t and e3f=e3v
         pe3vw(:,:,:) = pe3w(:,:,:)          !    ==>>  no need of lbc_lnk calls

         ! IF(lwp) WRITE(numout,*) 'e3u unchanged'
         ! pe3u (:,:,:) = pe3t(:,:,:)
         !
#endif
      ENDIF
      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z1d  ***
      !!
      !! ** Purpose :   set the depth of model levels and the resulting
      !!      vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ]
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!
      !!            ===    Here constant vertical resolution   ===
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:), INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zt, zw   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z1d : Reference vertical z-coordinates: uniform dz = ', rn_dz
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF
      !
      ! Reference z-coordinate (depth - scale factor at T- and W-points)   ! Madec & Imbard 1996 function
      ! ----------------------
      DO jk = 1, jpk
!!an zw and zt useless
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
!!an
         pdepw_1d(jk) =    rn_dz *   REAL( jk-1 , wp )
         pdept_1d(jk) =    rn_dz * ( REAL( jk-1 , wp ) + 0.5_wp )
         pe3w_1d (jk) =    rn_dz
         pe3t_1d (jk) =    rn_dz
      END DO
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z1d

   SUBROUTINE zgr_pse( ki, kj, kk,                  &
      &                plam,pdepth,prpo,            &
      &                cpoint                       )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_pse  ***
      !!
      !! ** Purpose :   Estimate the porosity field within the cell
      !!
      !! ** Method  :   find the intersection points with the boundaries
      !!                calculate the area land/water
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in   ) ::   ki, kj, kk    ! coordinate of the center of the cell
      REAL, DIMENSION(:,:)      , INTENT(in   ) ::   plam          ! horizontal position    [m]
      REAL, DIMENSION(:)        , INTENT(in   ) ::   pdepth        ! depth array            [m]
      REAL, DIMENSION(:,:,:)    , INTENT(inout) ::   prpo          ! porosity field
      INTEGER                   , INTENT(in   ) ::   cpoint        ! type of point in the C-grid
                                                                   ! T : 1, U : 2, V : 3, F : 4
      INTEGER  ::  ji                              ! dummy loop variables
      REAL     ::  z1d, zxd, zf1                   ! dummy variable
      REAL, DIMENSION(2) ::   zA, zB, zC, zD, zM   ! coordinate array M(x,y)
      REAL               ::   zhA, zhC, zhtt       ! dummy variable
      !!----------------------------------------------------------------------
      !
      !                                      !==  Preparatory work  ==!
      SELECT CASE ( cpoint )                     !* Defining vertices
      CASE ( nT )                                               ! UW coordinates
          zA(1) = plam(ki-1,kj) ; zA(2) = pdepth(kk+1)
          zB(1) = plam(ki  ,kj) ; zB(2) = pdepth(kk+1)
          zC(1) = plam(ki  ,kj) ; zC(2) = pdepth(kk  )
          zD(1) = plam(ki-1,kj) ; zD(2) = pdepth(kk  )
        CASE ( nU )                                             ! W coordinates
          ! not working
          zA(1) = plam(ki-1,kj) ; zA(2) = pdepth(kk+1)
          zB(1) = plam(ki  ,kj) ; zB(2) = pdepth(kk+1)
          zC(1) = plam(ki  ,kj) ; zC(2) = pdepth(kk  )
          zD(1) = plam(ki-1,kj) ; zD(2) = pdepth(kk  )
        CASE ( nW )                                             ! U coordinates
          ! not working
          zA(1) = plam(ki-1,kj) ; zA(2) = pdepth(kk+1)
          zB(1) = plam(ki  ,kj) ; zB(2) = pdepth(kk+1)
          zC(1) = plam(ki  ,kj) ; zC(2) = pdepth(kk  )
          zD(1) = plam(ki-1,kj) ; zD(2) = pdepth(kk  )
      END SELECT
      !
      !      A --------- B
      !      |           |
      !      |     +     |
      !      |  (ki,kk)  |
      !      D --------- C
      !
      ! True height given by the profile
      zhA = profilz(zA(1)) ; zhC = profilz(zC(1)) ;  zhtt = profilz(plam(ki,kj))
      !
      IF      ( zhC < zC(2) )  THEN   ! full land
        z1d = 0._wp
      ELSE IF ( zhA > zA(2) )  THEN   ! full water
        z1d = 1._wp
      ELSE                            ! porous land
         z1d = 0._wp                   ! rectangle integration method
         !
         !  -- + -----------o------------ + --   nn_abp = 1
         !  -- + ---------->| dx/2
         !  -- + -----o-----------o------ + --   nn_abp = 2
         !  -- + ---->| dx/4
         !  -- + --o--------o--------o--- + --   nn_abp = 3
         !  -- + ->| dx/6
         !    zA(1)                     zB(1)
         !
         zxd = zA(1) + 1._wp / ( REAL(nn_abp, wp) * 2._wp )
         ! zf1 is normalised in x and z
         ! warning rn_dx is not applied yet
         DO ji = 1,nn_abp
               zf1 = MIN(1._wp, MAX( 0._wp, (profilz(zxd) - zC(2))/rn_dz) )
               ! IF(lwp) WRITE(numout,*) '               zf1 =',zf1
               ! IF(lwp) WRITE(numout,*) '               xA =',zA(1)
               ! IF(lwp) WRITE(numout,*) '               zxd =',zxd
               ! IF(lwp) WRITE(numout,*) '               z =',profilz(zxd)
               ! IF(lwp) WRITE(numout,*) '               zC =',zC(2)
               z1d = z1d + zf1   / REAL(nn_abp, wp)
               zxd = zxd + 1._wp / REAL(nn_abp, wp)
         END DO
         ! IF(lwp) WRITE(numout,*) '               porous z1d=',z1d
         ! z1d = -1._wp
         ! IF(lwp) WRITE(numout,*) '               porous (ki,kj,kk)=',ki,kj,kk
         ! as profilz is downward, the integral does represent the water fraction
      ENDIF
      !
      prpo(ki,kj,kk) =  MAX(z1d,rn_abp)

   END SUBROUTINE zgr_pse

   FUNCTION profilz(x)  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE profilz  ***
      !!
      !! ** Purpose : topographic slope
      !!
      !! ** Method  :
      !!
      !!----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(in) :: x
      REAL             :: f   ! result
      !!----------------------------------------------------------------------
      !
      f = + (  500. + 0.5 * 1500. * ( 1.0 + tanh( (x - 40.) / 7. ) )  )
      !
   END FUNCTION profilz

   ! FUNCTION profilx(z)  RESULT(f)
   !    !!----------------------------------------------------------------------
   !    !!                 ***  ROUTINE profilz  ***
   !    !!
   !    !! ** Purpose : topographic slope
   !    !!
   !    !! ** Method  :
   !    !!
   !    !!----------------------------------------------------------------------
   !    IMPLICIT NONE
   !    REAL, INTENT(in) :: z
   !    REAL             :: f   ! result
   !    !!----------------------------------------------------------------------
   !    !
   !    f = atanh(2._wp * ( z - 500._wp ) - 1._wp ) * 7._wp + 40._wp
   !    !
   ! END FUNCTION profilx
    !
    FUNCTION profilz_(x,gridw)  RESULT(f)
       !!----------------------------------------------------------------------
       !!                 ***  ROUTINE profilz  ***
       !!
       !! ** Purpose : topographic slope
       !!
       !! ** Method  :
       !!
       !!----------------------------------------------------------------------
       IMPLICIT NONE
       REAL, INTENT(in) :: x
       REAL, DIMENSION(jpk), INTENT(in) :: gridw
       REAL             :: f   ! result
       !!----------------------------------------------------------------------
       !
       f = + (  500. + 0.5 * 1500. * ( 1.0 + tanh( (x - 40.) / 7. ) )  )
       ! f = MIN(f,)
       !
    END FUNCTION profilz_
   !!======================================================================
END MODULE usrdef_zgr

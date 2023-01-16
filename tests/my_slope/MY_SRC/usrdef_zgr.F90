MODULE usrdef_zgr
   !!======================================================================
   !!                     ***  MODULE usrdef_zgr  ***
   !!
   !!                       ===  SLOPE case  ===
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
   USE dom_oce ,  ONLY: glamt                    ! ocean space and time domain
   USE usrdef_nam     ! User defined : namelist variables
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

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
      INTEGER  ::   ji, jj, jk        ! dummy indices
      INTEGER  ::   ik                ! local integers
      REAL(wp) ::   zfact, z1_jpkm1   ! local scalar
      REAL(wp) ::   ze3min, z1d,zx0, z1b, z1c, zl10, zl20, zl11,zl21           ! local scalar
      REAL(wp) ::   za1, za2, za3, za, zb1, zb2, zb, zc1, zc2, zc, zd1, zd, zx,zx1,zx_d, zz_d,zap          ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zht, zhu, z2d   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : SLOPE configuration (z(ps)- or s-coordinate closed box ocean without cavities)'
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
      !     !== SENSIBILITE ==!
      zht(:,:) = 500._wp
      !     !== FULL SLOPE - works well in linear ==!
      ! zht(:,:) = 2000._wp
      ! DO ji=1,jpi
      !   DO jj=1,jpj
      !     ! + 2*rn_dx
      !     z1d = (glamt(ji,jj) + rn_dx*1E-3 ) *rn_p
      !     zht(ji,jj) = MIN(2000._wp,MAX(z1d, 0._wp))  ! la surface doit etre libre
      !     !IF (jj==1) WRITE(*,*) "ji,z1d",ji,z1d,zht(ji,jj)
      !   END DO
      ! END DO
      !     !== SLOPE DROITE ==!
      ! zht(:,:) = 2000._wp
      ! DO ji=1,jpi
      !   DO jj=1,jpj
      !   IF(      glamt(ji,jj) <=  rn_x0 )   THEN
      !     zht(ji,jj) = 500._wp
      !   ELSE IF( glamt(ji,jj) > rn_x0 )   THEN
      !     z1d = 500._wp + (glamt(ji,jj)-rn_x0) *rn_p
      !     zht(ji,jj) = MAX(z1d, 0._wp)
      !   END IF
      !   END DO
      ! END DO
      !     !== SLOPE SPLINE ==!
      ! zl10 = 1._wp ; zl20 = 1._wp ; zx1 = rn_x0 + 1500._wp / rn_p
      ! zl11 = - zl10; zl21 = -zl20
      ! zht(:,:) = 2000._wp
      ! DO ji=1,jpi
      !   DO jj=1,jpj
      !   IF( glamt(ji,jj) <  rn_x0 - zl10 )   THEN
      !     zht(ji,jj) = 500._wp
      !   ELSE IF( glamt(ji,jj) >=  rn_x0 - zl10 .AND. glamt(ji,jj) <  rn_x0 + zl20)   THEN
      !     za1 = (zl10-zl20)/(zl10+zl20)                  ;    zd1 = 2._wp*(zl10*zl10*zl20*zl20)/(zl10+zl20)
      !     za2 = 1._wp/(zl10*zl10 + zl20*zl20 - zl10*zl20)  ;    zd = 500._wp + zd1*za2*za3
      !     za3 = rn_p/4.
      !     za = za3 * za2 * za1
      !     !
      !     zb1 = rn_p/(2._wp*(zl10+zl20))           ;   zc1 = rn_p*zl10/(zl10+zl20)
      !     zb2 = 3._wp*za*(zl10-zl20)               ;   zc2 = 3._wp*za*zl10*zl20
      !     zb = zb1 + zb2                         ;   zc = zc1 - zc2
      !     !
      !     zx         = glamt(ji,jj) - rn_x0
      !     z1d        =        za   * zx
      !     z1d        = (z1d + zb ) * zx
      !     z1d        = (z1d + zc ) * zx
      !     zht(ji,jj) =  z1d + zd
      !     !
      !   ELSE IF( glamt(ji,jj) >= rn_x0 + zl20 .AND. glamt(ji,jj) <  zx1 + zl21 )   THEN
      !     z1d = 500._wp + (glamt(ji,jj)-rn_x0) *rn_p
      !     zht(ji,jj) = MAX(z1d, 0._wp)
      !   ELSE IF( glamt(ji,jj) >=  zx1 + zl21 .AND. glamt(ji,jj) <  zx1 - zl11)   THEN
      !     za1 = (zl11-zl21)/(zl11+zl21)                  ;    zd1 = 2._wp*(zl11*zl11*zl21*zl21)/(zl11+zl21)
      !     za2 = 1._wp/(zl11*zl11 + zl21*zl21 - zl11*zl21)  ;    zd = 2000._wp + zd1*za2*za3
      !     za3 = rn_p/4.
      !     za = za3 * za2 * za1
      !     !
      !     zb1 = rn_p/(2._wp*(zl11+zl21))           ;   zc1 = rn_p*zl11/(zl11+zl21)
      !     zb2 = 3._wp*za*(zl11-zl21)               ;   zc2 = 3._wp*za*zl11*zl21
      !     zb = zb1 + zb2                         ;   zc = zc1 - zc2
      !     !
      !     zx         = glamt(ji,jj) - zx1
      !     z1d        =        za   * zx
      !     z1d        = (z1d + zb ) * zx
      !     z1d        = (z1d + zc ) * zx
      !     zht(ji,jj) =  z1d + zd
      !     !
      !   END IF
      !   END DO
      ! END DO
      !
      ! at u-point: averaging zht
      DO ji = 1, jpim1
         zhu(ji,:) = 0.5_wp * ( zht(ji,:) + zht(ji+1,:) )
      END DO
      CALL lbc_lnk( 'usrdef_zgr', zhu, 'U', 1. )     ! boundary condition: this mask the surrouding grid-points
      !                                ! ==>>>  set by hand non-zero value on first/last columns & rows
      DO ji = mi0(1), mi1(1)              ! first row of global domain only
         zhu(ji,2) = zht(ji,2)
      END DO
       DO ji = mi0(jpiglo), mi1(jpiglo)   ! last row of global domain only
         zhu(ji,2) = zht(ji,2)
      END DO
      zhu(:,1) = zhu(:,2)
      zhu(:,3) = zhu(:,2)
      !
      CALL zgr_z1d( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
#if defined key_bvp
      ! 1) Definition of the porosity field
      rpot(:,:,:) = 1._wp
      DO ji = 1, jpi
        DO jj = 1, jpj
          DO jk = 1, jpk
            !     !== RECTANGLE / FULLSLOPE ==!
            ! gradient normal a zht
            zx0 = pdept_1d(jk)/rn_p - rn_dx*1E-3
            z1d = 0.5_wp * (1._wp+rn_abp) + (1._wp-rn_abp)*(glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp)*rn_dx*1e-3)
            !
            rpot(ji,jj,jk) = MAX(rn_abp*0.5_wp,MIN(1._wp, z1d ) )
            !
            IF (pdept_1d(jk) >= 500._wp .AND. pdept_1d(jk) < 2000._wp) THEN
              !
              !     !== RECTANGLE / slope droite ==!
              ! ! gradient normal a zht
              ! zx0 = (pdept_1d(jk) - 500._wp)/rn_p + rn_x0
              ! z1d = 0.5_wp * (1._wp+rn_abp) + (1._wp-rn_abp)*(glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp)*rn_dx*1e-3)
              ! !
              ! rpot(ji,jj,jk) = MAX(rn_abp*0.5_wp,MIN(1._wp, z1d ) )
              !
              ! !     !== RECTANGLE + SMOOTH / slope droite ==!
              ! zx0  = (pdept_1d(jk) - 500._wp)/rn_p + rn_x0
              ! zx_d = (glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp)*rn_dx*1e-3)  ! vaut -1/2 jusqua'à 1/2
              ! zz_d = (600._wp - pdept_1d(jk))/100._wp                   ! vaut 1 jusqua'à 0
              ! !
              ! zap = rn_abp
              ! ! pas réellement nécessaire au final... mais joli
              ! ! IF (pdept_1d(jk) <= 600._wp .AND. zx_d >= -0.5_wp) zap = (zz_d + rn_abp)/(1._wp + rn_abp)
              ! z1d = 0.5_wp * (1._wp+zap) + (1._wp-zap)*zx_d
              ! rpot(ji,jj,jk) = MAX(rn_abp*0.5_wp,MIN(1._wp, z1d ) )
              !
              !     !== CONE (pic vers le fond)==!
              ! ! phi plus abrupte vers le fond = résolution variant avec la profondeur
              ! z1d = 1._wp
              ! z1c = (pdept_1d(jk)-500._wp)/1500._wp
              ! z1b = 1._wp*(1._wp - z1c) + 0.01_wp*z1c
              !
              ! gradient normal a zht
              ! z1b = 1._wp
              ! zx0 = (pdept_1d(jk) - 500._wp)/20._wp + rn_x0
              ! z1d = 0.5_wp * (1._wp+rn_abp) + (1._wp-rn_abp)*(glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp)*z1b)
              !
              ! rpot(ji,jj,jk) = MAX(rn_abp*0.5_wp,MIN(1._wp, z1d))
              !
              !     !== BRUSH PENTE (pic vers le fond)==!
              ! phi plus abrupte vers le fond = résolution variant avec la profondeur
              ! z1c = (pdept_1d(jk)-500._wp)/1500._wp
              ! z1b = 1._wp*(1._wp - z1c) + 0.01_wp*z1c
              ! !
              ! zx0 = (pdept_1d(jk) - 500._wp)/20._wp + rn_x0
              ! z1d = 0.5_wp * (1._wp+rn_abp*z1b) + (1._wp-rn_abp*z1b)*(glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp))
              ! ! !
              ! rpot(ji,jj,jk) = MAX(rn_abp*0.5*z1b,MIN(1._wp, z1d))
              !
              !     !== CONE (pic vers le haut) ==!
              ! z1d = 0.5_wp * (1._wp+rn_abp) + (1._wp-rn_abp)*(glamt(ji,jj) - zx0)/(REAL(nn_cnp,wp)) ! nice too
              ! rpot(ji,jj,jk) = MAX(rn_abp*0.5_wp,MIN(1._wp, z1d*z1b))
            ELSEIF (pdept_1d(jk) >= 2000._wp) THEN
              rpot(ji,jj,jk) = rn_abp*0.1_wp
            END IF
         END DO
       END DO
      END DO
      CALL lbc_lnk( 'usrdef_zgr', rpot, 'T', 1.)
      !
      DO ji = 1, jpim1
        DO jj = 1, jpj
           DO jk = 1, jpkm1
             rpou(ji,jj,jk) = 0.5_wp * ( rpot(ji,jj,jk) + rpot(ji+1,jj,jk  ) )
             rpow(ji,jj,jk) = 0.5_wp * ( rpot(ji,jj,jk) + rpot(ji  ,jj,jk+1) )
           END DO
        END DO
      END DO
      !! no real need to call for communication, vertical decomposition
      CALL lbc_lnk_multi( 'usrdef_zgr', rpou, 'U', 1., rpow, 'W', 1. )
      !
      DO ji = 1, jpi
        DO jj = 1, jpj
           DO jk = 1, jpk
             r1_rpot(ji,jj,jk) = 1._wp / rpot(ji,jj,jk)
             r1_rpou(ji,jj,jk) = 1._wp / rpou(ji,jj,jk)
             r1_rpow(ji,jj,jk) = 1._wp / rpow(ji,jj,jk)
           END DO
        END DO
      END DO
#endif
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
         !                                !  SLOPE case : identical with j-index (T=V, U=F)
         z1_jpkm1 = 1._wp / REAL( jpkm1 , wp)
         DO jk = 1, jpk
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
        !
        ! 2) definition mask
         k_bot(:,:) = jpkm1 ! avant dernier niveau (defini au point T)
         DO jk = jpkm1, 1, -1
            WHERE( rpot(:,:,jk) <= rn_abp )   k_bot(:,:) = jk-1
         END DO
         !
         ! 2) definition mask (avec pente)
         ! k_bot(:,:) = jpkm1 ! avant dernier niveau (defini au point T)
         ! DO jk = jpkm1, 1, -1
         !   z1c = (pdept_1d(jk)-500._wp)/1500._wp
         !   z1b = 1._wp*(1._wp - z1c) + 0.01_wp*z1c
         !   WHERE( rpot(:,:,jk) <= rn_abp*0.5*z1b )   k_bot(:,:) = jk-1
         ! END DO
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
        k_bot(:,:) = jpkm1 ! avant dernier niveau (defini au point T)
        DO jk = jpkm1, 1, -1
           WHERE( rpot(:,:,jk) <= rn_abp )   k_bot(:,:) = jk-1
        END DO
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
        ! classic definition without penalisation
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
         zw = REAL( jk , wp )
         zt = REAL( jk , wp ) + 0.5_wp
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

   !!======================================================================
END MODULE usrdef_zgr

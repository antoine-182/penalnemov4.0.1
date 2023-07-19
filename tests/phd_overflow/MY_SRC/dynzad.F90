MODULE dynzad
   !!======================================================================
   !!                       ***  MODULE  dynzad  ***
   !! Ocean dynamics : vertical advection trend
   !!======================================================================
   !! History :  OPA  ! 1991-01  (G. Madec) Original code
   !!   NEMO     0.5  ! 2002-07  (G. Madec) Free form, F90
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zad       : vertical advection momentum trend
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition: ocean
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PUBLIC, PARAMETER :: gamma1v = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS

   PUBLIC   dyn_zad       ! routine called by dynadv.F90
   PUBLIC   dyn_zad_up3    ! routine called by dynadv.F90

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynzad.F90 10068 2018-08-28 14:09:04Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_zad ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad  ***
      !!
      !! ** Purpose :   Compute the now vertical momentum advection trend and
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The now vertical advection of momentum is given by:
      !!         w dz(u) = ua + 1/(e1e2u*e3u) mk+1[ mi(e1e2t*wn) dk(un) ]
      !!         w dz(v) = va + 1/(e1e2v*e3v) mk+1[ mj(e1e2t*wn) dk(vn) ]
      !!      Add this trend to the general trend (ua,va):
      !!         (ua,va) = (ua,va) + w dz(u,v)
      !!
      !! ** Action  : - Update (ua,va) with the vert. momentum adv. trends
      !!              - Send the trends to trddyn for diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step inedx
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), DIMENSION(jpi,jpj)     ::   zww
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zwuw, zwvw
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_zad')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zad : 2nd order vertical advection scheme'
      ENDIF

      IF( l_trddyn )   THEN         ! Save ua and va trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = ua(:,:,:)
         ztrdv(:,:,:) = va(:,:,:)
      ENDIF
      DO jk = 2, jpkm1              ! Vertical momentum advection at level w and u- and v- vertical
         DO jj = 2, jpj                   ! vertical fluxes
            DO ji = fs_2, jpi             ! vector opt.
#if defined key_bvp
               zww(ji,jj) = 0.25_wp * e1e2t(ji,jj) * wn(ji,jj,jk) * rpow(ji,jj,jk)
#else
               zww(ji,jj) = 0.25_wp * e1e2t(ji,jj) * wn(ji,jj,jk)
#endif
            END DO
         END DO
         DO jj = 2, jpjm1                 ! vertical momentum advection at w-point
            DO ji = fs_2, fs_jpim1        ! vector opt.
               zwuw(ji,jj,jk) = ( zww(ji+1,jj  ) + zww(ji,jj) ) * ( un(ji,jj,jk-1) - un(ji,jj,jk) )
               zwvw(ji,jj,jk) = ( zww(ji  ,jj+1) + zww(ji,jj) ) * ( vn(ji,jj,jk-1) - vn(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      ! Surface and bottom advective fluxes set to zero
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1           ! vector opt.
            zwuw(ji,jj, 1 ) = 0._wp
            zwvw(ji,jj, 1 ) = 0._wp
            zwuw(ji,jj,jpk) = 0._wp
            zwvw(ji,jj,jpk) = 0._wp
         END DO
      END DO
      !
      DO jk = 1, jpkm1              ! Vertical momentum advection at u- and v-points
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1       ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zwuw(ji,jj,jk) + zwuw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zwvw(ji,jj,jk) + zwvw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO

      IF( l_trddyn ) THEN           ! save the vertical advection trends for diagnostic
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zad, kt )
         DEALLOCATE( ztrdu, ztrdv )
      ENDIF
      !                             ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' zad  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_zad')
      !
   END SUBROUTINE dyn_zad



   SUBROUTINE dyn_zad_up3( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad_up3  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfu_t, zfu_f, zfu_uw, zfu
      REAL(wp), DIMENSION(jpi,jpj,jpk)   ::   zfv_t, zfv_f, zfv_vw, zfv, zfw
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlu_uu, zlu_uv
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlv_vv, zlv_vu
      REAL(wp) ::   zfuk, zfvk, z1d
      REAL(wp), DIMENSION(jpi,jpj,jpk,2) ::   zlu_uw, zlv_vw,zlw_uw, zlw_vw
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_adv_up3 : UP3 flux form momentum advection on all directions'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      zfu_t(:,:,:) = 0._wp
      zfv_t(:,:,:) = 0._wp
      zfu_f(:,:,:) = 0._wp
      zfv_f(:,:,:) = 0._wp
      !
      zlu_uu(:,:,:,:) = 0._wp
      zlv_vv(:,:,:,:) = 0._wp
      zlu_uv(:,:,:,:) = 0._wp
      zlv_vu(:,:,:,:) = 0._wp
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = ua(:,:,:)
         zfv_vw(:,:,:) = va(:,:,:)
      ENDIF
      !
      !
#if defined key_bvp
      DO jk = 1, jpk            ! vertical volume fluxes penalised !
         zfw(:,:,jk) = 0.25_wp * e1e2t(:,:) * wn(:,:,jk) * rpow(:,:,jk)
      END DO
#else
      DO jk = 1, jpk            ! vertical volume fluxes
         zfw(:,:,jk) = 0.25_wp * e1e2t(:,:) * wn(:,:,jk)
      END DO
#endif
      !
     zlu_uw(:,:, 1   ,1) = 0._wp   ;   zlw_uw(:,:, 1   ,2) = 0._wp
     zlv_vw(:,:, 1   ,1) = 0._wp   ;   zlw_vw(:,:, 1   ,2) = 0._wp
     !                           ! =========== !
     DO jk = 2, jpkm1            !  Laplacian  !
       !                         ! =========== !
       !
         zlu_uw(:,:,jk,1) = ( ub(:,:,jk-1) - 2._wp*ub (:,:,jk) + ub (:,:,jk+1) ) * umask(:,:,jk)   ! U
         zlv_vw(:,:,jk,1) = ( vb(:,:,jk-1) - 2._wp*vb (:,:,jk) + vb (:,:,jk+1) ) * vmask(:,:,jk)   ! V
         !
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = fs_2, fs_jpim1   ! vector opt.
               zlw_uw(ji,jj,jk,2) = ( zfw (ji+1,jj,jk) - 2.*zfw (ji,jj,jk) + zfw (ji-1,jj,jk) ) * wmask(ji,jj,jk)   ! W
               zlw_vw(ji,jj,jk,2) = ( zfw (ji,jj+1,jk) - 2.*zfw (ji,jj,jk) + zfw (ji,jj-1,jk) ) * wmask(ji,jj,jk)   ! W
          END DO
        END DO
     END DO
      ! no need to call  lbc_lnk
      !
      !z1d = 1._wp/32._wp
      DO jk = 2, jpkm1
         DO jj = 1, jpjm1
            DO ji = 1, fs_jpim1   ! vector opt.
               !
               zfuk = ( zfw(ji+1,jj,jk) + zfw(ji,jj,jk) )  ! vertical transport en ji+1/2,jj    ,jk+1/2 - UW
               zfvk = ( zfw(ji,jj+1,jk) + zfw(ji,jj,jk) )  !                       ji    ,jj+1/2,jk+1/2 - VW
               IF( zfuk > 0._wp ) THEN   ;    zl_u = zlu_uw( ji,jj,jk  ,1)
               ELSE                      ;    zl_u = zlu_uw( ji,jj,jk-1,1)
               ENDIF
               IF( zfvk > 0._wp ) THEN   ;    zl_v = zlv_vw( ji,jj,jk  ,1)
               ELSE                      ;    zl_v = zlv_vw( ji,jj,jk-1,1)
               ENDIF
               !                          C2 - upstream
               zfu_uw(ji  ,jj  ,jk) = zfuk * ( un(ji,jj,jk) + un(ji,jj,jk-1) - gamma1v * zl_u ) ! UW nodes
               zfv_vw(ji  ,jj  ,jk) = zfvk * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) - gamma1v * zl_v ) ! VW nodes
            END DO
         END DO
      END DO
      !
      DO jk = 1, jpkm1                          ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.
      !!an RHS ecrit en vector form, les e3u seront re-multiplie derriere pour avoir un transport
               ua(ji,jj,jk) =  ua(ji,jj,jk) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj) / e3u_n(ji,jj,jk)
               va(ji,jj,jk) =  va(ji,jj,jk) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj) / e3v_n(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                       ! save the vertical advection trend for diagnostic
         zfu_t(:,:,:) = ua(:,:,:) - zfu_t(:,:,:)
         zfv_t(:,:,:) = va(:,:,:) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt )
      ENDIF
      !                                         ! Control print
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' ubs2 adv - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_zad_up3


   !!======================================================================
END MODULE dynzad

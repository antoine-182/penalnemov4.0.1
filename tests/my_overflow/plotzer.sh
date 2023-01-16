#!/bin/bash
# FROM the config

# Implicit
pgear rhoOverflow -n EXP_ref/zco_imp/OVF_zco_imp_grid_T.nc EXP_ref/zps_10_imp/OVF_zps_imp_grid_T.nc -j 0.15 -s
pgear rhoOverflow -n EXP_bvp/bvp_imp/OVF_bvp_imp_grid_T.nc EXP_ref/zps_10_imp/OVF_zps_imp_grid_T.nc -j 0.15 -s
#
#Nominal
pgear rhoOverflow -n EXP_ref/zco/OVF_ref_zco_grid_T.nc EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
pgear rhoOverflow -n EXP_bvp/bvp/OVF_bvp_grid_T.nc     EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s

# list_rot="imp"
list_rot="x1z0 x2z0 x3z0 x1z1 x2z2 x3z3"
# list_rot="imp"
for name in $list_rot;
do
  directory=bvp_${name}
  # no need the corresponding mesh_mask
  pgear rhoOverflow -n EXP_bvp/${directory}/OVF_${directory}_grid_T.nc EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
done

#Param - ABBL
pgear rhoOverflow -n EXP_ref/zps_trabbl/OVF_zps_trabbl_grid_T.nc         EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
#Param - ABBLx2
pgear rhoOverflow -n EXP_ref/zps_trabblx2/OVF_trabblx2_grid_T.nc     EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
#Param - NPCA
pgear rhoOverflow -n EXP_ref/zps_npc/OVF_zps_npc_grid_T.nc               EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
#Param - ABBL+NPCA
pgear rhoOverflow -n EXP_ref/zps_trabbl_npc/OVF_zps_trabbl_npc_grid_T.nc EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s

#test
name="x3z3"
directory=bvp_${name}
pgear rhoOverflow -n EXP_ref/zps_trabbl_npc/OVF_zps_trabbl_npc_grid_T.nc EXP_bvp/${directory}/OVF_${directory}_grid_T.nc -j 0.15 -s

pgear rhoOverflow -n EXP_bvp/bvp_x3z0/OVF_bvp_x3z0_grid_T.nc EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s

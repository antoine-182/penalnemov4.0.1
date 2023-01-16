#!/bin/bash
# FROM the config

# Implicit
# pgear rhoOverflow -n EXP_ref/zco_imp/OVF_zco_imp_grid_T.nc EXP_ref/zps_10_imp/OVF_zps_imp_grid_T.nc -j 0.15 -s
# pgear rhoOverflow -n EXP_bvp/bvp_imp/OVF_bvp_imp_grid_T.nc EXP_ref/zps_10_imp/OVF_zps_imp_grid_T.nc -j 0.15 -s
# #
# #Nominal
# pgear rhoOverflow -n EXP_ref/zco/OVF_ref_zco_grid_T.nc EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
# pgear rhoOverflow -n EXP_bvp/bvp/OVF_bvp_grid_T.nc     EXP_ref/zps_10/OVF_ref_zps_grid_T.nc -j 0.15 -s
#
# list_rot="5m 10m 20m 50m 100m 200m 250m"
# for name in $list_rot:
# do
#   directory=zco_${name}
#   echo $directory
#   cd $directory   # no need the corresponding mesh_mask
#   pgear toceOverflow SLO_${directory}_*grid_T.nc -f
#   pgear sfOverflow SLO_${directory}_*grid_U.nc -f
#   cd ../
# done

# list_rot="5m 10m 20m 50m 100m"
# serie=""
# for name in $list_rot:
# do
#   directory=zco_${name}
#   echo $directory
#   serie+=" EXP_ref/${directory}/SLO_${directory}*grid_T.nc"   # no need the corresponding mesh_mask
# done
# pgear rhoOverflow -n ${serie} -s -e

cd EXP_ref
# list_rot="5m 10m 20m 50m 100m"
list_rot="100m"
for name in $list_rot:
do
  directory="fullslope/${name}"
  echo $directory
  cd $directory   # no need the corresponding mesh_mask
  pgear toceOverflow SLO*grid_T.nc -f
  pgear sfOverflow SLO*grid_U.nc -f
  cd ../../
done

#!/bin/bash

list_rot="x1z0 x2z0 x3z0 x1z1 x2z2 x3z3"
for name in $list_rot;
do
  directory=bvp_${name}
  #
  echo
  echo ${directory}
  cd ${directory}
  pgear plothistrhox OVF_${directory}_grid_T.nc -m ../meshmask/mesh_mask_${name}.nc -s
  cd ../
done

#!/bin/bash

# list_rot="imp"
list_rot="imp x1z0 x2z0 x3z0 x1z1 x2z2 x3z3"
# list_rot="imp"
for name in $list_rot;
do
  directory=bvp_${name}
  #
  echo ${directory}
  cd ${directory}
  # no need the corresponding mesh_mask
  pgear plotoverflow OVF_${directory}_grid_T.nc -m ../meshmask/mesh_mask.nc -s -f
  # pgear rpoover OVF_${directory}_grid_T.nc -m ../meshmask/mesh_mask.nc -s
  cd ../
done

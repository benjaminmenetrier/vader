#!/bin/bash

# Loop over all files in the specified directory
for file in ./*.nc; do
  echo $file
  ncks $file ${file}_out.nc
  ncrename -v potential_temperature,air_potential_temperature -d nz_potential_temperature,nz_air_potential_temperature -h ${file}_out.nc
  ncrename -v exner_levels_minus_one,dimensionless_exner_function_levels_minus_one -d nz_exner_levels_minus_one,nz_dimensionless_exner_function_levels_minus_one -h ${file}_out.nc
  ncrename -v exner,dimensionless_exner_function -d nz_exner,nz_dimensionless_exner_function -h ${file}_out.nc
  mv ${file}_out.nc $file
done

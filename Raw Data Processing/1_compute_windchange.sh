#!/bin/bash

# Take the Wind intensity computed before - NOT Deseasonalized of course 

# Take the period 1999-2010 : Before Bloom

cdo selyear,1999/2010 Modulo_of_wind_era5_1993-2022_f64.nc Modulo_of_wind_era5_1993-2022_f64_before.nc

# compute multi-year seasonal means before bloom 
cdo yseasmean Modulo_of_wind_era5_1993-2022_f64_before.nc yseasm_before.nc

# Take the period 2011-2022 : After Bloom 

cdo selyear,2011/2022 Modulo_of_wind_era5_1993-2022_f64.nc Modulo_of_wind_era5_1993-2022_f64_after.nc 

# compute the multi-year seasonal means after bloom
cdo yseasmean Modulo_of_wind_era5_1993-2022_f64_after.nc yseasm_after.nc 


# Finally compute the Changes in Multi-year seasonal means (after - before)

cdo sub yseasm_after.nc yseasm_before.nc changes_yseasm_wind_after-before.nc

# Comment/Uncomment based on what you want to save or not 
rm Modulo_of_wind_era5_1993-2022_f64_after.nc
rm Modulo_of_wind_era5_1993-2022_f64_before.nc 
rm yseasm_after.nc
#rm yseasm_before.nc

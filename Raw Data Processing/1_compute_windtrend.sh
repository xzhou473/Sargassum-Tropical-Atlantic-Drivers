#!/bin/bash 

# copy data as downloaded from there to here: 
cp ../data_wind/wind_data_era5/adaptor.mars.internal-1717099791.002421-20723-10-94a46487-c083-4c6e-b6f7-4bc1a5fe70ab.nc . 


# change precision and use a meaningful name

cdo -b 64 -copy adaptor.mars.internal-1717099791.002421-20723-10-94a46487-c083-4c6e-b6f7-4bc1a5fe70ab.nc wind_era5_1993-2022_f64_unmasked.nc

# rm useless to save space 
rm adaptor.mars.internal-1717099791.002421-20723-10-94a46487-c083-4c6e-b6f7-4bc1a5fe70ab.nc 

# Mask out the Pacific Ocean grid points and remove temporary useles files 
cdo maskregion,mymasknew.dat wind_era5_1993-2022_f64_unmasked.nc wind_era5_1993-2022_f64_notbox.nc 
rm wind_era5_1993-2022_f64_unmasked.nc

# Select lon lat area of interest 

cdo sellonlatbox,-90,0,0,20 wind_era5_1993-2022_f64_notbox.nc wind_era5_1993-2022_f64.nc
rm wind_era5_1993-2022_f64_notbox.nc 

# First, Separate between the two compnents u and v 

cdo selvar,u wind_era5_1993-2022_f64.nc Uwind_era5_1993-2022_f64.nc
cdo selvar,v wind_era5_1993-2022_f64.nc Vwind_era5_1993-2022_f64.nc
rm wind_era5_1993-2022_f64.nc

# Second, Remove seasonal cycle from each component 
cdo -L -ymonsub Uwind_era5_1993-2022_f64.nc -ymonmean Uwind_era5_1993-2022_f64.nc Uwind_era5_1993-2022_f64_deseas.nc

cdo -L -ymonsub Vwind_era5_1993-2022_f64.nc -ymonmean Vwind_era5_1993-2022_f64.nc Vwind_era5_1993-2022_f64_deseas.nc

# Third, compute trend maps of the deseasonalized components, separately first 

cdo trend Uwind_era5_1993-2022_f64_deseas.nc a TrendMap_of_Uwind_era5_1993-2022_f64_deseas.nc
cdo trend Vwind_era5_1993-2022_f64_deseas.nc b TrendMap_of_Vwind_era5_1993-2022_f64_deseas.nc

rm a 
rm b 

# Additionally, compute the trend of the deseasonalized Wind Intensity, as follows: 

# 1. compute the module of the Wind using U and V (non deseasonalized)
# mod = sqrt(u^2 + v^2) 

# faccio il quadrato:
cdo sqr Uwind_era5_1993-2022_f64.nc u_alla2.nc 
cdo sqr Vwind_era5_1993-2022_f64.nc v_alla2.nc

# sommo i quadrati
cdo add u_alla2.nc v_alla2.nc u2_piu_v2.nc

# faccio la radice e trovo il modulo:
cdo sqrt u2_piu_v2.nc Modulo_of_wind_era5_1993-2022_f64.nc

# 2. Remove the seasonal cycle from the Wind module: 

cdo -L -ymonsub Modulo_of_wind_era5_1993-2022_f64.nc -ymonmean Modulo_of_wind_era5_1993-2022_f64.nc Modulo_of_wind_era5_1993-2022_f64_deseas.nc

#rm u2_piu_v2.nc 
#rm u_alla2.nc 
#rm v_alla2.nc 

# 3. Trend of the deseasonalized wind modulus: 
cdo trend Modulo_of_wind_era5_1993-2022_f64_deseas.nc a TrendMap_WindMod_era5_1993-2022_f64_deseas.nc


rm a 
#rm  Modulo_of_wind_era5_1993-2022_f64_deseas.nc
#rm  Uwind_era5_1993-2022_f64_deseas.nc
#rm  Vwind_era5_1993-2022_f64_deseas.nc
rm wind_era5_1993-2022_f64.nc

# Additional useful calculations: 


# =============================================================
#   TRENDS OF SEASONAL MEANS, in case one wants to check them:
# =============================================================

# First I compute seasonal means of deseasonalized data, so that for each year I have 4 timepoints
cdo seasmean Modulo_of_wind_era5_1993-2022_f64_deseas.nc seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas.nc

# Now I separate the seasons 
cdo splitseas seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas.nc seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_

# Now, because points are equally spaced, I can compute trends in each season separately 

cdo trend seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_DJF.nc a TrendMap_seasm_DJF_wind.nc 
cdo trend seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_MAM.nc a TrendMap_seasm_MAM_wind.nc
cdo trend seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_JJA.nc a TrendMap_seasm_JJA_wind.nc
cdo trend seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_SON.nc a TrendMap_seasm_SON_wind.nc

# Comment/Uncomment based on what you want to remove or save
rm seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas.nc
rm seasm_of_Modulo_of_wind_era5_1993-2022_f64_deseas_*
rm a 

rm Uwind_era5_1993-2022_f64_deseas.nc
rm Vwind_era5_1993-2022_f64_deseas.nc
#rm Uwind_era5_1993-2022_f64.nc
#rm Vwind_era5_1993-2022_f64.nc
rm TrendMap_of_Uwind_era5_1993-2022_f64_deseas.nc
rm TrendMap_of_Vwind_era5_1993-2022_f64_deseas.nc





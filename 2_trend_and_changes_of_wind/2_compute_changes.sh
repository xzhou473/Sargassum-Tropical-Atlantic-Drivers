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


# ---- Added part below ---- 
# Compute % changes in multi-year seasonal standard deviation of detrended wind
# 

# Take the period 1999-2010 : Before Bloom
cdo selyear,1999/2010 Modulo_of_wind_era5_1993-2022_f64.nc Modulo_of_wind_era5_1993-2022_f64_before.nc

# remove trend of this subperiod
cdo detrend Modulo_of_wind_era5_1993-2022_f64_before.nc Modulo_of_wind_era5_1993-2022_f64_before_detr.nc

# compute multi-year seasonal standard dev before 
cdo yseasstd Modulo_of_wind_era5_1993-2022_f64_before_detr.nc yseasstd_of_detr_wind_before.nc


# Take the period 2011-2022 : After Bloom
cdo selyear,2011/2022 Modulo_of_wind_era5_1993-2022_f64.nc Modulo_of_wind_era5_1993-2022_f64_after.nc

# remove trend of this subperiod
cdo detrend Modulo_of_wind_era5_1993-2022_f64_after.nc Modulo_of_wind_era5_1993-2022_f64_after_detr.nc

# compute mulit-year seasonal standard dev after 
cdo yseasstd Modulo_of_wind_era5_1993-2022_f64_after_detr.nc yseasstd_of_detr_wind_after.nc 

# Percentual Changes of Variability season by season: 

cdo sub yseasstd_of_detr_wind_after.nc yseasstd_of_detr_wind_before.nc var_after-bef.nc
cdo div var_after-bef.nc yseasstd_of_detr_wind_before.nc fraction.nc 
cdo mulc,100 fraction.nc perc_changes_variability_wind.nc

# Percentual change of variability all together, i.e. not separate by season. 
# (i.e. considering the period before bloom detrended and after bloom detrended)

cdo timstd Modulo_of_wind_era5_1993-2022_f64_before_detr.nc std_whole_wind_before.nc
cdo timstd Modulo_of_wind_era5_1993-2022_f64_after_detr.nc std_whole_wind_after.nc 

cdo sub std_whole_wind_after.nc std_whole_wind_before.nc diff_std_whole_after-bef.nc
cdo div diff_std_whole_after-bef.nc std_whole_wind_before.nc fraction_whole.nc 
cdo mulc,100 fraction_whole.nc perc_changes_variability_whole_wind.nc


# remove useless stuff
rm fraction
rm var_after-bef.nc
rm yseasstd_of_detr_wind_before.nc
rm yseasstd_of_detr_wind_after.nc
rm Modulo_of_wind_era5_1993-2022_f64_after_detr.nc
rm Modulo_of_wind_era5_1993-2022_f64_before_detr.nc
rm Modulo_of_wind_era5_1993-2022_f64_after.nc
rm Modulo_of_wind_era5_1993-2022_f64_before.nc


# ----- Another added part ------
# Changes in extremes 
# 

# Take the period Before Bloom
cdo selyear,1999/2010 Modulo_of_wind_era5_1993-2022_f64.nc wind_before.nc

# Take the period 2011-2022 : After Bloom
cdo selyear,2011/2022 Modulo_of_wind_era5_1993-2022_f64.nc wind_after.nc

# seasonal maximums in the Before period 

cdo yseasmax wind_before.nc yseasmax_before.nc 

# Select in the after period, in each season separately, only the points that exceed the threshold just computed. Then count them, i.e. find how many times it is exceeded

cdo splitseas yseasmax_before.nc yseasmax_bef_

cdo splitseas wind_after.nc wind_after_

# mask with 1 if primo input > secondo input
cdo gt wind_after_DJF.nc yseasmax_bef_DJF.nc mask_DJF.nc
cdo gt wind_after_MAM.nc yseasmax_bef_MAM.nc mask_MAM.nc
cdo gt wind_after_JJA.nc yseasmax_bef_JJA.nc mask_JJA.nc
cdo gt wind_after_SON.nc yseasmax_bef_SON.nc mask_SON.nc


# now count 
cdo timsum mask_DJF.nc summed_DJF_tmp.nc
cdo timsum mask_MAM.nc summed_MAM_tmp.nc
cdo timsum mask_JJA.nc summed_JJA_tmp.nc
cdo timsum mask_SON.nc summed_SON_tmp.nc

# fix dates because they mess up when doing the sum, bu I know the seasons are correct so: 
cdo setdate,2022-01-01 summed_DJF_tmp.nc summed_DJF.nc
cdo setdate,2022-04-01 summed_MAM_tmp.nc summed_MAM.nc
cdo setdate,2022-06-01 summed_JJA_tmp.nc summed_JJA.nc
cdo setdate,2022-10-01 summed_SON_tmp.nc summed_SON.nc

# what percentage of the total number of months in the After period in each season this summed numbers represent 
# percent_occurr = 100(summed/36) because each period is 12 years, each season has 3 months so 12*3= 36 months total in each season in the after period.

cdo divc,36 summed_DJF.nc fraction_DJF.nc 
cdo divc,36 summed_MAM.nc fraction_MAM.nc
cdo divc,36 summed_JJA.nc fraction_JJA.nc
cdo divc,36 summed_SON.nc fraction_SON.nc

cdo mulc,100 fraction_DJF.nc changes_extremes_wind_DJF.nc 
cdo mulc,100 fraction_MAM.nc changes_extremes_wind_MAM.nc
cdo mulc,100 fraction_JJA.nc changes_extremes_wind_JJA.nc
cdo mulc,100 fraction_SON.nc changes_extremes_wind_SON.nc

#rm fraction*
#rm mask_*
#rm yseasmax_bef*
#rm summed_*
rm wind_after_*
rm wind_before_*



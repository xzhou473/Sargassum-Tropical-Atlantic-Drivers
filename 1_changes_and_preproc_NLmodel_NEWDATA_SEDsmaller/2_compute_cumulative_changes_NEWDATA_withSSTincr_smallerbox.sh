#!/bin/bash


# I remap the dust on the 025°grid, just because calculations need to be on the same grid. 
cdo remapbil,grid025new.dat yseasm_changes_2011-2019_vs_2003-2011_dust.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc

# A. MLT ad Dust: 
# Select a reduced lonlat box for moth mlt and dust, with less lat extent: 

cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc
cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2022_vs_1999-2010_mlt.nc yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc

#rm yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc
#rm yseasm_changes_2011-2022_vs_1999-2010_mlt.nc

# I want only the grid point where mlt and dust changes in yseasm are positive 

# So I create a mask: 1 where greater than 0; put 0 elsewhere. 

cdo gtc,0 yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc maskdust.nc 
cdo gtc,0 yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc maskmlt.nc  

# B. SST: I only want the grid points where the changes are POSITIVE : This is if we hyp. that Sargassum thrives with increasing average Temp


# Again, first select a smaller lonlat box, with less latitudinal extent also for SST: 

# Celsius version
cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2022_vs_1999-2010_sst.nc yseasm_changes_2011-2022_vs_1999-2010_sst_smaller.nc
rm yseasm_changes_2011-2022_vs_1999-2010_sst.nc

# Then create the mask: put a 1 when sst changes are positive
cdo gtc,0 yseasm_changes_2011-2022_vs_1999-2010_sst_smaller.nc masksst.nc


# Kelvin version 

cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2022_vs_1999-2010_sst_Kelvin.nc yseasm_changes_2011-2022_vs_1999-2010_sst_smaller_Kelvin.nc
#rm yseasm_changes_2011-2022_vs_1999-2010_sst_Kelvin.nc

# Then create the mask: put a 1 when sst changes are positive
cdo gtc,0 yseasm_changes_2011-2022_vs_1999-2010_sst_smaller_Kelvin.nc masksst_Kelvin.nc


# C. Only keep the points where the mask is 1 for all the variables. 

cdo mul maskdust.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc dust_changes_positive_2011-2019_vs_2003-2011.nc
cdo mul maskmlt.nc yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc mlt_changes_positive_2011-2022_vs_1999-2010.nc
cdo mul masksst.nc yseasm_changes_2011-2022_vs_1999-2010_sst_smaller.nc sst_changes_positive_2011-2022_vs_1999-2010.nc
cdo mul masksst_Kelvin.nc yseasm_changes_2011-2022_vs_1999-2010_sst_smaller_Kelvin.nc sst_changes_positive_2011-2022_vs_1999-2010_Kelvin.nc
# Compute a measure for cumulative change, in each season separately. 
# Similar to SED but do not sum over seasons

# winter_change = sqrt(sum_over_indicators(ind_i_winter/p95(|ind_i_winter|))^2
# spring_change = sqrt(sum_over_indicators(ind_i_spring/p95(|ind_i_winter|))^2
# summer_change = sqrt(sum_over_indicators(ind_i_summer/p95(|ind_i_winter|))^2
# fall_change = sqrt(sum_over_indicators(ind_i_fall/p95(|ind_i_winter|))^2

cdo griddes yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc > grid025new_smaller.dat

# Compute p95 of the positive-only changes in each variable in each season. 
# Since everything is positive, no need to compute the absolute value before calculating the 95th percentile. 

cdo fldpctl,95 dust_changes_positive_2011-2019_vs_2003-2011.nc p95_dust.nc
cdo enlarge,grid025new_smaller.dat p95_dust.nc gridded_p95_dust.nc

cdo fldpctl,95 mlt_changes_positive_2011-2022_vs_1999-2010.nc p95_mlt.nc
cdo enlarge,grid025new_smaller.dat p95_mlt.nc gridded_p95_mlt.nc

cdo fldpctl,95 sst_changes_positive_2011-2022_vs_1999-2010.nc p95_sst.nc
cdo enlarge,grid025new_smaller.dat p95_sst.nc gridded_p95_sst.nc

cdo fldpctl,95 sst_changes_positive_2011-2022_vs_1999-2010_Kelvin.nc p95_sst_Kelvin.nc
cdo enlarge,grid025new_smaller.dat p95_sst_Kelvin.nc gridded_p95_sst_Kelvin.nc

cdo div dust_changes_positive_2011-2019_vs_2003-2011.nc gridded_p95_dust.nc dust_changes_norm.nc
cdo div mlt_changes_positive_2011-2022_vs_1999-2010.nc gridded_p95_mlt.nc mlt_changes_norm.nc
cdo div sst_changes_positive_2011-2022_vs_1999-2010.nc gridded_p95_sst.nc sst_changes_norm.nc
cdo div sst_changes_positive_2011-2022_vs_1999-2010_Kelvin.nc gridded_p95_sst.nc sst_changes_norm_Kelvin.nc

cdo sqr dust_changes_norm.nc dust_changes_norm_alla2.nc
cdo sqr mlt_changes_norm.nc mlt_changes_norm_alla2.nc
cdo sqr sst_changes_norm.nc sst_changes_norm_alla2.nc
cdo sqr sst_changes_norm_Kelvin.nc sst_changes_norm_alla2_Kelvin.nc

# Celsius version
cdo add dust_changes_norm_alla2.nc mlt_changes_norm_alla2.nc sum_partial_alla2.nc 
cdo add sum_partial_alla2.nc sst_changes_norm_alla2.nc sum_alla_2.nc 

cdo sqrt sum_alla_2.nc cumulative_positive_change_in_each_season_with_SSTincrease.nc

# Kelvin version
cdo add sum_partial_alla2.nc sst_changes_norm_alla2_Kelvin.nc sum_alla_2_Kelvin.nc

cdo sqrt sum_alla_2_Kelvin.nc cumulative_positive_change_in_each_season_with_SSTincrease_Kelvin.nc


# Uncomment to clear some data - leave it commented to save them instead
rm *_norm_*
rm sum_partial*
rm sum_alla*
rm p95_*



#!/bin/bash

# I remap the dust on the 0.25°grid, just because calculations need to be on the same grid (and mlt is already on that grid)
cdo remapbil,grid025new.dat yseasm_changes_2011-2019_vs_2003-2011_dust.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc

# Select a smaller lonlat box (with smaller latitudinal extent) for both dust and mlt:
cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc
cdo sellonlatbox,-90,0,0,20 yseasm_changes_2011-2022_vs_1999-2010_mlt.nc yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc

#rm  yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid.nc 
#rm  yseasm_changes_2011-2022_vs_1999-2010_mlt.nc

# I want only the grid point where mlt and dust changes in yseasm are positive 

# create a mask: 1 where greater than 0; put 0 elsewhere. 

cdo gtc,0 yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc maskdust.nc 
cdo gtc,0 yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc maskmlt.nc  

# Only keep the points where the mask is 1. 

cdo mul maskdust.nc yseasm_changes_2011-2019_vs_2003-2011_dust_on025grid_smaller.nc dust_changes_positive_2011-2019_vs_2003-2011.nc
cdo mul maskmlt.nc yseasm_changes_2011-2022_vs_1999-2010_mlt_smaller.nc mlt_changes_positive_2011-2022_vs_1999-2010.nc

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

cdo div dust_changes_positive_2011-2019_vs_2003-2011.nc gridded_p95_dust.nc dust_changes_norm.nc
cdo div mlt_changes_positive_2011-2022_vs_1999-2010.nc gridded_p95_mlt.nc mlt_changes_norm.nc

cdo sqr dust_changes_norm.nc dust_changes_norm_alla2.nc
cdo sqr mlt_changes_norm.nc mlt_changes_norm_alla2.nc

cdo add dust_changes_norm_alla2.nc mlt_changes_norm_alla2.nc sum_alla_2.nc 

cdo sqrt sum_alla_2.nc cumulative_positive_change_in_each_season_without_SST.nc



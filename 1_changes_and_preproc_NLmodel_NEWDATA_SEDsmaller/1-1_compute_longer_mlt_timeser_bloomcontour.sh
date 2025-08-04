#!/bin/bash

# ------------------------------------------------
# ONLY preprocessing of MLT data:
# ------------------------------------------------

# 1. mlt extended from Jan 1993 to Dec 2024, to compute the timeseries. 

# Copy data at original resolution from there to here

# until june 2021:
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc .

# interim, from july 2021 to dec 2024: 
cp ../data_new/glorys_1-12deg/cmems_mlt_interim_until_dec_2024.nc . 

# 1.1 Concatenate mlt 1/12deg in time and remove useless tmp data
cdo mergetime cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc cmems_mlt_interim_until_dec_2024.nc mlt_1993-dec2024.nc

rm cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc
rm cmems_mlt_interim_until_dec_2024.nc

# 1.3 Convert to F64 precision and remove useless tmp data
 cdo -b 64 copy mlt_1993-dec2024.nc mlt_1993-dec2024_F64.nc
 rm mlt_1993-dec2024.nc

# 1.4 Mask-out Pacific Ocean grid points and remove useless tmp data
 cdo maskregion,mymasknew.dat mlt_1993-dec2024_F64.nc mlt_1993-dec2024_F64_masked_tmp.nc
 rm mlt_1993-dec2024_F64.nc

# AND FOCUS ON THE LON LAT AREA SPECIFIED HERE, FOR THE CHANGES CALCULATION
  cdo sellonlatbox,-90,0,0,20 mlt_1993-dec2024_F64_masked_tmp.nc mlt_1993-dec2024_F64_masked.nc

 rm mlt_1993-dec2024_F64_masked_tmp.nc 


# 1.5 Prepare MLT timeseries for the NL Regression Model.

# At the ORIGINAL resolution for the whole available period, subsect a smaller area to retrieve the time series necessary for the NL regression:
# The non linear regression timeseries are spatial means of deseasonalized (point-by-point) data at high resolution
# over the follwing area

         
        # Deseasonalize the 1993-2022 part and the 2023-2024 part separately, point by point. 
	# This is because we want to do prediction only on the 2023-2024 part, and the regr model only on the 1993-2022 part 

        #cdo -L -ymonsub mlt_1993-dec2024_F64_masked.nc -ymonmean mlt_1993-dec2024_F64_masked.nc mlt_1993-dec2024_F64_masked_des.nc


	# SELECT PART 1 (JAN 1993-DEC 2022) AND DESEASONALIZE THIS PART ALONE (remove seasonal cycle of 1993-2022 only)
        
        cdo selyear,1993/2022 mlt_1993-dec2024_F64_masked.nc mlt_1993-dec2022_F64_masked.nc 
        cdo -L ymonsub mlt_1993-dec2022_F64_masked.nc -ymonmean mlt_1993-dec2022_F64_masked.nc mlt_1993-dec2022_F64_masked_deseas.nc


        # Store the 1993-2022 ymonmean climatology in a file 

        cdo ymonmean mlt_1993-dec2022_F64_masked.nc ymonmean_of_mlt_1993-dec2022_F64_masked.nc  #This will have only 12 timesteps, which I will have to repeat and subctratc month by month in the following 2023-2024 period 

        # Now repeat the ymonmean timesteps to reach the same length of the monthly mlt data over 2023-2024 (i.e. 24 timesteps in total)
        cdo copy ymonmean_of_mlt_1993-dec2022_F64_masked.nc ymonmean_repeat.nc

        for i in $(seq 1 )
        do
        shift="${i}year"
        cdo -shifttime,${shift} ymonmean_of_mlt_1993-dec2022_F64_masked.nc next_year.nc
        cdo -cat ymonmean_repeat.nc next_year.nc tmp.nc
        mv tmp.nc ymonmean_repeat.nc
        done

        # SELECT PART 2 (JAN 2023-DEC 2024) AND DESEASONALIZE THIS PART ALONE (remove seasonal cycle of 2023-2024 only)

        cdo selyear,2023/2024 mlt_1993-dec2024_F64_masked.nc mlt_2023-dec2024_F64_masked.nc
        cdo sub mlt_2023-dec2024_F64_masked.nc ymonmean_repeat.nc mlt_2023-dec2024_F64_masked_deseas.nc


        # RECOMPOSE TEH FIELD, I.E. MERGE IN TIME THE DESEASONALIZED 1993-2022 PART AND THE DESEASONALIZED 2023-2024 PART. THESE ARE STILL MAPS. 
        cdo mergetime  mlt_1993-dec2022_F64_masked_deseas.nc  mlt_2023-dec2024_F64_masked_deseas.nc  mlt_1993-dec2024_F64_masked_splitted-deseas.nc # <=== deseas point-by-point map, deseasonalized in chunks

	# Subset smaller area OF THE RECOMPOSED FIELD (PART 1 + PART2): select ONLY the points within the Bloom Contour obtained by the digitization
        
    cp ../utils/mask_bloom.dat .

        cdo maskregion,mask_bloom.dat mlt_1993-dec2024_F64_masked_splitted-deseas.nc mlt_1993-dec2024_F64_masked_des_BLOOM-AREA.nc 

        # Space average over smaller area: THE RESULTING TIMESERIES IS THE INPUT MLT FOR THE NL REGRESSION MODEL:
        cdo fldmean mlt_1993-dec2024_F64_masked_des_BLOOM-AREA.nc fldm_of_1993_dec2024_glorys_mlt_BLOOM-AREA_masked_deseas_1-12deg.nc # <=== This is for NL regression model, if extending mlt until Dec 2024


        # remove useless stuff
        rm  mlt_1993-dec2024_F64_masked_des.nc
	rm  mlt_1993-dec2024_F64_masked.nc
#	rm  mlt_1993-dec2024_F64_masked_des_reducedGASB.nc
    

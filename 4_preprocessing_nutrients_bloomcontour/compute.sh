#!/bin/bash

# ------------------------------------------------
# 1) Preprocessing of MLT data:
# ------------------------------------------------

# 1. mlt whole period Jan 1993- Nov 2024: do a little initial preprocessing on the mlt data: 

# Copy data at original resolution from there to here
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc .
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc .

# 1.1 Concatenate mlt 1/12deg in time and remove useless tmp data
cdo mergetime cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc mlt_1-12deg_1993_01-2024_11.nc

rm cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc
rm cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc

# 1.3 Convert to F64 precision and remove useless tmp data
 cdo -b 64 copy mlt_1-12deg_1993_01-2024_11.nc mlt_1-12deg_1993_01-2024_11_F64.nc
 rm mlt_1-12deg_1993_01-2024_11.nc

# 1.4 Mask-out Pacific Ocean grid points and remove useless tmp data
 cdo maskregion,mymasknew.dat mlt_1-12deg_1993_01-2024_11_F64.nc mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc
 rm mlt_1-12deg_1993_01-2024_11_F64.nc

# AND FOCUS ON THE LON LAT AREA SPECIFIED HERE, FOR THE CHANGES CALCULATION 
  cdo sellonlatbox,-90,0,0,20 mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc mlt_1-12deg_1993_01-2024_11_F64_masked.nc

 rm mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc


# 1.5 Remap at 0.25deg to get rid of submesoscales for the changes calculation
cdo remapbil,grid025new_smaller.dat mlt_1-12deg_1993_01-2024_11_F64_masked.nc mlt_025deg_1993_01-2024_11_F64_masked.nc


# 1.6 Select from January 1993 to December 2022 only 
#cdo selyear,1993/2022 mlt_025deg_1993_01-2024_11_F64_masked.nc mlt_1993-2022.nc
#rm mlt_1-12deg_1993_01-2024_11_F64_masked.nc 
#rm mlt_025deg_1993_01-2024_11_F64_masked.nc 


# 1.7 SELECT ONLY MONTHLY MLT DATA AFTER BLOOM, I.E. FROM 2011 TO 2022
cdo selyear,2011/2022 mlt_025deg_1993_01-2024_11_F64_masked.nc MLT_2011-2022_f64_025_masked.nc # <--- This is one of the files to use in matlab 

# 1.8 COMPUTE MONTHLY CLIMATOLOGY OVER ONLY THE THE PRE-BLOOM PERIOD, I.E. 1999-2010 
cdo selyear,1999/2010 mlt_025deg_1993_01-2024_11_F64_masked.nc mlt_1999-2010.nc 
cdo ymonmean mlt_1999-2010.nc ymonmean_MLT_1999-2010_f64_025_masked.nc # this will have only 12 timesteps 


# 1.9 Finally, repeat the ymonmean timesteps to reach the same length of the monthly mlt data after bloom (for calcultion reasons, see matlab)

cdo copy ymonmean_MLT_1999-2010_f64_025_masked.nc ymonmean_mlt_1999-2010_f64_025_masked_repeat.nc

for i in $(seq 1 11)
do
  shift="${i}year"
  cdo -shifttime,${shift} ymonmean_MLT_1999-2010_f64_025_masked.nc next_year.nc
  cdo -cat ymonmean_mlt_1999-2010_f64_025_masked_repeat.nc next_year.nc tmp.nc
  mv tmp.nc ymonmean_mlt_1999-2010_f64_025_masked_repeat.nc
done


rm ymonmean_MLT_1999-2010_f64_025_masked.nc

rm next_year.nc

# Remove useless stuff 
rm mlt_025deg_1993_01-2024_11_F64_masked.nc 
rm mlt_1-12deg_1993_01-2024_11_F64_masked.nc 
rm mlt_1993-2022.nc 
rm mlt_1999-2010.nc



# ------ End of MLT preprocessing -----

# ------------------------------------
# 2) Preprocessing of NUTRIENTS data 
# ------------------------------------

# 2.1 copy data from there to here, put them to f64 precision. 
cp ../data_nutrients/cmems_mod_glo_bgc_my_0.25_P1M-m_no3-po4_90.00W-0.00E_0.00N-25.00N_0.51-163.16m_1993-01-01-2022-12-01.nc .

cdo -b 64 copy cmems_mod_glo_bgc_my_0.25_P1M-m_no3-po4_90.00W-0.00E_0.00N-25.00N_0.51-163.16m_1993-01-01-2022-12-01.nc nutri_1993-2022_f64.nc

# rm useless 
rm cmems_mod_glo_bgc_my_0.25_P1M-m_no3-po4_90.00W-0.00E_0.00N-25.00N_0.51-163.16m_1993-01-01-2022-12-01.nc

# 2.2. make sure th grid is teh same, i.e. remap just to have teh same structure. This will not change the resolution because it's already at 0.25. But maybe different netCDF have different data structure, attributes ecc. That's why I want to have everything on teh same 0.25 grid. 


cdo remapbil,grid025new_smaller.dat nutri_1993-2022_f64.nc nutri_1993-2022_f64_025.nc

rm nutri_1993-2022_f64.nc

# 2.3 Mask out points in the Pacific Ocean 
cdo maskregion,mymasknew.dat nutri_1993-2022_f64_025.nc nutri_1993-2022_f64_025_masked.nc

# 2.4 Select monthly nutrients before bloom
cdo selyear,1999/2010 nutri_1993-2022_f64_025_masked.nc nutri_1999-2010.nc

rm nutri_1993-2022_f64_025_masked.nc

# 2.5 NUTRIENT MONTHLY CLIMATOLOGY BEFORE BLOOM

cdo ymonmean nutri_1999-2010.nc ymonmean_NUTRI_1999-2010.nc 

# 2.6 Finally, repeat the same 12 months over and over to reach the same length as the mlt monthly data (for calculation reasons, see matlab script) 

cdo copy ymonmean_NUTRI_1999-2010.nc ymonmean_NUTRI_1999-2010_repeat.nc

for i in $(seq 1 11)
do
  shift="${i}year"
  cdo -shifttime,${shift} ymonmean_NUTRI_1999-2010.nc next_year.nc
  cdo -cat ymonmean_NUTRI_1999-2010_repeat.nc next_year.nc tmp.nc
  mv tmp.nc ymonmean_NUTRI_1999-2010_repeat.nc
done

# Remove useless stuff
rm ymonmean_NUTRI_1999-2010.nc
rm next_year.nc
rm nutri_1993-2022_f64_025.nc 
rm nutri_1999-2010.nc

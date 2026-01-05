#!/bin/bash


# Note: In this script the dust is at its native resolution: 0.75°. 
# Everything else is remapped at 0.25° to get rid of small scale structures


# ------------------------------------------------
# Part A: Preliminary preprocessing of MLT data:
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

# 1.5 Prepare MLT timeseries for the NL Regression Model. 

# At the ORIGINAL resolution for the whole available period, subsect a smaller area to retrieve the time series necessary for the NL regression: 
# The non linear regression timeseries are spatial means of deseasonalized (point-by-point) data at high resolution
# over the follwing area 

	# Select until Dec 2022 - the latest availabe date of the NL regression model
	cdo selyear,1993/2022 mlt_1-12deg_1993_01-2024_11_F64_masked.nc mlt_1-12deg_1993-2022_F64_masked.nc
	# Deseas point by point
	cdo -L -ymonsub mlt_1-12deg_1993-2022_F64_masked.nc -ymonmean mlt_1-12deg_1993-2022_F64_masked.nc mlt_1-12deg_1993-2022_F64_masked_des.nc
	# Subset smaller area 
	cdo sellonlatbox,-89,-15,1,15 mlt_1-12deg_1993-2022_F64_masked_des.nc mlt_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc
	
	# Space average over smaller area: THE RESULTING TIMESERIES IS THE INPUT MLT FOR THE NL REGRESSION MODEL: 
	cdo fldmean mlt_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc fldm_of_1993_2022_glorys_mlt_reducedGASB_masked_deseas_1-12deg_Lat_1-15.nc # <=== In case you want the NonLinear Regression model over that box and at 1/12°.

	# remove useless stuff 
	rm mlt_1-12deg_1993-2022_F64_masked_des.nc
	rm mlt_1-12deg_1993-2022_11_F64_masked.nc


# 1.6 Remap at 0.25deg to get rid of submesoscales for the changes calculation
cdo remapbil,grid025new_smaller.dat mlt_1-12deg_1993_01-2024_11_F64_masked.nc mlt_025deg_1993_01-2024_11_F64_masked.nc

# ------ End of MLT preliminary preprocessing -----



# ------------------------------------------------
# Part B: Preliminary preprocessing of SST data:
# ------------------------------------------------

# 2. thetao whole period Jan 1993- Nov 2024: do a little initial preprocessing on the sst data:

# Copy data at original resolution from there to here
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc .
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc .

# 2.1 Concatenate sst 1/12deg in time and remove useless tmp data
cdo mergetime cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc sst_1-12deg_1993_01-2024_11.nc

rm cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc 
rm cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc

# 2.3 Convert to F64 precision and remove useless tmp data
 cdo -b 64 copy sst_1-12deg_1993_01-2024_11.nc sst_1-12deg_1993_01-2024_11_F64.nc
 rm sst_1-12deg_1993_01-2024_11.nc

# 2.4 Mask-out Pacific Ocean grid points and remove useless tmp data
 cdo maskregion,mymasknew.dat sst_1-12deg_1993_01-2024_11_F64.nc sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc
 rm sst_1-12deg_1993_01-2024_11_F64.nc

# AND FOCUS ON THE LON LAT AREA SPECIFIED HERE, FOR THE CHANGES CALCULATION 
cdo sellonlatbox,-90,0,0,20 sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc sst_1-12deg_1993_01-2024_11_F64_masked.nc 
rm sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

# 2.5: Prepare SST timeseries for the NL Regression Model
 
# At the ORIGINAL resolution for the whole available period, subsect a smaller area to retrieve the time series necessary for the NL regression:
# The non linear regression timeseries are spatial means of deseasonalized (point-by-point) data at high resolution
# over the follwing area

        # Select until Dec 2022 - the latest availabe date of the NL regression model
        cdo selyear,1993/2022 sst_1-12deg_1993_01-2024_11_F64_masked.nc sst_1-12deg_1993-2022_F64_masked.nc
        # Deseas point by point
        cdo -L -ymonsub sst_1-12deg_1993-2022_F64_masked.nc -ymonmean sst_1-12deg_1993-2022_F64_masked.nc sst_1-12deg_1993-2022_F64_masked_des.nc
        # Subset smaller area                 
	cdo sellonlatbox,-89,-15,1,15 sst_1-12deg_1993-2022_F64_masked_des.nc sst_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc                                                           

	# Space average over smaller area: THE RESULTING TIMESERIES IS THE INPUT MLT FOR THE NL REGRESSION MODEL:      
	cdo fldmean sst_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc fldm_of_1993_2022_glorys_sst_reducedGASB_masked_deseas_1-12deg_Lat_1-15.nc # <=== In case you want the NonLinear Regression model over a box -89,-15,5,25 and at 1/12°.

        # remove useless stuff
        rm sst_1-12deg_1993-2022_F64_masked_des.nc
        rm sst_1-12deg_1993-2022_11_F64_masked.nc


# 2.6 Remap at 0.25deg to get rid of submesoscales for the changes calculations
cdo remapbil,grid025new_smaller.dat sst_1-12deg_1993_01-2024_11_F64_masked.nc sst_025deg_1993_01-2024_11_F64_masked.nc

# ------ End of SST preliminary preprocessing -----

# ------------------------------------------------
# Part C1: Preliminary preprocessing of SALINITY data:
# ------------------------------------------------

# 1. Salinity whole period Jan 1993 - Nov 2024: do a little initial preprocessing on the mlt data: 

# Copy data at original resolution from there to here
cp ../data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1765297853330.nc .

# 1.1 no need to concatenate, already whole period. 

# 1.2 Convert to F64 precision and remove useless tmp data

 cdo -b 64 copy cmems_mod_glo_phy_my_0.083deg_P1M-m_1765297853330.nc sal_1-12deg_1993_01-2024_11_F64.nc
 rm cmems_mod_glo_phy_my_0.083deg_P1M-m_1765297853330.nc

# 1.4 Mask-out Pacific Ocean grid points and remove useless tmp data

 cdo maskregion,mymasknew.dat sal_1-12deg_1993_01-2024_11_F64.nc sal_1-12deg_1993_01-2024_11_F64_masked_tmp.nc
 rm sal_1-12deg_1993_01-2024_11_F64.nc

# AND FOCUS ON THE LON LAT AREA SPECIFIED HERE, FOR THE CHANGES CALCULATION 
 cdo sellonlatbox,-90,0,0,20 sal_1-12deg_1993_01-2024_11_F64_masked_tmp.nc sal_1-12deg_1993_01-2024_11_F64_masked.nc

 rm sal_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

# 1.5 Prepare SAL timeseries for the NL Regression Model - Review section. 

# At the ORIGINAL resolution for the whole available period, subset a smaller area to retrieve the time series necessary for the NL regression: 
# The non linear regression timeseries are spatial means of deseasonalized (point-by-point) data at high resolution
# over the following area 

        # Select until Dec 2022 - the latest available date of the NL regression model
        cdo selyear,1993/2022 sal_1-12deg_1993_01-2024_11_F64_masked.nc sal_1-12deg_1993-2022_F64_masked.nc

        # Deseas point by point
        cdo -L -ymonsub sal_1-12deg_1993-2022_F64_masked.nc -ymonmean sal_1-12deg_1993-2022_F64_masked.nc sal_1-12deg_1993-2022_F64_masked_des.nc

        # Subset smaller area 
        cdo sellonlatbox,-89,-15,1,15 sal_1-12deg_1993-2022_F64_masked_des.nc sal_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc

        # Space average over smaller area: THE RESULTING TIME SERIES IS THE INPUT SAL FOR THE NL REGRESSION MODEL (REVISION): 
        cdo fldmean sal_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc fldm_of_1993_2022_glorys_sal_reducedGASB_masked_deseas_1-12deg_Lat_1-15.nc # <=== In case you want the NonLinear Regression model over that box and at 1/12°.

        # remove useless stuff 
        rm sal_1-12deg_1993-2022_F64_masked_des.nc
        rm sal_1-12deg_1993-2022_11_F64_masked.nc


# 1.6 Remap at 0.25deg to get rid of submesoscales for the changes calculation
cdo remapbil,grid025new_smaller.dat sal_1-12deg_1993_01-2024_11_F64_masked.nc sal_025deg_1993_01-2024_11_F64_masked.nc

# --- End of Salinity preprocessing --- 

# ------------------------------------------------------
# Part C2: Preliminary preprocessing of Dust data
# ------------------------------------------------------

# 3.1 Dust original resolution is 0.75°: copy data from there to here
cp ../data_dust/duaod550_f64_res075_2003-2022.nc .

# 3.2 Select a wider area aroung the GASB just for the dust data that need it (GoM not included), 
# as the mlt and sst data are already over that area - this is for the changes calculations. 
cdo sellonlatbox,-90,0,0,20 duaod550_f64_res075_2003-2022.nc duaod550_f64_res075_2003-2022_GASBtmp.nc

# 3.3. Increase dust data precision to F64 
cdo -b 64 -copy duaod550_f64_res075_2003-2022_GASBtmp.nc duaod550_f64_res075_2003-2022_GASB_unmasked.nc

# 3.4 Mask out the Pacific Ocean points 

cdo maskregion,mymasknew.dat duaod550_f64_res075_2003-2022_GASB_unmasked.nc duaod550_f64_res075_2003-2022_GASB.nc
rm duaod550_f64_res075_2003-2022_GASB_unmasked.nc

# 3.5 Prepare DUST timeseries for NL Regression Model: 

# At the ORIGINAL resolution for the whole available period, subsect a smaller area to retrieve the time series necessary for the NL regression:
# The non linear regression timeseries are spatial means of deseasonalized (point-by-point) data at high resolution
# over the follwing area

	# Remove seasonal cycle 
        cdo -L -ymonsub duaod550_f64_res075_2003-2022_GASB.nc -ymonmean duaod550_f64_res075_2003-2022_GASB.nc duaod550_f64_res075_2003-2022_GASB_deseas.nc
	# Select Area 
        cdo sellonlatbox,-89,-15,1,15 duaod550_f64_res075_2003-2022_GASB_deseas.nc duaod550_f64_res075_2003-2022_reducedGASB_deseas_Lat_1-15.nc	
	
	rm duaod550_f64_res075_2003-2022_GASB_deseas.nc

        # Get spatial average: the resultimg timeseries is the dust input for the NL model
        cdo fldmean duaod550_f64_res075_2003-2022_reducedGASB_deseas_Lat_1-15.nc fldm_of_duaod550_f64_res075_2003-2022_deseas_reducedGASB_Lat_1-15.nc # <=== NL model input for dust
        rm duaod550_f64_res075_2003-2022_reducedGASB_deseas.nc

# Uncomment to remove useless data to save space
#rm *remapped.nc
#rm *masked.nc
#rm *masked.nc


# ------------------------------------------------
# Part D: Multi year seasonal means calculation 
# ------------------------------------------------

# 4.1 Remove seasonal cycle point by point from dust, mlt and sst data
cdo -L -ymonsub duaod550_f64_res075_2003-2022_GASB.nc -ymonmean duaod550_f64_res075_2003-2022_GASB.nc duaod550_f64_res075_2003-2022_GASB_deseas.nc

cdo -L -ymonsub mlt_025deg_1993_01-2024_11_F64_masked.nc -ymonmean mlt_025deg_1993_01-2024_11_F64_masked.nc mlt_025deg_1993_01-2024_11_F64_masked_deseas.nc

cdo -L -ymonsub sst_025deg_1993_01-2024_11_F64_masked.nc -ymonmean sst_025deg_1993_01-2024_11_F64_masked.nc sst_025deg_1993_01-2024_11_F64_masked_deseas.nc

mv ./*deseas.nc ./tmpfiles2/

# -----------------------------------------------------------------------
#  PART E: Additional calculations: 
#  TREND MAPS OF DESEASONALIZED FIELDS IN SUBPERIODS AND WHOLE PERIOD
# -----------------------------------------------------------------------

export SKIP_SAME_TIME=1


# ================
# Period BEFORE
# ================

cdo selyear,1999/2010 mlt_025deg_1993_01-2024_11_F64_masked.nc mlt_025deg_1999-2010_F64_masked.nc
cdo selyear,1999/2010 sst_025deg_1993_01-2024_11_F64_masked.nc sst_025deg_1999-2010_F64_masked.nc
cdo selyear,1999/2010 sal_025deg_1993_01-2024_11_F64_masked.nc sal_025deg_1999-2010_F64_masked.nc

# Dust are available only from 2003. Before for dust is 2003-2011 (or 2003-2010, up to the user- but make sure "before" and "after" have the same length)
cdo selyear,2003/2011 duaod550_f64_res075_2003-2022_GASB.nc duaod550_f64_res075_2003-2011_GASB.nc

# Remove seasonal cycle computed over just this subperiod:

cdo -L -ymonsub duaod550_f64_res075_2003-2011_GASB.nc -ymonmean duaod550_f64_res075_2003-2011_GASB.nc duaod550_f64_res075_2003-2011_GASB_des.nc

cdo -L -ymonsub mlt_025deg_1999-2010_F64_masked.nc -ymonmean mlt_025deg_1999-2010_F64_masked.nc mlt_025deg_1999-2010_F64_masked_des.nc
cdo -L -ymonsub sst_025deg_1999-2010_F64_masked.nc -ymonmean sst_025deg_1999-2010_F64_masked.nc sst_025deg_1999-2010_F64_masked_des.nc


# Compute trend maps of deseasonalized field, in this subperiod:

cdo trend duaod550_f64_res075_2003-2011_GASB_des.nc a trendmap_duaod550_2003-2011.nc
cdo trend mlt_025deg_1999-2010_F64_masked_des.nc a trendmap_mlt_1999-2010.nc
cdo trend sst_025deg_1999-2010_F64_masked_des.nc a trendmap_sst_1999-2010.nc
rm a

# ==============
# Period AFTER
# ==============

# "After" for Dust is different (less data available and "before" and "after" must have the same length for a same variable)
cdo selyear,2011/2019 duaod550_f64_res075_2003-2022_GASB.nc duaod550_f64_res075_2011-2019_GASB.nc
cdo selyear,2011/2022 duaod550_f64_res075_2003-2022_GASB.nc duaod550_f64_res075_2011-2022_GASB.nc
cdo selyear,2011/2022 mlt_025deg_1993_01-2024_11_F64_masked.nc mlt_025deg_2011-2022_F64_masked.nc
cdo selyear,2011/2022 sst_025deg_1993_01-2024_11_F64_masked.nc sst_025deg_2011-2022_F64_masked.nc
cdo selyear,2011/2022 sal_025deg_1993_01-2024_11_F64_masked.nc sal_025deg_2011-2022_F64_masked.nc

# Remove seasonal cycle computed over this subperiod: 
cdo -L -ymonsub duaod550_f64_res075_2011-2019_GASB.nc -ymonmean duaod550_f64_res075_2011-2019_GASB.nc duaod550_f64_res075_2011-2019_GASB_des.nc

cdo -L -ymonsub mlt_025deg_2011-2022_F64_masked.nc -ymonmean mlt_025deg_2011-2022_F64_masked.nc mlt_025deg_2011-2022_F64_masked_des.nc
cdo -L -ymonsub sst_025deg_2011-2022_F64_masked.nc -ymonmean sst_025deg_2011-2022_F64_masked.nc sst_025deg_2011-2022_F64_masked_des.nc

# Compute trend maps of deseasonalized fields, in period 3:
cdo trend duaod550_f64_res075_2011-2019_GASB_des.nc a trendmap_duaod550_2011-2019.nc
cdo trend mlt_025deg_2011-2022_F64_masked_des.nc a trendmap_mlt_2011-2022.nc
cdo trend sst_025deg_2011-2022_F64_masked_des.nc a trendmap_sst_2011-2022.nc
rm a 

# *******************************************
# PART F: SIMPLIFIED COMPARISON OF CHANGES
# *******************************************

# ========================================
#       After vs Before
# ========================================

# ------------------------------------------
#  Changes in multi-year seasonal means 
# ------------------------------------------


# Compute multi-year seasonal means 

# before:
cdo yseasmean duaod550_f64_res075_2003-2011_GASB.nc yseasm_2003-2011_dust.nc 
cdo yseasmean mlt_025deg_1999-2010_F64_masked.nc yseasm_1999-2010_mlt.nc
cdo yseasmean sst_025deg_1999-2010_F64_masked.nc yseasm_1999-2010_sst.nc
cdo yseasmean sal_025deg_1999-2010_F64_masked.nc yseasm_1999-2010_sal.nc

# provide also a version in Kelvin for yseasm SST (so that we can compute changes in absolute units too - useful if interested in computing percentages in the postprocessing part later on) 
cdo addc,273.15 yseasm_1999-2010_sst.nc yseasm_1999-2010_sst_Kelvin.nc

# after:
cdo yseasmean duaod550_f64_res075_2011-2019_GASB.nc yseasm_2011-2019_dust.nc 
cdo yseasmean duaod550_f64_res075_2011-2022_GASB.nc yseasm_2011-2022_dust.nc 
cdo yseasmean mlt_025deg_2011-2022_F64_masked.nc yseasm_2011-2022_mlt.nc
cdo yseasmean sst_025deg_2011-2022_F64_masked.nc yseasm_2011-2022_sst.nc
cdo yseasmean sal_025deg_2011-2022_F64_masked.nc yseasm_2011-2022_sal.nc

# provide also a version in Kelvin for yseasm SST (so that we can compute changes in absolute units too - useful if interested in computing percentages in the postprocessing part later on)
cdo addc,273.15 yseasm_2011-2022_sst.nc yseasm_2011-2022_sst_Kelvin.nc

# Compute yseasm after - yseasm before: 

cdo sub yseasm_2011-2019_dust.nc  yseasm_2003-2011_dust.nc yseasm_changes_2011-2019_vs_2003-2011_dust.nc
cdo sub yseasm_2011-2022_mlt.nc  yseasm_1999-2010_mlt.nc yseasm_changes_2011-2022_vs_1999-2010_mlt.nc
cdo sub yseasm_2011-2022_sst.nc  yseasm_1999-2010_sst.nc yseasm_changes_2011-2022_vs_1999-2010_sst.nc
cdo sub yseasm_2011-2022_sst_Kelvin.nc  yseasm_1999-2010_sst_Kelvin.nc yseasm_changes_2011-2022_vs_1999-2010_sst_Kelvin.nc
cdo sub yseasm_2011-2022_sal.nc  yseasm_1999-2010_sal.nc yseasm_changes_2011-2022_vs_1999-2010_sal.nc

#mv yseasm* ./tmp_files/
#mv ./tmp_files/yseasm_changes* . 

mv trend* tmp_files/
mv tmp_files/trend_changes* . 
mv tmp_files/trendmap* .
mv tmp_files/fldm* . 
mv tmp_files/fldsum . 

#rm ./tmp_files/*
rm *tmp.nc


mv ./tmpfiles2/* .




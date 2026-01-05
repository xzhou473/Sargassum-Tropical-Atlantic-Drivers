#!/bin/bash
set -e

#### This script is used to extract different domains
lon1=-55
lon2=-20
lat1=0
lat2=10
outputname=ITCZ

rm -f *.nc

echo "------------------------------------------------"
echo "# Part A: Preliminary preprocessing of MLT data"
echo "------------------------------------------------"

# 1. Copy original data
ln -s /localdata1/shared_ljuba_Sargassum/TropicalAtlantic_to_share_with_Regr_Lat_1_15/data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc .
ln -s /localdata1/shared_ljuba_Sargassum/TropicalAtlantic_to_share_with_Regr_Lat_1_15/data_new/glorys_1-12deg/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc .

# 1.1 Merge
cdo mergetime \
    cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc \
    cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc \
    mlt_1-12deg_1993_01-2024_11.nc

rm -f cmems_mod_glo_phy_my_0.083deg_P1M-m_1737743544477.nc
rm -f cmems_mod_glo_phy_myint_0.083deg_P1M-m_1740598002242.nc

# 1.3 Convert precision
cdo -b 64 copy mlt_1-12deg_1993_01-2024_11.nc mlt_1-12deg_1993_01-2024_11_F64.nc
rm -f mlt_1-12deg_1993_01-2024_11.nc

# 1.4 Mask-out Pacific
cdo maskregion,mymasknew.dat \
    mlt_1-12deg_1993_01-2024_11_F64.nc \
    mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

rm -f mlt_1-12deg_1993_01-2024_11_F64.nc

# Restrict domain (0–20N, -90–0W)
cdo sellonlatbox,-90,0,0,20 \
    mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc \
    mlt_1-12deg_1993_01-2024_11_F64_masked.nc

rm -f mlt_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

# 1.5 Prepare MLT for regression model
cdo selyear,1993/2022 \
    mlt_1-12deg_1993_01-2024_11_F64_masked.nc \
    mlt_1-12deg_1993-2022_F64_masked.nc

# Deseasonalize
cdo -L -ymonsub \
    mlt_1-12deg_1993-2022_F64_masked.nc \
    -ymonmean mlt_1-12deg_1993-2022_F64_masked.nc \
    mlt_1-12deg_1993-2022_F64_masked_des.nc

# Extract smaller domain
cdo sellonlatbox,$lon1,$lon2,$lat1,$lat2 \
    mlt_1-12deg_1993-2022_F64_masked_des.nc \
    mlt_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc

# Spatial mean
cdo fldmean \
    mlt_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc \
    fldm_of_1993_2022_glorys_mlt_reducedGASB_masked_deseas_$outputname.nc

rm -f mlt_1-12deg_1993-2022_F64_masked_des.nc
rm -f mlt_1-12deg_1993-2022_F64_masked.nc

echo "------------------------------------------------"
echo "# Part B: Preliminary preprocessing of SST data"
echo "------------------------------------------------"

# 2. Copy SST data
ln -s /localdata1/shared_ljuba_Sargassum/TropicalAtlantic_to_share_with_Regr_Lat_1_15/data_new/glorys_1-12deg/cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc .
ln -s /localdata1/shared_ljuba_Sargassum/TropicalAtlantic_to_share_with_Regr_Lat_1_15/data_new/glorys_1-12deg/cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc .

# 2.1 Merge
cdo mergetime \
    cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc \
    cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc \
    sst_1-12deg_1993_01-2024_11.nc

rm -f cmems_mod_glo_phy_my_0.083deg_P1M-m_1739994093893.nc
rm -f cmems_mod_glo_phy_myint_0.083deg_P1M-m_1739994256211.nc

# Convert precision
cdo -b 64 copy sst_1-12deg_1993_01-2024_11.nc sst_1-12deg_1993_01-2024_11_F64.nc
rm -f sst_1-12deg_1993_01-2024_11.nc

# Mask Pacific
cdo maskregion,mymasknew.dat \
    sst_1-12deg_1993_01-2024_11_F64.nc \
    sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

rm -f sst_1-12deg_1993_01-2024_11_F64.nc

# Restrict domain
cdo sellonlatbox,-90,0,0,20 \
    sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc \
    sst_1-12deg_1993_01-2024_11_F64_masked.nc

rm -f sst_1-12deg_1993_01-2024_11_F64_masked_tmp.nc

# Select 1993–2022
cdo selyear,1993/2022 \
    sst_1-12deg_1993_01-2024_11_F64_masked.nc \
    sst_1-12deg_1993-2022_F64_masked.nc

# Deseasonalize
cdo -L -ymonsub \
    sst_1-12deg_1993-2022_F64_masked.nc \
    -ymonmean sst_1-12deg_1993-2022_F64_masked.nc \
    sst_1-12deg_1993-2022_F64_masked_des.nc

# Extract box
cdo sellonlatbox,$lon1,$lon2,$lat1,$lat2 \
    sst_1-12deg_1993-2022_F64_masked_des.nc \
    sst_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc

# Spatial mean
cdo fldmean \
    sst_1-12deg_1993-2022_F64_masked_des_reducedGABS_Lat_1-15.nc \
    fldm_of_1993_2022_glorys_sst_reducedGASB_masked_deseas_$outputname.nc

rm -f sst_1-12deg_1993-2022_F64_masked_des.nc
rm -f sst_1-12deg_1993-2022_F64_masked.nc

echo "------------------------------------------------"
echo "# Part C: Dust preprocessing"
echo "------------------------------------------------"

# Dust data
ln -s /localdata1/shared_ljuba_Sargassum/TropicalAtlantic_to_share_with_Regr_Lat_1_15/data_dust/duaod550_f64_res075_2003-2022.nc .

# Subset domain
cdo sellonlatbox,-90,0,0,20 \
    duaod550_f64_res075_2003-2022.nc \
    duaod550_f64_res075_2003-2022_GASBtmp.nc

# Convert precision
cdo -b 64 copy \
    duaod550_f64_res075_2003-2022_GASBtmp.nc \
    duaod550_f64_res075_2003-2022_GASB_unmasked.nc

rm -f duaod550_f64_res075_2003-2022_GASBtmp.nc

# Mask Pacific
cdo maskregion,mymasknew.dat \
    duaod550_f64_res075_2003-2022_GASB_unmasked.nc \
    duaod550_f64_res075_2003-2022_GASB.nc

rm -f duaod550_f64_res075_2003-2022_GASB_unmasked.nc

# Deseasonalize
cdo -L -ymonsub \
    duaod550_f64_res075_2003-2022_GASB.nc \
    -ymonmean duaod550_f64_res075_2003-2022_GASB.nc \
    duaod550_f64_res075_2003-2022_GASB_deseas.nc

# Extract box
cdo sellonlatbox,$lon1,$lon2,$lat1,$lat2 \
    duaod550_f64_res075_2003-2022_GASB_deseas.nc \
    duaod550_f64_res075_2003-2022_reducedGASB_deseas_Lat_1-15.nc

rm -f duaod550_f64_res075_2003-2022_GASB_deseas.nc

# Spatial mean
cdo fldmean \
    duaod550_f64_res075_2003-2022_reducedGASB_deseas_Lat_1-15.nc \
    fldm_of_duaod550_f64_res075_2003-2022_deseas_reducedGASB_$outputname.nc


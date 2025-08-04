#!/bin/bash

# Comment/Uncomment based on which version of 2. you want to run. Do not deactivate 1. 

rm *.nc 
cp ../utils/grid025new_smaller.dat .

./1_compute_yseasm_changes_NEWDATA_bloomcontour.sh
./1-1_compute_longer_mlt_timeser_bloomcontour.sh
#./2_compute_cumulative_changes_NEWDATA_withoutSST_smallerbox.sh
#./2_compute_cumulative_changes_NEWDATA_withSSTincr_smallerbox.sh


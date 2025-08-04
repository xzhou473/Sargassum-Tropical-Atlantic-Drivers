#!/bin/bash

rm *.nc

./1_compute_trend.sh  
./2_compute_changes.sh


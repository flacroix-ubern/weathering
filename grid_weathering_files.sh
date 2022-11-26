#!/bin/sh
#load modules
#load module python
#load module cdo
# Define inputs

ECHAM_FILE="echam_inputs/echam_test_data.nc"
ECHAM_AREA="echam_inputs/echam_area.nc"
#Run weathering model
python weathering.py $ECHAM_FILE $ECHAM_AREA

# Regrid P and Alk release from ECHAM to HD model grid for calculations of transformations
cdo setunit,"g P m-2 yr-1" -setgrid,$ECHAM_FILE test_prelease.nc prelease_wthgrd.nc
cdo setunit,"g C m-2 yr-1 of HCO3-" -setgrid,$ECHAM_FILE test_alkrelease.nc alkrelease_wthgrd.nc
cdo setunit,"g C m-2 yr-1 of CO2" -setgrid,$ECHAM_FILE test_co2drawdown.nc co2drawdown_wthgrd.nc
cdo gridarea co2drawdown_wthgrd.nc weathering_gridarea.nc


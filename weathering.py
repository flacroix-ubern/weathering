#Compute global weathering
# Script by Fabrice Lacroix
# 10.10.2019
# Updated 23/04/20 
# Clean up 15/11/20 
#
# The script computes the global p weathering release derived using a spatiotemporal model as described in Hartmann et al. 2014
#  and Lacroix et al. (2020).

#  The P weathering release (Pr) is calculated:                
#     Pr = b * q *f_t * s * glac_covr
#     where
#     q: annual runoff [mm/a]
#     s: factor to correct for soil shielding
#     b: empirical lithology dependent factors 
#     glac_covr: Glacier coverage [-]
#     temperature response term f_t is given by:
#     f_t = exp(-Ea/R*(1/T-1/Tref))
#     with
#     Ea:  apparent activation energy (depends on lithology)
#     R: gas constant [J/mol /K] 
#     T: annual mean 2m air temperature in [K]
#     Tref: 284.15K


# Use at least multiple year means of atmospheric data for smooth changes(!)

# Script contains and runs:
#
# readlithology()     Reads input data for the weathering model   ! called for Palmod
# readconstants()     Defines model parameters                    ! called for Palmod
# weatheringmain()    Weathering equations                        ! called for Palmod
# land2ocean_trans()  Parametrizes global land-to-ocean transformations  ! turned off for Palmod
# riv_spatial_distr() River input spatial distribution            ! turned off for Palmod
# summary()           Prints summary in konsole                   ! called for Palmod
# writeresults2nc()   Write nc file output                        ! called for Palmod
# main()              Calls everything                            ! called for Palmod
#
# option to change echam input files, call as following:
# python weatherin.py $ECHAM_INPUTFILE $ECHAM_GRIDINFO


import netCDF4 as nc4
from netCDF4 import Dataset
import numpy as np
from array import array
import math
import sys

def main():

    # echam input options
    if   len(sys.argv) == 1:
        echam_filename = "echam_inputs/echam_test_data.nc"
        echam_area_file = "echam_inputs/echam_area.nc"
    elif len(sys.argv) == 2:      
        echam_filename  =  str(sys.argv[1])
        echam_area_file = "echam_inputs/echam_area.nc"
    elif len(sys.argv) == 3:
        echam_filename  =  str(sys.argv[1])
        echam_area_file =  str(sys.argv[2])

    print( 'Run weathering model with the following echam input file:', echam_filename)
         
    # run program
    # read input data files
    lithologydata,soil_shdata,runoff_yearly,drainage_yearly, \
    t2m_yearly,glac_covr,area = readlithology(echam_filename)

    # Define weathering model parameters
    Ea,bS_emp_P,bC_emp_P,bS_emp_alk,bC_emp_alk = constantssoil()

    # Run weathering model from weatheringmain() function
    F_temp,\
    p_release_lith,\
    alk_release_lith,\
    co2_release_lith,\
    p_release_global,\
    alk_release_global, \
    co2_release_global \
                     = weatheringmain(lithologydata,\
                     soil_shdata,\
                     runoff_yearly,\
                     drainage_yearly,\
                     t2m_yearly,\
                     glac_covr,\
                     area,\
                     Ea,\
                     bS_emp_P,\
                     bC_emp_P,\
                     bS_emp_alk,\
                     bC_emp_alk)

    # Print results summary
    summary(t2m_yearly,drainage_yearly,runoff_yearly,F_temp, \
    p_release_global,alk_release_global,co2_release_global)
    # Write nc
    writeresults2nc(lithologydata,runoff_yearly,p_release_lith,alk_release_lith,co2_release_lith)

def readlithology(echam_filename):

# Define global variables:

# Read lithology file: lithology and soil shield (SS) from Hartmann et al. (2009) and Hartmann et al. (2014)
    lithology_nc_file = 'wthr_model_inputs/lithologyres_lith.nc' #save this file somewhere public and use this link for dynamic runs
    fh = Dataset(lithology_nc_file, mode='r')
    lithology = fh.variables['lith_frac'][:,:,:]
    lithology[lithology < 0] = np.nan

    soil_shield_nc_file = 'wthr_model_inputs/lithologyres_soil_shield.nc'
    fh = Dataset(soil_shield_nc_file, mode='r')
    soil_sh = fh.variables['soil_shield'][:,:]
    soil_sh[soil_sh < 0] = np.nan
    fh.close()
# flip y axis from South (i=0) to North (i=95) for consistency with ECHAM output
    lithologydata = lithology.data[:,::-1,:]
    soil_shdata=soil_sh[:,:]


# Read computed 100-year climatologies for runoff, drainage and temperature from ECHAM CMIP5 output 
#This needs to be changed with the relative path of PALMOD echam output (average over t10 years in shell script before)
    echam_nc_file  = echam_filename
    fh = Dataset(echam_nc_file, mode='r')
    runoff_s = fh.variables['var160'][:,:,:]
    drainage_s = fh.variables['var161'][:,:,:]
    t2m = fh.variables['var167'][:,:,:]
    t2m[t2m > 9E36] = np.nan
    runoff_s[runoff_s > 9E36] = np.nan
    drainage_s[drainage_s > 9E36] = np.nan
    fh.close()
    
# Glacier ice coverage
# The input stream for this would need to be added (!) 

    glac_covr = np.zeros((96,192))  

# T63 grid information 

    area_nc_file= 'wthr_model_inputs/test.nc'
    fh = Dataset(area_nc_file, mode='r')
    area = fh.variables['cell_area'][:,:]
    area[area>10E11] = np.nan
    fh.close()

#Compute annual means
#Runoff is in mm/s for every model output timestep so first needs to be calculated in mm/a then averaged over the time series 
    #scaling constant to account for runoff underestimation in ECHAM (in Lacroix et al., 2020, f_sclr = 1.59)
    f_sclr = 1.10   
    # seconds in year
    s2yr            = 3.1536e+07
    runoff_yearly   = f_sclr * runoff_s.mean(axis=0) * s2yr #calculate to [mm/yr] and average over all time steps
    drainage_yearly = f_sclr * drainage_s.mean(axis=0) * s2yr # calculate to [mm/yr] and average over all time steps
    t2m_yearly      = t2m.mean(axis=0) # [K] average over all time steps

    return lithologydata,soil_shdata,runoff_yearly,drainage_yearly,t2m_yearly,glac_covr,area


def constantssoil():
  
# Define global constants for activation temperatures of lithologies Ea, and weathering parameters (see Hartmann et al. 2014)
# Lithologies are in the same array order as in the lithology nc file (from Hartmann et al., 2019)
#!    ! 1  - evaporites
#!    ! 2  - ice & glaciers
#!    ! 3  - metamorphics
#!    ! 4  - no data
#!    ! 5  - acid plutonic rocks
#!    ! 6  - basic plutonic rocks
#!    ! 7  - intermediat plutonic rocks
#!    ! 8  - pyroclastics
#!    ! 9  - carbonate sedimentary rocks 
#!    ! 10 - mixed sedimentary rocks 
#!    ! 11 - siliciclastic sedimentary rocks
#!    ! 12 - unconsolidated sediments
#!    ! 13 - acid volcanic rocks
#!    ! 14 - basic volcanic rocks
#!    ! 15 -0 intermediate volcanic rocks
#!    ! 16 - water bodies

    Ea         = np.array([0.0, 0.0, 6.0E4, 0, 6.0E4, 5.0E04, 6.0E4, 4.6E4, 6.0E4,\
                           6.0E4, 6.0E4, 6.0E4, 6.0E4, 5.0E4, 5.0E4, 0.0])
# Empirical factor to relate P release by silicate weathering to runoff
    CWSi_P     = np.array([0.0, 0.0, 0.01, 1.25E-5, 0.019307, 0.04054, 0.19307, \
                           0.076876, 0.0, 0.21806, 0.019699, 0.021333, 0.020762, \
                           0.04054, 0.04054, 0.0])
    frac_P     = np.array([0.000624, 0.0, 0.001048, 0.0010, 0.000915, 0.003504, \
                           0.002, 0.002, 0.001143, 0.000814, 0.000766, 0.000814, \
                           0.000382, 0.002102, 0.001, 0.0])
# Empirical factor to relate P release by carbonate weathering to runoff
    CWCa_P     = np.array([0.151243, 0.0, 0.020118, 0.0, 0.007603, 0.007603, 0.0, \
                           0.0, 0.151243, 0.035915, 0.0081, 0.0, 0.0, 0.0, 0.0, 0.0])

    SiO2_abs   = np.array([15.87, 0.0, 0.0, 0.0, 67.78, 46.39, 57.09, 0.9, 7.52, \
                           52.3, 52.3, 60.32, 72.7, 51.26, 61.98, 0.0])
    Cat_SiO2   =  np.array([41.98, 0.0, 0.0, 0.0, 76.28, 64.77, 71.00, 1, 42, \
                            64.36, 68.41, 64.36, 79.94, 64.35, 72.145, 0.0])    
# b value for alkalinity
    b_lith_alk = np.array([0.0, 0.0, 0.0076, 0.0, 0.005095, 0.007015,0.007015,0.0061, \
                           0.03804,0.012481,0.005341,0.003364,0.002455,0.007015,0.007015,0.0])
    CWSi_alk   = np.array([0.0, 0.0, 0.25, 0.0, 0.58, 1.0, 0.58, 1.0, 0.0, 2.4, 0.64, 1.0, \
                           1.0,1.0,1.0,0.0])
    CWCa_alk   =  np.array([0.0,0.0,0.75,0.0,0.42,0.0,0.42,0.0,1.0,0.76,0.36,0.0,0,0,0,0]) 

# Empirical factor to calculate weathered P and Alk from amount of rock weathering for each lithology
    bS_emp_P = CWSi_P * frac_P
    bC_emp_P = CWCa_P * frac_P
    bS_emp_alk = b_lith_alk * CWSi_alk
    bC_emp_alk = b_lith_alk * CWCa_alk 
#   SiO2_rel = SiO2_abs/Cat_SiO2;  #  to calculate DSi, does not work, do not use

    return Ea,bS_emp_P,bC_emp_P,bS_emp_alk,bC_emp_alk
   
def weatheringmain(lithologydata,soil_shdata,runoff_yearly,drainage_yearly,t2m_yearly,glac_covr,area,Ea,bS_emp_P,bC_emp_P,bS_emp_alk,bC_emp_alk):
    
    #Initialize & set dimensions
    kdim0          = np.shape(lithologydata)[0]
    kdim1          = np.shape(runoff_yearly)[0]
    kdim2          = np.shape(runoff_yearly)[1]
    F_temp         = np.zeros((kdim0,kdim1,kdim2))
    p_release_lith = np.zeros((kdim0,kdim1,kdim2))
    alk_release_lith = np.zeros((kdim0,kdim1,kdim2))
    co2_release_lith = np.zeros((kdim0,kdim1,kdim2))

    #local constants
    R    = 9.81        # Gas constant [J/mol/K]
    Tref = 284.15      # Hartmann et al. (2014) [K]
    km2m = 10E6        # convert km2 to m2
    t2g  = 10E6        # convert tons to g      
    # loop over lithologies and 2-d field to compute weathering fluxes
    for dsg0 in range(0,kdim0-1):
        for dsg1 in range(0,kdim1-1):
            for dsg2 in range(0,kdim2-1):
                # Temperature dependence function of weatherin
                if Ea[dsg0] != 0.0:
                    F_temp[dsg0,dsg1,dsg2] = \
                    math.exp(-Ea[dsg0]/R *(1.0/t2m_yearly[dsg1,dsg2]-1.0/Tref))
                else:
                    F_temp[dsg0,dsg1,dsg2] = 0.0 

                #P Weathering release in [mm/km2/yr]
                p_release_lith[dsg0,dsg1,dsg2] = \
                (F_temp[dsg0,dsg1,dsg2]*bS_emp_P[dsg0] + 
                bC_emp_P[dsg0]) * \
                lithologydata[dsg0,dsg1,dsg2] * \
                (runoff_yearly[dsg1,dsg2] + drainage_yearly[dsg1,dsg2]) * \
                soil_shdata[dsg1,dsg2] #* \
                #glac_covr[dsg1,dsg2]

                #Alk Weathering release
                alk_release_lith[dsg0,dsg1,dsg2] = \
                (F_temp[dsg0,dsg1,dsg2]*bS_emp_alk[dsg0] +
                bC_emp_alk[dsg0]) * \
                lithologydata[dsg0,dsg1,dsg2] * \
                (runoff_yearly[dsg1,dsg2] + drainage_yearly[dsg1,dsg2]) * \
                soil_shdata[dsg1,dsg2] #* \
                #glac_covr[dsg1,dsg2]

                #CO2 Weathering release
                co2_release_lith[dsg0,dsg1,dsg2] = \
                (F_temp[dsg0,dsg1,dsg2]*bS_emp_alk[dsg0] +
                1/2 * bC_emp_alk[dsg0]) * \
                lithologydata[dsg0,dsg1,dsg2] * \
                (runoff_yearly[dsg1,dsg2] + drainage_yearly[dsg1,dsg2]) * \
                soil_shdata[dsg1,dsg2] #* \


    p_release_global   =  np.nansum(p_release_lith * area / km2m) * t2g  # multiply with area convert from km-2 to m2, and from t yr-1 for output in g P/yr-1
    alk_release_global =  np.nansum(alk_release_lith * area / km2m) * t2g  # multiply with area, convert from km-2 to m2 and from t yr-1 for output in g P/yr-1
    co2_release_global =  np.nansum(co2_release_lith * area / km2m) * t2g
    return F_temp,p_release_lith,alk_release_lith,co2_release_lith,p_release_global,alk_release_global,co2_release_global
 
def summary(t2m_yearly,drainage_yearly,runoff_yearly,F_temp,p_release_global,alk_release_global,co2_release_global):
    print("land t2m mean [K] = ",np.nanmean(t2m_yearly))
    print("land surface runoff mean [mm] = ",np.nanmean(runoff_yearly+drainage_yearly))
    print("weathering F_temp mean = ", np.nanmean(F_temp))
    print("Global P weathering release averaged = ", p_release_global, "[g/yr]")
    print("Global Alk weathering release averaged = ", alk_release_global, "[g/yr]")
    print("Global CO2 weathering release averaged = ", co2_release_global, "[g/yr]")
#    print("Global DIP river load  = ", DIP_load_global, "[g/yr]")
#    print("Global DIN river load  = ", DIN_load_global, "[g/yr]")
#    print("Global DSi river load  = ", DSi_load_global, "[g/yr]")
#    print("Global DFe river load  = ", DFe_load_global, "[g/yr]")
#    print("Global Alk river load  = ", Alk_load_global, "[g/yr]")
#    print("Global DIC river load  = ", DIC_load_global, "[g/yr]") 
#    print("Global tDOM river load  = ", tDOM_load_global, "[g/yr]")
#    print("Global POM river load  = ", tPOM_load_global, "[g/yr]")  

def writeresults2nc(lithologydata,runoff_yearly,p_release_lith,alk_release_lith,co2_release_lith):

    kdim0 = np.shape(lithologydata)[0]
    kdim1 = np.shape(runoff_yearly)[0]
    kdim2 = np.shape(runoff_yearly)[1]

    #kdim_oc0 = np.shape(riv_runoff_yearly)[0]
    #kdim_oc1 = np.shape(riv_runoff_yearly)[1]

    f = nc4.Dataset('test_runoff.nc','w', format='NETCDF4') #'w' stands for write
    tempgrp = f.createGroup('Temp_data')
    tempgrp.createDimension('lith', kdim0)
    tempgrp.createDimension('y', kdim1)
    tempgrp.createDimension('x', kdim2)
    lithologyaxis = tempgrp.createVariable('lith', 'i4','lith')
    xaxis=tempgrp.createVariable('x', 'i4','x')
    yaxis=tempgrp.createVariable('y', 'i4','y')
    temp = tempgrp.createVariable('runoff_yearly', 'f4', ('y', 'x'))

    lithdata = np.arange(kdim0)
    xdata = np.arange(kdim2)
    ydata = np.arange(kdim1)

    xaxis[:] =  xdata
    yaxis[:] =  ydata
    temp[:,:] = runoff_yearly
    f.close

    f = nc4.Dataset('test_prelease.nc','w', format='NETCDF4') #'w' stands for write
    f.createDimension('y', kdim1)
    f.createDimension('x', kdim2)
    yaxis=f.createVariable('y', 'i4','y')
    xaxis=f.createVariable('x', 'i4','x')    
    temp = f.createVariable('P_release', 'f4', ('y', 'x')) 
  
    temp[:,:] = np.nansum(p_release_lith,axis=0)  
    f.close 

    f = nc4.Dataset('test_alkrelease.nc','w', format='NETCDF4') #'w' stands for write
    f.createDimension('y', kdim1)
    f.createDimension('x', kdim2)
    yaxis = f.createVariable('y', 'i4','y')
    xaxis = f.createVariable('x', 'i4','x')
    temp = f.createVariable('weathering_hco3_release', 'f4', ('y', 'x'))

    xaxis[:] =  xdata
    yaxis[:] =  ydata
    temp[:,:] = np.nansum(alk_release_lith,axis=0)
    f.close

    f = nc4.Dataset('test_co2drawdown.nc','w', format='NETCDF4') #'w' stands for write
    f.createDimension('y', kdim1)
    f.createDimension('x', kdim2)
    yaxis = f.createVariable('y', 'i4','y')
    xaxis = f.createVariable('x', 'i4','x')
    temp = f.createVariable('weathering_co2_drawdown', 'f4', ('y', 'x'))

    xaxis[:] =  xdata
    yaxis[:] =  ydata
    temp[:,:] = np.nansum(co2_release_lith,axis=0)
    f.close

main()




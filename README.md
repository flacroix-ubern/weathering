# weathering
Script to compute spatial and global mean uptake of CO2 and release of HCO3- and P through weathering

The script computes the global p weathering release derived using a spatiotemporal model as described in Hartmann et al. 2014
and Lacroix et al. (2020).

  The P weathering release (Pr) is calculated:                
     Pr = b * q *f_t * s * glac_covr
     where
     q: annual runoff [mm/a]
     s: factor to correct for soil shielding
     b: empirical lithology dependent factors 
     glac_covr: Glacier coverage [-]
     temperature response term f_t is given by:
     f_t = exp(-Ea/R*(1/T-1/Tref))
     with
     Ea:  apparent activation energy (depends on lithology)
     R: gas constant [J/mol /K] 
     T: annual mean 2m air temperature in [K]
     Tref: 284.15K


 Use at least multiple year means of atmospheric data for smooth changes(!)

 Script contains and runs:

 readlithology()     Reads input data for the weathering model   ! called for Palmod
 readconstants()     Defines model parameters                    ! called for Palmod
 weatheringmain()    Weathering equations                        ! called for Palmod
 land2ocean_trans()  Parametrizes global land-to-ocean transformations  ! turned off for Palmod
 riv_spatial_distr() River input spatial distribution            ! turned off for Palmod
 summary()           Prints summary in konsole                   ! called for Palmod
 writeresults2nc()   Write nc file output                        ! called for Palmod
 main()              Calls everything                            ! called for Palmod

 option to change echam input files, call as following:
 python weathering.py $ECHAM_INPUTFILE $ECHAM_GRIDINFO
 
 The weathering script can also be called through the shell script weathering.sh:#
  sh weathering.sh $ECHAM_INPUTFILE $ECHAM_GRIDINFO

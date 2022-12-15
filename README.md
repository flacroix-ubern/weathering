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

 call as following:
 
 python weathering.py $ECHAM_INPUTFILE $ECHAM_GRIDINFO
 
 The weathering script can also be called through the shell script weathering.sh:
 
  sh weathering.sh $ECHAM_INPUTFILE $ECHAM_GRIDINFO

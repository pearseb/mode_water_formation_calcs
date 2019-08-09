#!/bin/bash


# This script calls the procedures required to calculate: 
#	(1) the formation rates of SAMW and AAIW
#	(2) the contribution of nitrification to SAMW and AAIW NO3 content
#	(3) the contribution of nitrification to the Si* of these water masses



# (1.1) regrid physical properties (T,S,MLD) to regular grid
ferret -script regrid_phys_properties_2x1.jnl

# (1.2) save mixedlayer and density files 
ncks -O -v somxl010 ../RUNDIR4/dyna_grid_T.nc mixedlayer.nc
ncrename -v somxl010,mld mixedlayer.nc
ferret -script create_rho.jnl

# (1.3) find isopycnal bounds of AAIW and SAMW
#	transfer the following files to a directory where python is working:
#		- salinity_regular2x1.nc
#		- potentialdensity_regular2x1.nc
#		- salinity_regular2x1_woa18.nc
#		- potentialdensity_regular2x1_woa18.nc
#		- mixedlayer.nc
#		- rho.nc
#	execute the following scripts in order:
#		- modewater_isopycnal_bounds.py
#		- properties_at_modewater_outcrops.py
#	and transfer the following files back:
#		- transports_at_modewater_outcrops_3D.nc  (lon x lat x time)
#		- properties_at_modewater_outcrops_3D.nc  (lon x lat x time)
#		- transports_at_modewater_outcrops_4D.nc  (lon x lat x dep x time)
#		- properties_at_modewater_outcrops_4D.nc  (lon x lat x dep x time)

# (1.4) regrid the 3D transports and properties onto regular grid 
ferret -script regrid_transports_2x1.jnl
ferret -script regrid_bgc_properties_2x1.jnl

# (1.5) execute the "calculate_modewater_formation.jnl" script
#	saves:
#		- ./results/modewater_formation_3Dfield.nc
#		- ./results/modewater_formation_timeseries.nc
#		- ./results/modewater_formation_3Dfield_winteraverage.nc
#	and lists the subduction and obduction rates in Sv for AAIW and SAMW	
ferret -script calculate_modewater_formation.jnl

# (2) execute the "calculate_nitrif_contrib_no3trans.jnl" script
#	saves:
#		- ./results/nitrate_transformations_timeseries.nc (mols/s)
#	and lists:
#		- mean NO3 concentration at AAIW and SAMW sites of outcropping (mmol/m3)
#		- integrated no3 transports at aaiw and samw sites (mmol/s) 
#		- mean nitrification of no3 at the sites of aaiw and samw subduction/obduction (mmol/s)
#		- integrated nitrification of no3 at the sites of aaiw and samw subduction/obduction (mmol/s)
#		- percentage of nitrified NO3 SUBDUCTED into AAIW and SAMW as a percentage of total NO3 (WINTER)
#               - percentage of nitrified NO3 OBDUCTED into AAIW and SAMW as a percentage of total NO3 (WINTER)
#		- percentage of nitrified NO3 SUBDUCTED into AAIW and SAMW as a percentage of total NO3 (SUMMER)
#               - percentage of nitrified NO3 OBDUCTED into AAIW and SAMW as a percentage of total NO3 (SUMMER)
#		- percentage of nitrified NO3 SUBDUCTED into AAIW and SAMW as a percentage of total NO3 (ANNUAL)
#               - percentage of nitrified NO3 OBDUCTED into AAIW and SAMW as a percentage of total NO3 (ANNUAL)
ferret -script calculate_nitrif_contrib_no3trans.jnl

# (3) execute the "calculate_sistar_contrib.jnl" script



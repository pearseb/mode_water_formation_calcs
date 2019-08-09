# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 15:50:23 2017


Purpose
-------
    Produce 2D fields of horizontal and bottom vertical velocities at 
    grid cells that represent the maximum mixed layer base


Inputs
------
    The depth of the maximum annual mixed layer
    The horizontal and vertical velocities through the water column at the 
        time of the maximum annual mixed layer
    The area of each grid box
    3D annual average density 


Outputs
-------
    2D (lat x lon) arrays of u, v, and w velocities at the mixed layer depth
    2D (lat x lon) array of total transports through base of mixed layer depth
    Total integrated transports into isopycnals corresponding with AAIW


Steps:
    1. Find the grid boxes (k-level) where the maximum mixing occurs
    2. Save the zonal, meridional and vertical transports at points and times
       when maximum mixing occurs
    3. Save these arrays to the existing netcdf files


@author: pearseb
"""


#%% imports

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


#%% get data

os.chdir("C://Users//pearseb//Dropbox//PostDoc//model output//calculate AAIW formation rates")

data = nc.Dataset('mixedlayer.nc','r')
mld = np.ma.getdata(data.variables['mld'][...])
data = nc.Dataset('rho.nc','r')
rho = data.variables['RHO'][...]
data = nc.Dataset('dyna_grid_U.nc','r')
u = data.variables['uocetr_eff'][...]
data = nc.Dataset('dyna_grid_V.nc','r')
v = data.variables['vocetr_eff'][...]
data = nc.Dataset('dyna_grid_W.nc','r')
w = data.variables['wocetr_eff'][...]

time = data.variables['time_centered'][:]/(86400*365)+1899
time_bnds = data.variables['time_centered_bounds'][:]/(86400*365)+1899
dep = data.variables['depthw'][:]
dep_bnds = data.variables['depthw_bounds'][:]
lats = data.variables['nav_lat'][:]
lons = data.variables['nav_lon'][:]

inputfile = 'TEST_5d_ptrc_Y1395.nc'
data = nc.Dataset(inputfile,'r')
nitri = data.variables['NITR'][...]
no3 = data.variables['NO3'][...]
si = data.variables['Si'][...]

data = nc.Dataset('ORCA2.0.full_grid.nc','r')
vol = data.variables['volume'][...]

data.close()



#%% get the isopycnal bounds

iso_bnds = np.load('modewater_isopycnal_bounds.npy')
print(iso_bnds)


#%% Find the depth levels of maximum mixing

klevel = np.zeros(np.shape(mld))
for i in np.arange(len(mld[0,0,:])):
    print(i)
    for j in np.arange(len(mld[0,:,0])):
        for t in np.arange(len(mld[:,0,0])):
            for k,z in enumerate(dep):
                if mld[t,j,i] == z:
                    klevel[t,j,i] = k 


#%% select only transports (and nitrification) at mixed layer depths where isopycnals of AAIW outcrop

u_mix_aaiw = np.zeros(np.shape(mld))
v_mix_aaiw = np.zeros(np.shape(mld))
w_mix_aaiw = np.zeros(np.shape(mld))

u_mix_samw = np.zeros(np.shape(mld))
v_mix_samw = np.zeros(np.shape(mld))
w_mix_samw = np.zeros(np.shape(mld))

nit_mix_aaiw = np.zeros(np.shape(mld))
nit_mix_samw = np.zeros(np.shape(mld))

no3_mix_aaiw = np.zeros(np.shape(mld))
no3_mix_samw = np.zeros(np.shape(mld))

si_mix_aaiw = np.zeros(np.shape(mld))
si_mix_samw = np.zeros(np.shape(mld))

vol_mix_aaiw = np.zeros(np.shape(mld))
vol_mix_samw = np.zeros(np.shape(mld))


for i in np.arange(len(u_mix_aaiw[0,0,:])):
    print(i)
    for j in np.arange(len(u_mix_aaiw[0,:,0])):
        for t in np.arange(len(u_mix_aaiw[:,0,0])):
            for k,z in enumerate(dep):
                #kp1 = k+1 # select the depth level below the maximum mixed layer depth
                if klevel[t,j,i] == k:
                    if rho[t,k,j,i] > iso_bnds[0] and rho[t,k,j,i] < iso_bnds[1]:
                        u_mix_aaiw[t,j,i] = u[t,k,j,i]
                        v_mix_aaiw[t,j,i] = v[t,k,j,i]
                        w_mix_aaiw[t,j,i] = w[t,k,j,i]
                        nit_mix_aaiw[t,j,i] = nitri[t,k,j,i]
                        no3_mix_aaiw[t,j,i] = no3[t,k,j,i]
                        si_mix_aaiw[t,j,i] = si[t,k,j,i]
                        vol_mix_aaiw[t,j,i] = vol[k,j,i]
                    else:
                        nit_mix_aaiw[t,j,i] = np.nan
                        no3_mix_aaiw[t,j,i] = np.nan
                        si_mix_aaiw[t,j,i] = np.nan
                        vol_mix_aaiw[t,j,i] = np.nan
                    
                    if rho[t,k,j,i] > iso_bnds[0]-0.5 and rho[t,k,j,i] < iso_bnds[0]:
                        u_mix_samw[t,j,i] = u[t,k,j,i]
                        v_mix_samw[t,j,i] = v[t,k,j,i]
                        w_mix_samw[t,j,i] = w[t,k,j,i]
                        nit_mix_samw[t,j,i] = nitri[t,k,j,i]
                        no3_mix_samw[t,j,i] = no3[t,k,j,i]
                        si_mix_samw[t,j,i] = si[t,k,j,i]
                        vol_mix_samw[t,j,i] = vol[k,j,i]
                    else:
                        nit_mix_samw[t,j,i] = np.nan
                        no3_mix_samw[t,j,i] = np.nan
                        si_mix_samw[t,j,i] = np.nan
                        vol_mix_samw[t,j,i] = np.nan


nit_mix_aaiw = np.ma.masked_where(np.isnan(nit_mix_aaiw), nit_mix_aaiw)
no3_mix_aaiw = np.ma.masked_where(np.isnan(no3_mix_aaiw), no3_mix_aaiw)
si_mix_aaiw = np.ma.masked_where(np.isnan(si_mix_aaiw), si_mix_aaiw)
vol_mix_aaiw = np.ma.masked_where(np.isnan(vol_mix_aaiw), vol_mix_aaiw)

nit_mix_samw = np.ma.masked_where(np.isnan(nit_mix_samw), nit_mix_samw)
no3_mix_samw = np.ma.masked_where(np.isnan(no3_mix_samw), no3_mix_samw)
si_mix_samw = np.ma.masked_where(np.isnan(si_mix_samw), si_mix_samw)
vol_mix_samw = np.ma.masked_where(np.isnan(vol_mix_samw), vol_mix_samw)


### NO3 subducted = average NO3 (mol/m3) * subduction rate (m3/s) = mol/s
### NO3 obducted = average NO3 (mol/m3) * obduction rate (m3/s) = mol/s
### NO3 nitrified = nitrification rate (mol/m3/s) * volume (m3) = mol/s


#%% check arrays

plt.figure()
plt.contourf(np.mean(u_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(v_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(w_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(nit_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(vol_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(no3_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(si_mix_aaiw,axis=0),corner_mask=False)
plt.colorbar()


#%% save as netcdfs

import shutil
os.remove('transports_at_modewater_outcrops_3D.nc')
shutil.copyfile('mixedlayer.nc', 'transports_at_modewater_outcrops_3D.nc')


data = nc.Dataset('transports_at_modewater_outcrops_3D.nc','a')

Uaaiw = data.createVariable('Utrans_AAIW','f8',('time_counter','y','x'))
Vaaiw = data.createVariable('Vtrans_AAIW','f8',('time_counter','y','x'))
Waaiw = data.createVariable('Wtrans_AAIW','f8',('time_counter','y','x'))
Usamw = data.createVariable('Utrans_SAMW','f8',('time_counter','y','x'))
Vsamw = data.createVariable('Vtrans_SAMW','f8',('time_counter','y','x'))
Wsamw = data.createVariable('Wtrans_SAMW','f8',('time_counter','y','x'))


Uaaiw.units = 'm3/s east'
Vaaiw.units = 'm3/s north'
Waaiw.units = 'm3/s up'
Usamw.units = 'm3/s east'
Vsamw.units = 'm3/s north'
Wsamw.units = 'm3/s up'

data.variables['Utrans_AAIW'][...] = u_mix_aaiw
data.variables['Vtrans_AAIW'][...] = v_mix_aaiw
data.variables['Wtrans_AAIW'][...] = w_mix_aaiw
data.variables['Utrans_SAMW'][...] = u_mix_samw
data.variables['Vtrans_SAMW'][...] = v_mix_samw
data.variables['Wtrans_SAMW'][...] = w_mix_samw


data.close()


#%%

os.remove('properties_at_modewater_outcrops_3D.nc')
shutil.copyfile(inputfile, 'properties_at_modewater_outcrops_3D.nc')

data = nc.Dataset('properties_at_modewater_outcrops_3D.nc','a')

Nitaaiw3d = data.createVariable('Nitrif_AAIW','f8',('time_counter','y','x'))
Nitsamw3d = data.createVariable('Nitrif_SAMW','f8',('time_counter','y','x'))
NO3aaiw3d = data.createVariable('Nitrat_AAIW','f8',('time_counter','y','x'))
NO3samw3d = data.createVariable('Nitrat_SAMW','f8',('time_counter','y','x'))
Silaaiw3d = data.createVariable('Silica_AAIW','f8',('time_counter','y','x'))
Silsamw3d = data.createVariable('Silica_SAMW','f8',('time_counter','y','x'))
Volaaiw3d = data.createVariable('Volume_AAIW','f8',('time_counter','y','x'))
Volsamw3d = data.createVariable('Volume_SAMW','f8',('time_counter','y','x'))

Nitaaiw3d.units = 'mol/m3/s'
Nitsamw3d.units = 'mol/m3/s'
NO3aaiw3d.units = 'mmol/m3'
NO3samw3d.units = 'mmol/m3'
Silaaiw3d.units = 'mmol/m3'
Silsamw3d.units = 'mmol/m3'
Volaaiw3d.units = 'm3'
Volsamw3d.units = 'm3'

data.variables['Nitrif_AAIW'][...] = nit_mix_aaiw
data.variables['Nitrif_SAMW'][...] = nit_mix_samw
data.variables['Nitrat_AAIW'][...] = no3_mix_aaiw
data.variables['Nitrat_SAMW'][...] = no3_mix_samw
data.variables['Silica_AAIW'][...] = si_mix_aaiw
data.variables['Silica_SAMW'][...] = si_mix_samw
data.variables['Volume_AAIW'][...] = vol_mix_aaiw
data.variables['Volume_SAMW'][...] = vol_mix_samw

data.close()


#%% select only transports (and nitrification) above the mixed layer depths where isopycnals of AAIW outcrop

u_mix_aaiw_wml = np.zeros(np.shape(u))
v_mix_aaiw_wml = np.zeros(np.shape(u))
w_mix_aaiw_wml = np.zeros(np.shape(u))

u_mix_samw_wml = np.zeros(np.shape(u))
v_mix_samw_wml = np.zeros(np.shape(u))
w_mix_samw_wml = np.zeros(np.shape(u))

nit_mix_aaiw_wml = np.zeros(np.shape(u))
nit_mix_samw_wml = np.zeros(np.shape(u))

no3_mix_aaiw_wml = np.zeros(np.shape(u))
no3_mix_samw_wml = np.zeros(np.shape(u))

si_mix_aaiw_wml = np.zeros(np.shape(u))
si_mix_samw_wml = np.zeros(np.shape(u))

vol_mix_aaiw_wml = np.zeros(np.shape(u))
vol_mix_samw_wml = np.zeros(np.shape(u))


for i in np.arange(len(u_mix_aaiw_wml[0,0,0,:])):
    print(i)
    for j in np.arange(len(u_mix_aaiw_wml[0,0,:,0])):
        for t in np.arange(len(u_mix_aaiw_wml[:,0,0,0])):
            for k,z in enumerate(dep):
                #kp1 = k+1 # select the depth level below the maximum mixed layer depth
                if k < klevel[t,j,i] and np.ma.is_masked(nit_mix_aaiw[t,j,i]) == False:
                    if rho[t,k,j,i] > iso_bnds[0] and rho[t,k,j,i] < iso_bnds[1]:
                        u_mix_aaiw_wml[t,k,j,i] = u[t,k,j,i]
                        v_mix_aaiw_wml[t,k,j,i] = v[t,k,j,i]
                        w_mix_aaiw_wml[t,k,j,i] = w[t,k,j,i]
                        nit_mix_aaiw_wml[t,k,j,i] = nitri[t,k,j,i]
                        no3_mix_aaiw_wml[t,k,j,i] = no3[t,k,j,i]
                        si_mix_aaiw_wml[t,k,j,i] = si[t,k,j,i]
                        vol_mix_aaiw_wml[t,k,j,i] = vol[k,j,i]
                    else:
                        nit_mix_aaiw_wml[t,k,j,i] = np.nan
                        no3_mix_aaiw_wml[t,k,j,i] = np.nan
                        si_mix_aaiw_wml[t,k,j,i] = np.nan
                        vol_mix_aaiw_wml[t,k,j,i] = np.nan
                else:
                    nit_mix_aaiw_wml[t,k,j,i] = np.nan
                    no3_mix_aaiw_wml[t,k,j,i] = np.nan
                    si_mix_aaiw_wml[t,k,j,i] = np.nan
                    vol_mix_aaiw_wml[t,k,j,i] = np.nan
                
                if k < klevel[t,j,i] and np.ma.is_masked(nit_mix_samw[t,j,i]) == False:
                    if rho[t,k,j,i] > iso_bnds[0]-0.5 and rho[t,k,j,i] < iso_bnds[0]:
                        u_mix_samw_wml[t,k,j,i] = u[t,k,j,i]
                        v_mix_samw_wml[t,k,j,i] = v[t,k,j,i]
                        w_mix_samw_wml[t,k,j,i] = w[t,k,j,i]
                        nit_mix_samw_wml[t,k,j,i] = nitri[t,k,j,i]
                        no3_mix_samw_wml[t,k,j,i] = no3[t,k,j,i]
                        si_mix_samw_wml[t,k,j,i] = si[t,k,j,i]
                        vol_mix_samw_wml[t,k,j,i] = vol[k,j,i]
                    else:
                        nit_mix_samw_wml[t,k,j,i] = np.nan
                        no3_mix_samw_wml[t,k,j,i] = np.nan
                        si_mix_samw_wml[t,k,j,i] = np.nan
                        vol_mix_samw_wml[t,k,j,i] = np.nan
                else:
                    nit_mix_samw_wml[t,k,j,i] = np.nan
                    no3_mix_samw_wml[t,k,j,i] = np.nan
                    si_mix_samw_wml[t,k,j,i] = np.nan
                    vol_mix_samw_wml[t,k,j,i] = np.nan



nit_mix_aaiw_wml = np.ma.masked_where(np.isnan(nit_mix_aaiw_wml), nit_mix_aaiw_wml)
no3_mix_aaiw_wml = np.ma.masked_where(np.isnan(no3_mix_aaiw_wml), no3_mix_aaiw_wml)
si_mix_aaiw_wml = np.ma.masked_where(np.isnan(si_mix_aaiw_wml), si_mix_aaiw_wml)
vol_mix_aaiw_wml = np.ma.masked_where(np.isnan(vol_mix_aaiw_wml), vol_mix_aaiw_wml)

nit_mix_samw_wml = np.ma.masked_where(np.isnan(nit_mix_samw_wml), nit_mix_samw_wml)
no3_mix_samw_wml = np.ma.masked_where(np.isnan(no3_mix_samw_wml), no3_mix_samw_wml)
si_mix_samw_wml = np.ma.masked_where(np.isnan(si_mix_samw_wml), si_mix_samw_wml)
vol_mix_samw_wml = np.ma.masked_where(np.isnan(vol_mix_samw_wml), vol_mix_samw_wml)


### NO3 subducted = average NO3 (mol/m3) * subduction rate (m3/s) = mol/s
### NO3 obducted = average NO3 (mol/m3) * obduction rate (m3/s) = mol/s
### NO3 nitrified = nitrification rate (mol/m3/s) * volume (m3) = mol/s


#%%


plt.figure()
plt.contourf(u_mix_aaiw_wml[50,16,:,:],corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(v_mix_aaiw_wml[50,10,:,:],corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(np.mean(w_mix_aaiw_wml[:,0,:,:],axis=0),corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(nit_mix_aaiw_wml[0,0,:,:],corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(vol_mix_aaiw_wml[50,0,:,:],corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(no3_mix_aaiw_wml[50,1,:,:],corner_mask=False)
plt.colorbar()

plt.figure()
plt.contourf(si_mix_aaiw_wml[50,1,:,:],corner_mask=False)
plt.colorbar()


#%% save as netcdf

import shutil
os.remove('transports_at_modewater_outcrops_4D.nc')
shutil.copyfile('mixedlayer.nc', 'transports_at_modewater_outcrops_4D.nc')


data = nc.Dataset('transports_at_modewater_outcrops_4D.nc','a')

Uaaiw = data.createVariable('Utrans_AAIW','f8',('time_counter','deptht','y','x'))
Vaaiw = data.createVariable('Vtrans_AAIW','f8',('time_counter','deptht','y','x'))
Waaiw = data.createVariable('Wtrans_AAIW','f8',('time_counter','deptht','y','x'))
Usamw = data.createVariable('Utrans_SAMW','f8',('time_counter','deptht','y','x'))
Vsamw = data.createVariable('Vtrans_SAMW','f8',('time_counter','deptht','y','x'))
Wsamw = data.createVariable('Wtrans_SAMW','f8',('time_counter','deptht','y','x'))

Uaaiw.units = 'm3/s east'
Vaaiw.units = 'm3/s north'
Waaiw.units = 'm3/s up'
Usamw.units = 'm3/s east'
Vsamw.units = 'm3/s north'
Wsamw.units = 'm3/s up'

data.variables['Utrans_AAIW'][...] = u_mix_aaiw_wml
data.variables['Vtrans_AAIW'][...] = v_mix_aaiw_wml
data.variables['Wtrans_AAIW'][...] = w_mix_aaiw_wml
data.variables['Utrans_SAMW'][...] = u_mix_samw_wml
data.variables['Vtrans_SAMW'][...] = v_mix_samw_wml
data.variables['Wtrans_SAMW'][...] = w_mix_samw_wml


data.close()


#%%

os.remove('properties_at_modewater_outcrops_3D.nc')
shutil.copyfile(inputfile, 'properties_at_modewater_outcrops_3D.nc')

data = nc.Dataset('properties_at_modewater_outcrops_3D.nc','a')

Nitaaiw4d = data.createVariable('Nitrif_AAIW','f8',('time_counter','deptht','y','x'))
Nitsamw4d = data.createVariable('Nitrif_SAMW','f8',('time_counter','deptht','y','x'))
NO3aaiw4d = data.createVariable('Nitrat_AAIW','f8',('time_counter','deptht','y','x'))
NO3samw4d = data.createVariable('Nitrat_SAMW','f8',('time_counter','deptht','y','x'))
Silaaiw4d = data.createVariable('Silica_AAIW','f8',('time_counter','deptht','y','x'))
Silsamw4d = data.createVariable('Silica_SAMW','f8',('time_counter','deptht','y','x'))
Volaaiw4d = data.createVariable('Volume_AAIW','f8',('time_counter','deptht','y','x'))
Volsamw4d = data.createVariable('Volume_SAMW','f8',('time_counter','deptht','y','x'))

Nitaaiw4d.units = 'mol/m3/s'
Nitsamw4d.units = 'mol/m3/s'
NO3aaiw4d.units = 'mmol/m3'
NO3samw4d.units = 'mmol/m3'
Silaaiw4d.units = 'mmol/m3'
Silsamw4d.units = 'mmol/m3'
Volaaiw4d.units = 'm3'
Volsamw4d.units = 'm3'

data.variables['Nitrif_AAIW'][...] = nit_mix_aaiw_wml
data.variables['Nitrif_SAMW'][...] = nit_mix_samw_wml
data.variables['Nitrat_AAIW'][...] = no3_mix_aaiw_wml
data.variables['Nitrat_SAMW'][...] = no3_mix_samw_wml
data.variables['Silica_AAIW'][...] = si_mix_aaiw_wml
data.variables['Silica_SAMW'][...] = si_mix_samw_wml
data.variables['Volume_AAIW'][...] = vol_mix_aaiw_wml
data.variables['Volume_SAMW'][...] = vol_mix_samw_wml

data.close()


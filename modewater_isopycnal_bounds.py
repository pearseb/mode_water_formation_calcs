# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 12:11:47 2016

Purpose
-------
    Script to find isopycnal bounds of SAMW and AAIW using the annual and zonal average of salinity at 30 deg S.

@author: pearseb
"""

#%% imports

import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc


#%% Model data

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\model output\\calculate AAIW formation rates')

data = nc.Dataset('salinity_regular2x1.nc','r')
sal = data.variables['SAL_REG'][:,:,:]
data = nc.Dataset('potentialdensity_regular2x1.nc','r')
rho = data.variables['RHO_REG'][:,:,:]

#mnth = data.variables['MONTH'][:]
deps = data.variables['DEPTHT'][:]
dep_bnds = data.variables['DEPTHT_bnds'][:]
lats = data.variables['YAX'][:]
lons = data.variables['XAX'][:]

# get zonal averages

sal_zo = np.ma.mean(sal,axis=2)
rho_zo = np.ma.mean(rho,axis=2)


#%% get obs

data = nc.Dataset('salinity_regular2x1_woa18.nc','r')
sal_woa = data.variables['SAL'][0,...]
data = nc.Dataset('potentialdensity_regular2x1_woa18.nc','r')
rho_woa = data.variables['RHO'][0,...]

sal_woa_zo = np.ma.mean(sal_woa,axis=2)
rho_woa_zo = np.ma.mean(rho_woa,axis=2)



#%% interpolate depth of minimum

from scipy.interpolate import interp1d

f0 = interp1d(deps, sal_woa_zo[:,60], kind='cubic')
f1 = interp1d(deps, sal_zo[:,60], kind='cubic')

deps2 = np.arange(12.5,2201,0.1)

obs_levmin = np.where(f0(deps2) == np.min(f0(deps2)))[0]
levmin = np.where(f1(deps2) == np.min(f1(deps2)))[0]


print("Obs salinity minimum and depth = ", np.min(f0(deps2)), deps2[obs_levmin])
print("NEMO salinity minimum and depth = ", np.min(f1(deps2)), deps2[levmin])


#depths = np.array([1148, 629, 942, 856, 1067, 1126])
#print "model mean depth", np.mean(depths), np.std(depths)


#%% find the full width, half maximum to define the depth range of intermediate water

# 1. find the salinity of the deep ocean 
obs_saldeep = f0(2200) 
saldeep = f1(2200)

# 2. ignore all salinities greater than the deep salinity
obs_sal = f0(deps2)
obs_sal[obs_sal > obs_saldeep] = np.nan
sal = f1(deps2)
sal[sal > saldeep] = np.nan

# 3. find half maximum (minimum in this case)
obs_sal_hm = np.nanmin(obs_sal) + (np.nanmin(obs_sal) - obs_saldeep)*-0.5
sal_hm = np.nanmin(sal) + (np.nanmin(sal) - saldeep)*-0.5

# 4. exclude all salinities outside of the half maximum on deep side of profile
obs_sal_fwhm = obs_sal
obs_sal_fwhm[obs_sal > obs_sal_hm] = np.nan
sal_fwhm = sal
sal_fwhm[sal > sal_hm] = np.nan


# 5. find the depth interval of the half maximum
obs_depths = deps2[~np.isnan(obs_sal_fwhm)]
depths = deps2[~np.isnan(sal_fwhm)]

print(obs_depths)
print(depths)


#%% generate colour scheme ###
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
 
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
for i in range(len(tableau20)):
    r, g, b = tableau20[i]
    tableau20[i] = (r / 255., g / 255., b / 255.)
    

#%% make the plots

fig = plt.figure(facecolor='w', figsize=(10,6))

ax1 = plt.subplot(1,1,1)
ax1.spines["top"].set_visible(True)
ax1.spines["bottom"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.spines["left"].set_visible(True)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.tick_params(axis="both", which="both", bottom="off", top="on",
                labelbottom="off", labeltop='on', left="on", right="off", labelleft="on")
plt.scatter(sal_woa_zo[:,60], deps, color='k', marker='o')
plt.scatter(sal_zo[:,60], deps, color=tableau20[0], marker='o')
plt.ylabel('Depth (m)', family='serif', fontsize=14)
plt.xlabel('Salinity (psu) at 30$^{\circ}$S', family='serif', fontsize=14)
ax1.xaxis.set_label_position('top')
plt.ylim(2200,0)
plt.yticks(np.arange(2000,0,-200), np.arange(2000,0,-200), family='serif')
plt.xticks(np.arange(34.2,35.8,0.4), np.arange(34.2,35.8,0.4), family='serif')

plt.plot(f0(deps2), deps2, color='k',linestyle='-',linewidth=1, label='WOA (932 m)')
plt.plot(f1(deps2), deps2, color=tableau20[0],linestyle='-',linewidth=1, label='NEMO (1106 m)')

plt.plot(obs_sal_fwhm, deps2, color='k',linestyle=':',linewidth=3)
plt.plot(sal_fwhm, deps2, color=tableau20[0],linestyle=':',linewidth=3)

plt.legend(loc='center right',ncol=1,labelspacing=0.05,fontsize='large',frameon=False)


#%%

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\collaborations\\Robyn & Raja - Nitrification\\figures')
plt.savefig('fig-salmin.png',dpi=300,bbox_inches='tight')
plt.savefig('fig-salmin.pdf',dpi=300,bbox_inches='tight')
os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\collaborations\\Robyn & Raja - Nitrification\\figures\\trans')
plt.savefig('fig-salmin_trans.png',dpi=300,bbox_inches='tight', transparent=True)


#%% now find the isopycnal bounds


f11 = interp1d(deps, rho_zo[:,60], kind='cubic')
f22 = interp1d(deps, rho_woa_zo[:,60], kind='cubic')

isos = f11(depths)
isos_obs = f22(obs_depths)

print("Obs isopycnals of AAIW", np.min(isos_obs), np.max(isos_obs))
print("NEMO isopycnals of AAIW", np.min(isos), np.max(isos))

isopycnals = np.array([np.min(isos), np.max(isos)])
# correct isopycnals of AAIW to be a maximum of 1027.55, so that we do not include UCDW
isopycnals[1] = 1027.55

os.chdir('C:\\Users\\pearseb\\Dropbox\\PostDoc\\model output\\calculate AAIW formation rates')
np.save('modewater_isopycnal_bounds', isopycnals)

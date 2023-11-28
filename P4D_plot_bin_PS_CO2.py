# Author: Kier von Konigslow
# Date: May 2017
#
# Description: The purpose of this file is to plot the saturated PS/CO2 binary mixture information
#			   based on both experimental and theoretical data.
#

import os,sys,math,csv,matplotlib.pyplot as plt,numpy as npy
from p_params import *
from s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from P4D_plotting_parameters import *
from calculatePureVariables import calculateNewMolecularParameters
from wrapperFunctions import calculateBinarySolubilitySwelling
from calculateBinaryResidual import calculateBinarySSQ

P0 = npy.linspace(min(P0_X),max(P0_X),20)

verbose = False

#===================================================================================
#Setting global plotting parameters.
#===================================================================================

#Defining axes
solubility_axes = [5,23,0.01,0.13]
swelling_axes = [5,23,1.0,1.17]

#Defining labels
dhv_label = r'Present theory'
tsl_label = '2-parameter SL-EOS'
osl_label = '1-parameter SL-EOS'

dhv_color = (0.0,0.0,0.0,1.0)
tsl_color = (0.0,0.0,0.0,1.0)
osl_color = (0.0,0.0,0.0,1.0)

T453_color = (0.0,0.0,0.0,1.0)
T473_color = (0.0,0.0,0.0,1.0)
T493_color = (0.0,0.0,0.0,1.0)

#Defining linetype
dhv_line = '-'
tsl_line = ':'
osl_line = '--'

#Markers
mark453K = 'o'
mark463K = '^'
mark473K = 's'

#Linestules
ls453K = '-'
ls463K = '--'
ls473K = ':'

#===================================================================================
#Determining theoretical solubility/swelling.
#===================================================================================

#Set = [ zps , delta ]

#Fit to all data, X only.
solu = [1.02127006,0.88102354]

zeta = solu[0]
delta = solu[1]
print('DHV mixture parameters zeta = {} and delta = {}.'.format(zeta,delta))

gammas,vhs,epsilons = calculateNewMolecularParameters(Psstar,Tsstar,Rsstar,Ms)
vh = delta*vhs/NA
print('The hole volume is vh = {}.'.format(vh))

Pmin = min(P0_X)
Pmax = max(P0_X)
Tmin = min(T0_X)
Tmax = max(T0_X)
print('The pressure range is {}-{}MPa and the temperature range is {}-{}K.'.format(Pmin,Pmax,Tmin,Tmax))

result = calculateBinarySolubilitySwelling('CHV',P0,403.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
ms_403_DHV = result[2]
Swn_403_DHV = result[3]
result = calculateBinarySolubilitySwelling('CHV',P0,423.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
ms_423_DHV = result[2]
Swn_423_DHV = result[3]
result = calculateBinarySolubilitySwelling('CHV',P0,463.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
ms_463_DHV = result[2]
Swn_463_DHV = result[3]

result = calculateBinarySolubilitySwelling('CHV',P0,403.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose,simplified=False)
Swa_403_DHV = result[3]
result = calculateBinarySolubilitySwelling('CHV',P0,423.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose,simplified=False)
Swa_423_DHV = result[3]
result = calculateBinarySolubilitySwelling('CHV',P0,463.0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose,simplified=False)
Swa_463_DHV = result[3]

#===================================================================================
#Plotting the PS/CO2 mixture results.
#===================================================================================

#Plotting the solubility of the PS+CO2 mixture.
figX = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_X_403,X0_X_403,'k',marker=mark453K,ls='',label='Saturated PS/CO2, 403K',ms=markersize)
plt.plot(P0_X_423,X0_X_423,'k',marker=mark463K,ls='',label='423K',ms=markersize)
plt.plot(P0_X_463,X0_X_463,'k',marker=mark473K,ls='',label='463K',ms=markersize)
plt.plot(P0,ms_403_DHV,'k',ls=ls453K,label='Present theory 403K',lw=linewidth)
plt.plot(P0,ms_423_DHV,'k',ls=ls463K,label='Present theory 423K',lw=linewidth)
plt.plot(P0,ms_463_DHV,'k',ls=ls473K,label='Present theory 463K',lw=linewidth)
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel('Solvent solubility $\chi_s$ (g/g)',fontsize=axis_size)
legend = plt.legend(loc=2,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis(solubility_axes)
plt.tight_layout()
figX.savefig('../'+output_folder+r'\bin_PS_CO2_Solubility'+img_extension,dpi=img_dpi)

#Plotting the swelling of the PS+CO2 mixture, neat polymer.
figSN = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_S_403,S0_S_403,'k',marker=mark453K,ls='',label='Saturated PS/CO2, 403K',ms=markersize)
plt.plot(P0_S_423,S0_S_423,'k',marker=mark463K,ls='',label='423K',ms=markersize)
plt.plot(P0_S_463,S0_S_463,'k',marker=mark473K,ls='',label='463K',ms=markersize)
plt.plot(P0,Swn_403_DHV,'k',ls=ls453K,label='Theory 403K',lw=linewidth)
plt.plot(P0,Swn_423_DHV,'k',ls=ls463K,label='Theory 423K',lw=linewidth)
plt.plot(P0,Swn_463_DHV,'k',ls=ls473K,label='Theory 463K',lw=linewidth)
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel('Swelling Ratio $S_W$ (cm$^3$/cm$^3$)',fontsize=axis_size)
legend = plt.legend(loc=2,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis(swelling_axes)
plt.tight_layout()
figSN.savefig('../'+output_folder+r'\bin_PS_CO2_Swelling_neat'+img_extension,dpi=img_dpi)

#Plotting the swelling of the PS+CO2 mixture, air compensation.
figSA = plt.figure(num=None, figsize=(img_vert, img_hori), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_S_403,S0_S_403,'k',marker=mark453K,ls='',label='Saturated PS/CO2, 403K',ms=markersize)
plt.plot(P0_S_423,S0_S_423,'k',marker=mark463K,ls='',label='423K',ms=markersize)
plt.plot(P0_S_463,S0_S_463,'k',marker=mark473K,ls='',label='463K',ms=markersize)
plt.plot(P0,Swa_403_DHV,'k',ls=ls453K,label='Theory 403K',lw=linewidth)
plt.plot(P0,Swa_423_DHV,'k',ls=ls463K,label='Theory 423K',lw=linewidth)
plt.plot(P0,Swa_463_DHV,'k',ls=ls473K,label='Theory 463K',lw=linewidth)
plt.xlabel('Pressure $P$ (MPa)',fontsize=axis_size)
plt.ylabel('Swelling Ratio $S_W$ (cm$^3$/cm$^3$)',fontsize=axis_size)
legend = plt.legend(loc=2,fontsize=size,numpoints=1)
if not showbox:
	legend.get_frame().set_facecolor('none')
	legend.get_frame().set_linewidth(0.0)
plt.axis(swelling_axes)
plt.tight_layout()
figSA.savefig('../'+output_folder+r'\bin_PS_CO2_Swelling_air'+img_extension,dpi=img_dpi)
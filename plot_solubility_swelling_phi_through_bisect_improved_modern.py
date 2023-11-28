# Author: Kier von Konigslow
# Date: November 2017
#
import os,sys,math,csv,matplotlib.pyplot as plt,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from calculatePureVariables import calculateNewMolecularParameters
# from wrapperFunctions import calculateBinarySolubilitySwelling
# from calculateBinaryResidual import calculateBinarySSQ
from checkResults import *
from sympy import *
import warnings
import cmath
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from Parameters_of_Different_Polymers import *
from All_Functions import *
from Parameters_for_Mixtures_and_Tg import *

###########################################################################
###########################################################################

Polymer_Type='PC'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}

Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)

P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

isotherm_label=1
number_of_isotherm=1
isotherm_scanner=T0_X[0]

for k in range(0,10):
	exec "len_of_isotherm%s=0" %(isotherm_label)
	isotherm_label=isotherm_label+1

isotherm_label=1

for k in range(0,len(T0_X)):
	if isotherm_scanner==T0_X[k]:
		exec "len_of_isotherm%s=len_of_isotherm%s+1" %(isotherm_label,isotherm_label)
	elif isotherm_scanner!=T0_X[k]:
		number_of_isotherm = number_of_isotherm+1
		isotherm_scanner=T0_X[k]
		isotherm_label=isotherm_label+1
		exec "len_of_isotherm%s=len_of_isotherm%s+1" %(isotherm_label,isotherm_label)


print 'number of isotherms are:',number_of_isotherm

isotherm_label=1
for i in range(0,number_of_isotherm):
	exec "P0_X_T%s=npy.zeros(len_of_isotherm%s)" %(isotherm_label,isotherm_label)
	exec "T0_X_T%s=npy.zeros(len_of_isotherm%s)" %(isotherm_label,isotherm_label)
	exec "X0_X_T%s=npy.zeros(len_of_isotherm%s)" %(isotherm_label,isotherm_label)
	isotherm_label=isotherm_label+1

isotherm_label=1
k=0
for i in range(0,len(T0_X)):
	exec "P0_X_T%s[k]=P0_X[i]" %(isotherm_label)
	exec "T0_X_T%s[k]=T0_X[i]" %(isotherm_label)
	exec "X0_X_T%s[k]=X0_X[i]" %(isotherm_label)
	k=k+1
	if i+2<=len(T0_X):
		if T0_X[i+1]!=T0_X[i]:
			k=0
			isotherm_label=isotherm_label+1

isotherm_label=1
for i in range(0,number_of_isotherm):
	exec "print P0_X_T%s" %(isotherm_label)
	exec "print T0_X_T%s" %(isotherm_label)
	exec "print X0_X_T%s" %(isotherm_label)
	isotherm_label=isotherm_label+1

Kier=False
Hassan=False  
Hassan_Var_Vol=True  
Condo=False  
Condo_Original=False 

kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Kier':Kier,'Hassan':Hassan,'Hassan_Var_Vol':Hassan_Var_Vol,'Condo':Condo,'Condo_Original':True}
cepsilon_s,cepsilon_p,cz,czeta,epsilon_p,g,x,delta,zeta=Parameters_for_Mixtures_and_Tg(**kwargs)
cdelta=100.0

xS_infty=x*(Ppstar/(Tpstar*Rpstar))*(1+ln(1+g))

Rescaling=True
if Rescaling:
	vh_pure_p=kB*Tpstar/Ppstar
	vh_pure_s=kB*Tsstar/Psstar
	v_r=vh_pure_p
	Tpstar,Tsstar=rescaling(Tpstar,Tsstar,vh_pure_p,vh_pure_s,v_r)

# alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
# alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)
# pfactor=Mp/alpha_p0
# sfactor=Ms/alpha_s0
# print pfactor
# print sfactor
# print alpha_p0
# print alpha_s0
# print vh_pure_p
# print vh_pure_s
#1.03 #1.20
#1.087, 1.115, 1.100
zeta=	2.0	#0.932 	#0.94 	#10.061	#1.068		#1.08820786	#PS zeta= 1.124 and delta=0.97423316 at 2.0 bar gives retrograde  #1.100
delta=	1.0		#0.54 	#1.3	#1.5	#0.97423316	#0.97423316

#Both PS and PC and PMMA are not showing retrograde for vhm=constant for both upper and lower isotherm fit:
# upperisothermfit=zeta=1.025 and lowerisothermfit=1.16 of PC with v_r=vhp and delta=1.0 and Hassan constant vhm=vhp
# upperisothermfit=zeta=1.18 and lowerisothermfit=1.30 of PS with v_r=vhp and delta=1.0 and Hassan constant vhm=vhp
# upperisothermfit=zeta=1.1 and lowerisothermfit=1.13 of PMMA with v_r=vhp and delta=1.0 and Hassan constant vhm=vhp

# Thus Hassan Var Volume is a flop method as it is giving retrograde in both PS and PC and PMMA for both upper and lower isotherm fit:
# upperisothermfit=zeta=1.025 and lowerisothermfit=1.16 of PC with v_r=vhp and delta=useless and Hassan Var Volume
# upperisothermfit=zeta=1.18 and lowerisothermfit=1.30 of PS with v_r=vhp and delta=useless and Hassan Var Volume
# upperisothermfit=zeta=1.1 and lowerisothermfit=1.13 of PMMA with v_r=vhp and delta=useless and Hassan Var Volume

# 1.10694978 Condo Solubility formula

Isotherms=False
Plot_Solubility=Isotherms

Isobars=True
Entropy=Isobars
Plot_Phi=Isobars

Plot_Swelling=False		#I have not corrected Swelling Error will come.
Plot_S2=True			

forward=True
backward=False 

if Isotherms:
	P0 = npy.linspace(min(P0_X),max(P0_X),4)
	# P0 = npy.linspace(0.301325,30,20)
	T1=T0_X_T1[0]	#403	#290
	T2=0.0#T0_X_T2[0]	#423	#304
	T3=0.0#T0_X_T3[0]	#463	#350
	T4=0.0#T0_X_T4[0]	#423	#304
	T5=0.0#T0_X_T5[0]	#463	#350
	T6=0.0#T0_X_T6[0]	#423	#304
	T7=0.0#T0_X_T7[0]	#463	#350
	T8=T0_X_T8[0]	#463	#350

if Isobars:
	T0 = npy.linspace(200,500,5)		#max: 1400000  #Small pressure ==> entropy max reaches at smaller temperature
	P1=0.1013250#0.101325
	P2=4.0#30.0
	P3=9.0#50.0


###############################################################################################
###############################################################################################

#===================================================================================
#Determining theoretical solubility/swelling.
#===================================================================================
# print('DHV mixture parameters zeta = {} and delta = {}.'.format(zeta,delta))
# gammas,vhs,epsilons = calculateNewMolecularParameters(Psstar,Tsstar,Rsstar,Ms)
# vh = delta*vhs/NA
# print('The hole volume is vh = {}.'.format(vh))
Pmin = min(P0_X)
Pmax = max(P0_X)
Tmin = min(T0_X)
Tmax = max(T0_X)
print('The pressure range is {}-{}MPa and the temperature range is {}-{}K.'.format(Pmin,Pmax,Tmin,Tmax))
############################################################################################################
if Kier or Hassan or Condo or Hassan_Var_Vol:
	if Isotherms:

		if T1!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T1,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T1_DHV = result[2]
			Sw_T1_DHV = result[3]
			phip_T1_DHV = result[4]
			phis_T1_DHV = result[5]
			Rtilde_T1_DHV = result[6]
			phip0_T1_DHV = result[7]
			phis0_T1_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T1,phip_T1_DHV,phis_T1_DHV,phip0_T1_DHV,phis0_T1_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T1_DHV = properties[2]
				S_2_T1_DHV = properties[3]
			######################################
		if T2!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T2,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T2_DHV = result[2]
			Sw_T2_DHV = result[3]
			phip_T2_DHV = result[4]
			phis_T2_DHV = result[5]
			Rtilde_T2_DHV = result[6]
			phip0_T2_DHV = result[7]
			phis0_T2_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T2,phip_T2_DHV,phis_T2_DHV,phip0_T2_DHV,phis0_T2_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T2_DHV = properties[2]
				S_2_T2_DHV = properties[3]
			######################################
		if T3!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T3,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T3_DHV = result[2]
			Sw_T3_DHV = result[3]
			phip_T3_DHV = result[4]
			phis_T3_DHV = result[5]
			Rtilde_T3_DHV = result[6]
			phip0_T3_DHV = result[7]
			phis0_T3_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T3,phip_T3_DHV,phis_T3_DHV,phip0_T3_DHV,phis0_T3_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T3_DHV = properties[2]
				S_2_T3_DHV = properties[3]
			#######################################
		if T4!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T4,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T4_DHV = result[2]
			Sw_T4_DHV = result[3]
			phip_T4_DHV = result[4]
			phis_T4_DHV = result[5]
			Rtilde_T4_DHV = result[6]
			phip0_T4_DHV = result[7]
			phis0_T4_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T4,phip_T4_DHV,phis_T4_DHV,phip0_T4_DHV,phis0_T4_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T4_DHV = properties[2]
				S_2_T4_DHV = properties[3]
			#######################################
		if T5!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T5,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T5_DHV = result[2]
			Sw_T5_DHV = result[3]
			phip_T5_DHV = result[4]
			phis_T5_DHV = result[5]
			Rtilde_T5_DHV = result[6]
			phip0_T5_DHV = result[7]
			phis0_T5_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T5,phip_T5_DHV,phis_T5_DHV,phip0_T5_DHV,phis0_T5_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T5_DHV = properties[2]
				S_2_T5_DHV = properties[3]
			#######################################
		if T6!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T6,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T6_DHV = result[2]
			Sw_T6_DHV = result[3]
			phip_T6_DHV = result[4]
			phis_T6_DHV = result[5]
			Rtilde_T6_DHV = result[6]
			phip0_T6_DHV = result[7]
			phis0_T6_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T6,phip_T6_DHV,phis_T6_DHV,phip0_T6_DHV,phis0_T6_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T6_DHV = properties[2]
				S_2_T6_DHV = properties[3]
			#######################################
		if T7!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T7,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T7_DHV = result[2]
			Sw_T7_DHV = result[3]
			phip_T7_DHV = result[4]
			phis_T7_DHV = result[5]
			Rtilde_T7_DHV = result[6]
			phip0_T7_DHV = result[7]
			phis0_T7_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T7,phip_T7_DHV,phis_T7_DHV,phip0_T7_DHV,phis0_T7_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T7_DHV = properties[2]
				S_2_T7_DHV = properties[3]
			#######################################
		if T8!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P0,T8,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_T8_DHV = result[2]
			Sw_T8_DHV = result[3]
			phip_T8_DHV = result[4]
			phis_T8_DHV = result[5]
			Rtilde_T8_DHV = result[6]
			phip0_T8_DHV = result[7]
			phis0_T8_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P0,T8,phip_T8_DHV,phis_T8_DHV,phip0_T8_DHV,phis0_T8_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_T8_DHV = properties[2]
				S_2_T8_DHV = properties[3]
			#######################################

	if Isobars:

		if P1!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P1,T0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_P1_DHV = result[2]
			Sw_P1_DHV = result[3]
			phip_P1_DHV = result[4]
			phis_P1_DHV = result[5]
			Rtilde_P1_DHV = result[6]
			phip0_P1_DHV = result[7]
			phis0_P1_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P1,T0,phip_P1_DHV,phis_P1_DHV,phip0_P1_DHV,phis0_P1_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_P1_DHV = properties[2]
				S_2_P1_DHV = properties[3]
			########################################
		if P2!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P2,T0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_P2_DHV = result[2]
			Sw_P2_DHV = result[3]
			phip_P2_DHV = result[4]
			phis_P2_DHV = result[5]
			Rtilde_P2_DHV = result[6]
			phip0_P2_DHV = result[7]
			phis0_P2_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P2,T0,phip_P2_DHV,phis_P2_DHV,phip0_P2_DHV,phis0_P2_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_P2_DHV = properties[2]
				S_2_P2_DHV = properties[3]
			########################################
		if P3!=0.0:
			result = calculateBinarySolubilitySwelling('CHV',P3,T0,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol, forward=forward,backward=backward)
			Xs_P3_DHV = result[2]
			Sw_P3_DHV = result[3]
			phip_P3_DHV = result[4]
			phis_P3_DHV = result[5]
			Rtilde_P3_DHV = result[6]
			phip0_P3_DHV = result[7]
			phis0_P3_DHV = result[8]
			
			if Entropy:
				properties=calculateThermodynamicVariables(P3,T0,phip_P3_DHV,phis_P3_DHV,phip0_P3_DHV,phis0_P3_DHV,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
				S_1_P3_DHV = properties[2]
				S_2_P3_DHV = properties[3]
			########################################

if Condo_Original:
	if Isotherms:
		if T1!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P0,T1,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
			Xs_T1_Condo = result[2]
			Sw_T1_Condo = result[3]
			phip_T1_Condo = result[4]
			phis_T1_Condo = result[5]
			Rtilde_T1_Condo = result[6]
			phip0_T1_Condo = result[7]
			phis0_T1_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P0,T1,phip_T1_Condo,phis_T1_Condo,phip0_T1_Condo,phis0_T1_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_T1_Condo = properties[2]
				S_2_T1_Condo = properties[3]
				######################################
		if T2!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P0,T2,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
			Xs_T2_Condo = result[2]
			Sw_T2_Condo = result[3]
			phip_T2_Condo = result[4]
			phis_T2_Condo = result[5]
			Rtilde_T2_Condo = result[6]
			phip0_T2_Condo = result[7]
			phis0_T2_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P0,T2,phip_T2_Condo,phis_T2_Condo,phip0_T2_Condo,phis0_T2_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_T2_Condo = properties[2]
				S_2_T2_Condo = properties[3]
			######################################
		if T3!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P0,T3,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
			Xs_T3_Condo = result[2]
			Sw_T3_Condo = result[3]
			phip_T3_Condo = result[4]
			phis_T3_Condo = result[5]
			Rtilde_T3_Condo = result[6]
			phip0_T3_Condo = result[7]
			phis0_T3_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P0,T3,phip_T3_Condo,phis_T3_Condo,phip0_T3_Condo,phis0_T3_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_T3_Condo = properties[2]
				S_2_T3_Condo = properties[3]
		#######################################

	if Isobars:

		if P1!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P1,T0,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
			Xs_P1_Condo = result[2]
			Sw_P1_Condo = result[3]
			phip_P1_Condo = result[4]
			phis_P1_Condo = result[5]
			Rtilde_P1_Condo = result[6]
			phip0_P1_Condo = result[7]
			phis0_P1_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P1,T0,phip_P1_Condo,phis_P1_Condo,phip0_P1_Condo,phis0_P1_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_P1_Condo = properties[2]
				S_2_P1_Condo = properties[3]
			########################################
		if P2!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P2,T0,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
			Xs_P2_Condo = result[2]
			Sw_P2_Condo = result[3]
			phip_P2_Condo = result[4]
			phis_P2_Condo = result[5]
			Rtilde_P2_Condo = result[6]
			phip0_P2_Condo = result[7]
			phis0_P2_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P2,T0,phip_P2_Condo,phis_P2_Condo,phip0_P2_Condo,phis0_P2_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_P2_Condo = properties[2]
				S_2_P2_Condo = properties[3]
			########################################
		if P3!=0.0:
			result = calculateBinarySolubilitySwelling('Condo_Original',P3,T0,Mp,Ms,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar, forward=forward,backward=backward)
			Xs_P3_Condo = result[2]
			Sw_P3_Condo = result[3]
			phip_P3_Condo = result[4]
			phis_P3_Condo = result[5]
			Rtilde_P3_Condo = result[6]
			phip0_P3_Condo = result[7]
			phis0_P3_Condo = result[8]
			
			if Entropy:
				properties=calculateCondoThermodynamicVariables_Original(P3,T0,phip_P3_Condo,phis_P3_Condo,phip0_P3_Condo,phis0_P3_Condo,Mp,Ms,z=cz,epsilon_s=cepsilon_s,epsilon_p=cepsilon_p,zeta=czeta,delta=cdelta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
				S_1_P3_Condo = properties[2]
				S_2_P3_Condo = properties[3]
			########################################

#===================================================================================
#Plotting the PS/CO2 mixture results.
#===================================================================================

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'Plot_Thesis'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Setting font size
axis_size = 20
title_size = 20
size = 8
label_size = 16
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Defining axes
solubility_axes = [5,25,0.0,0.13]
swelling_axes = [5,25,1.0,1.17]
phi_axes = [5,25,0.0,1.0]
TD_axes = [5,25,0.0,5.0]

#Defining labels
dhv_label = r'Present theory'
dhv_color = (0.0,0.0,0.0,1.0)
T1_color = (0.0,0.0,0.0,1.0)
T2_color = (0.0,0.0,0.0,1.0)
T3_color = (0.0,0.0,0.0,1.0)

#Defining linetype
dhv_line = '-'

#Markers
markT1K = 'o'
markT2K = 's'
markT3K = '>'
markT4K = '<'
markT5K = 'D'
markT6K = 'H'
markT7K = 'P'
markT8K = 'X'

#Linestyles
lsT1K = '-'
lsT2K = '--'
lsT3K = ':'

#General line properties.
linewidth = 2
markersize = 8
if Kier or Hassan or Condo or Hassan_Var_Vol:
	if Isotherms:

		if Plot_Phi:
			#Plotting the phi's of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(P0,phip_T1_DHV,'r',ls=lsT1K,label='phi_p_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phip_T2_DHV,'r',ls=lsT2K,label='phi_p_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phip_T3_DHV,'r',ls=lsT3K,label='phi_p_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phis_T1_DHV,'m',ls=lsT1K,label='phi_s_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phis_T2_DHV,'m',ls=lsT2K,label='phi_s_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phis_T3_DHV,'m',ls=lsT3K,label='phi_s_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,Rtilde_T1_DHV,'b',ls=lsT1K,label='Rtilde_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Rtilde_T2_DHV,'b',ls=lsT2K,label='Rtilde_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Rtilde_T3_DHV,'b',ls=lsT3K,label='Rtilde_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phip0_T1_DHV,'k',ls=lsT1K,label='phi_p0_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phip0_T2_DHV,'k',ls=lsT2K,label='phi_p0_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phip0_T3_DHV,'k',ls=lsT3K,label='phi_p0_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phis0_T1_DHV,'y',ls=lsT1K,label='phi_s0_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phis0_T2_DHV,'y',ls=lsT2K,label='phi_s0_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phis0_T3_DHV,'y',ls=lsT3K,label='phi_s0_{} K'.format(T3),lw=linewidth)

			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('phi',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(phi_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Plot_Solubility:
			#Plotting the solubility of the PS+CO2 mixture.
			figX = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(P0_X_T1,X0_X_T1,'r',marker=markT1K,ls='',label='Saturated mixture at {}K'.format(T0_X_T1[0]),ms=markersize)
			plt.plot(P0_X_T2,X0_X_T2,'b',marker=markT2K,ls='',label='{}K'.format(T0_X_T2[0]),ms=markersize)
			plt.plot(P0_X_T3,X0_X_T3,'g',marker=markT3K,ls='',label='{}K'.format(T0_X_T3[0]),ms=markersize)
			plt.plot(P0_X_T4,X0_X_T4,'m',marker=markT4K,ls='',label='{}K'.format(T0_X_T4[0]),ms=markersize)
			plt.plot(P0_X_T5,X0_X_T5,'c',marker=markT5K,ls='',label='{}K'.format(T0_X_T5[0]),ms=markersize)
			plt.plot(P0_X_T6,X0_X_T6,'y',marker=markT6K,ls='',label='{}K'.format(T0_X_T6[0]),ms=markersize)
			plt.plot(P0_X_T7,X0_X_T7,'k',marker=markT7K,ls='',label='{}K'.format(T0_X_T7[0]),ms=markersize)
			plt.plot(P0_X_T8,X0_X_T8,'r',marker=markT8K,ls='',label='{}K'.format(T0_X_T8[0]),ms=markersize)

			plt.plot(P0,Xs_T1_DHV,'r',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Xs_T2_DHV,'b',ls=lsT2K,label='{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Xs_T3_DHV,'g',ls=lsT3K,label='{} K'.format(T3),lw=linewidth)
			if T4!=0.0:
				plt.plot(P0,Xs_T4_DHV,'m',ls=lsT3K,label='{} K'.format(T4),lw=linewidth)
			if T5!=0.0:
				plt.plot(P0,Xs_T5_DHV,'c',ls=lsT3K,label='{} K'.format(T5),lw=linewidth)
			if T6!=0.0:
				plt.plot(P0,Xs_T6_DHV,'y',ls=lsT3K,label='{} K'.format(T6),lw=linewidth)
			if T7!=0.0:
				plt.plot(P0,Xs_T7_DHV,'k',ls=lsT3K,label='{} K'.format(T7),lw=linewidth)
			if T8!=0.0:
				plt.plot(P0,Xs_T8_DHV,'r',ls=lsT3K,label='{} K'.format(T8),lw=linewidth)
			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Gas solubility $\chi_s$ (g/g)',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(solubility_axes)
			plt.subplots_adjust(bottom=0.3)
			figX.savefig('./'+output_folder+r'\bin_PC_CO2_Solubility'+img_extension,dpi=240)

		if Plot_Swelling:
			#Plotting the swelling of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(P0_S_T1,S0_S_T1,'k',marker=markT1K,ls='',label='Saturated PS/CO2 mixture at T1K',ms=markersize)
			plt.plot(P0_S_T2,S0_S_T2,'k',marker=markT2K,ls='',label='T2K',ms=markersize)
			plt.plot(P0_S_T3,S0_S_T3,'k',marker=markT3K,ls='',label='T3K',ms=markersize)
			plt.plot(P0,Sw_T1_DHV,'k',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Sw_T2_DHV,'k',ls=lsT2K,label='Present theory {} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Sw_T3_DHV,'k',ls=lsT3K,label='Present theory {} K'.format(T3),lw=linewidth)
			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Swelling ratio $S_W$ $(cm^3/cm^3)$',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			plt.axis(swelling_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Entropy:
			#Plotting the phi's of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(P0,S_1_T1_DHV,'r',ls=lsT1K,label='S_1_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,S_1_T2_DHV,'r',ls=lsT2K,label='S_1_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,S_1_T3_DHV,'r',ls=lsT3K,label='S_1_{} K'.format(T3),lw=linewidth)

			if Plot_S2:
			
				plt.plot(P0,S_2_T1_DHV,'m',ls=lsT1K,label='S_2_{} K'.format(T1),lw=linewidth)
				if T2!=0.0:
					plt.plot(P0,S_2_T2_DHV,'m',ls=lsT2K,label='S_2_{} K'.format(T2),lw=linewidth)
				if T3!=0.0:
					plt.plot(P0,S_2_T3_DHV,'m',ls=lsT3K,label='S_2_{} K'.format(T3),lw=linewidth)

			# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
			plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')
			# S_max=npy.max(S_1_P1)
			# print 'S_max is:', S_max
			Tg_line=xS_infty#0.317*0.98268#0.310707*0.8708171#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361
			plt.axhline(y=Tg_line,lw=0.5,color='k', linestyle='-.')

			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Entropy',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(TD_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

	if Isobars:

		if Plot_Phi:
			#Plotting the phi's of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,phip_P1_DHV,'r',ls=lsT1K,label='phi_p_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phip_P2_DHV,'r',ls=lsT2K,label='phi_p_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phip_P3_DHV,'r',ls=lsT3K,label='phi_p_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phis_P1_DHV,'m',ls=lsT1K,label='phi_s_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phis_P2_DHV,'m',ls=lsT2K,label='phi_s_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phis_P3_DHV,'m',ls=lsT3K,label='phi_s_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,Rtilde_P1_DHV,'b',ls=lsT1K,label='Rtilde_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Rtilde_P2_DHV,'b',ls=lsT2K,label='Rtilde_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Rtilde_P3_DHV,'b',ls=lsT3K,label='Rtilde_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phip0_P1_DHV,'k',ls=lsT1K,label='phi_p0_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phip0_P2_DHV,'k',ls=lsT2K,label='phi_p0_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phip0_P3_DHV,'k',ls=lsT3K,label='phi_p0_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phis0_P1_DHV,'y',ls=lsT1K,label='phi_s0_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phis0_P2_DHV,'y',ls=lsT2K,label='phi_s0_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phis0_P3_DHV,'y',ls=lsT3K,label='phi_s0_{} MPa'.format(P3),lw=linewidth)

			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('phi',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(phi_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Plot_Solubility:
			#Plotting the solubility of the PS+CO2 mixture.
			figX = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,Xs_P1_DHV,'k',ls=lsT1K,label='Present theory {} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Xs_P2_DHV,'k',ls=lsT2K,label='Present theory {} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Xs_P3_DHV,'k',ls=lsT3K,label='Present theory {} MPa'.format(P3),lw=linewidth)
			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('Gas solubility $\chi_s$ (g/g)',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			# plt.axis(solubility_axes)
			# figX.savefig('./'+output_folder+r'\bin_PS_CO2_Solubility'+img_extension,dpi=img_dpi)

		if Plot_Swelling:
			#Plotting the swelling of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(T0,Sw_P1_DHV,'k',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Sw_P2_DHV,'k',ls=lsT2K,label='Present theory {} K'.format(T2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Sw_P3_DHV,'k',ls=lsT3K,label='Present theory {} K'.format(T3),lw=linewidth)
			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Swelling ratio $S_W$ $(cm^3/cm^3)$',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			plt.axis(swelling_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Entropy:
			#Plotting the phi's of the PS+CO2 mixture.
			figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,S_1_P1_DHV,'r',ls=lsT1K,label='S_1_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,S_1_P2_DHV,'r',ls=lsT2K,label='S_1_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,S_1_P3_DHV,'r',ls=lsT3K,label='S_1_{} MPa'.format(P3),lw=linewidth)

			if Plot_S2:
				
				plt.plot(T0,S_2_P1_DHV,'m',ls=lsT1K,label='S_2_{} MPa'.format(P1),lw=linewidth)
				if P2!=0.0:
					plt.plot(T0,S_2_P2_DHV,'m',ls=lsT2K,label='S_2_{} MPa'.format(P2),lw=linewidth)
				if P3!=0.0:
					plt.plot(T0,S_2_P3_DHV,'m',ls=lsT3K,label='S_2_{} MPa'.format(P3),lw=linewidth)
			
			# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
			plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')
			# S_max=npy.max(S_1_P1)
			# print 'S_max is:', S_max
			Tg_line=xS_infty#0.317*0.98268 #0.310707*0.8708171#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361
			plt.axhline(y=Tg_line,lw=0.5,color='k', linestyle='-.')

			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('Entropy',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(TD_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

if Condo_Original:
	if Isotherms:

		if Plot_Phi:
			#Plotting the phi's of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(P0,phip_T1_Condo,'r',ls=lsT1K,label='phi_p_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phip_T2_Condo,'r',ls=lsT2K,label='phi_p_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phip_T3_Condo,'r',ls=lsT3K,label='phi_p_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phis_T1_Condo,'m',ls=lsT1K,label='phi_s_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phis_T2_Condo,'m',ls=lsT2K,label='phi_s_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phis_T3_Condo,'m',ls=lsT3K,label='phi_s_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,Rtilde_T1_Condo,'b',ls=lsT1K,label='Rtilde_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Rtilde_T2_Condo,'b',ls=lsT2K,label='Rtilde_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Rtilde_T3_Condo,'b',ls=lsT3K,label='Rtilde_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phip0_T1_Condo,'k',ls=lsT1K,label='phi_p0_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phip0_T2_Condo,'k',ls=lsT2K,label='phi_p0_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phip0_T3_Condo,'k',ls=lsT3K,label='phi_p0_{} K'.format(T3),lw=linewidth)

			plt.plot(P0,phis0_T1_Condo,'y',ls=lsT1K,label='phi_s0_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,phis0_T2_Condo,'y',ls=lsT2K,label='phi_s0_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,phis0_T3_Condo,'y',ls=lsT3K,label='phi_s0_{} K'.format(T3),lw=linewidth)

			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('phi',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(phi_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Plot_Solubility:
			#Plotting the solubility of the PS+CO2 mixture.
			figX_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(P0_X_T1,X0_X_T1,'k',marker=markT1K,ls='',label='Saturated mixture at {}K'.format(T0_X_T1[0]),ms=markersize)
			plt.plot(P0_X_T2,X0_X_T2,'k',marker=markT2K,ls='',label='{}K'.format(T0_X_T2[0]),ms=markersize)
			plt.plot(P0_X_T3,X0_X_T3,'k',marker=markT3K,ls='',label='{}K'.format(T0_X_T3[0]),ms=markersize)
			plt.plot(P0_X_T4,X0_X_T4,'k',marker=markT4K,ls='',label='{}K'.format(T0_X_T4[0]),ms=markersize)
			plt.plot(P0_X_T5,X0_X_T5,'k',marker=markT5K,ls='',label='{}K'.format(T0_X_T5[0]),ms=markersize)
			plt.plot(P0_X_T6,X0_X_T6,'k',marker=markT6K,ls='',label='{}K'.format(T0_X_T6[0]),ms=markersize)

			plt.plot(P0,Xs_T1_Condo,'k',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Xs_T2_Condo,'k',ls=lsT2K,label='Present theory {} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Xs_T3_Condo,'k',ls=lsT3K,label='Present theory {} K'.format(T3),lw=linewidth)

			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Gas solubility $\chi_s$ (g/g)',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(solubility_axes)
			# figX.savefig('./'+output_folder+r'\bin_PS_CO2_Solubility'+img_extension,dpi=img_dpi)

		if Plot_Swelling:
			#Plotting the swelling of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(P0_S_T1,S0_S_T1,'k',marker=markT1K,ls='',label='Saturated mixture at {}K'.format(T0_X_T1[0]),ms=markersize)
			plt.plot(P0_S_T2,S0_S_T2,'k',marker=markT2K,ls='',label='{}K'.format(T0_X_T2[0]),ms=markersize)
			plt.plot(P0_S_T3,S0_S_T3,'k',marker=markT3K,ls='',label='{}K'.format(T0_X_T3[0]),ms=markersize)
			# plt.plot(P0_S_T4,S0_S_T4,'k',marker=markT1K,ls='',label='{}K'.format(T0_X_T4[0]),ms=markersize)
			# plt.plot(P0_S_T5,S0_S_T5,'k',marker=markT2K,ls='',label='{}K'.format(T0_X_T5[0]),ms=markersize)
			# plt.plot(P0_S_T6,S0_S_T6,'k',marker=markT3K,ls='',label='{}K'.format(T0_X_T6[0]),ms=markersize)

			plt.plot(P0,Sw_T1_Condo,'k',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,Sw_T2_Condo,'k',ls=lsT2K,label='Present theory {} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,Sw_T3_Condo,'k',ls=lsT3K,label='Present theory {} K'.format(T3),lw=linewidth)
			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Swelling ratio $S_W$ $(cm^3/cm^3)$',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			plt.axis(swelling_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Entropy:
			#Plotting the phi's of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(P0,S_1_T1_Condo,'r',ls=lsT1K,label='S_1_{} K'.format(T1),lw=linewidth)
			if T2!=0.0:
				plt.plot(P0,S_1_T2_Condo,'r',ls=lsT2K,label='S_1_{} K'.format(T2),lw=linewidth)
			if T3!=0.0:
				plt.plot(P0,S_1_T3_Condo,'r',ls=lsT3K,label='S_1_{} K'.format(T3),lw=linewidth)

			if Plot_S2:
				
				plt.plot(P0,S_2_T1_Condo,'m',ls=lsT1K,label='S_2_{} K'.format(T1),lw=linewidth)
				if T2!=0.0:
					plt.plot(P0,S_2_T2_Condo,'m',ls=lsT2K,label='S_2_{} K'.format(T2),lw=linewidth)
				if T3!=0.0:
					plt.plot(P0,S_2_T3_Condo,'m',ls=lsT3K,label='S_2_{} K'.format(T3),lw=linewidth)

			# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
			plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')
			# S_max=npy.max(S_1_P1)
			# print 'S_max is:', S_max
			Tg_line=xS_infty#0.310707*0.8708171#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361
			plt.axhline(y=Tg_line,lw=0.5,color='k', linestyle='-.')

			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Entropy',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(TD_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

	if Isobars:

		if Plot_Phi:
			#Plotting the phi's of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,phip_P1_Condo,'r',ls=lsT1K,label='phi_p_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phip_P2_Condo,'r',ls=lsT2K,label='phi_p_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phip_P3_Condo,'r',ls=lsT3K,label='phi_p_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phis_P1_Condo,'m',ls=lsT1K,label='phi_s_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phis_P2_Condo,'m',ls=lsT2K,label='phi_s_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phis_P3_Condo,'m',ls=lsT3K,label='phi_s_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,Rtilde_P1_Condo,'b',ls=lsT1K,label='Rtilde_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Rtilde_P2_Condo,'b',ls=lsT2K,label='Rtilde_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Rtilde_P3_Condo,'b',ls=lsT3K,label='Rtilde_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phip0_P1_Condo,'k',ls=lsT1K,label='phi_p0_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phip0_P2_Condo,'k',ls=lsT2K,label='phi_p0_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phip0_P3_Condo,'k',ls=lsT3K,label='phi_p0_{} MPa'.format(P3),lw=linewidth)

			plt.plot(T0,phis0_P1_Condo,'y',ls=lsT1K,label='phi_s0_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,phis0_P2_Condo,'y',ls=lsT2K,label='phi_s0_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,phis0_P3_Condo,'y',ls=lsT3K,label='phi_s0_{} MPa'.format(P3),lw=linewidth)

			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('phi',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(phi_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Plot_Solubility:
			#Plotting the solubility of the PS+CO2 mixture.
			figX_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,Xs_P1_Condo,'k',ls=lsT1K,label='Present theory {} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Xs_P2_Condo,'k',ls=lsT2K,label='Present theory {} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Xs_P3_Condo,'k',ls=lsT3K,label='Present theory {} MPa'.format(P3),lw=linewidth)
			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('Gas solubility $\chi_s$ (g/g)',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			# plt.axis(solubility_axes)
			# figX.savefig('./'+output_folder+r'\bin_PS_CO2_Solubility'+img_extension,dpi=img_dpi)

		if Plot_Swelling:
			#Plotting the swelling of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
			plt.plot(T0,Sw_P1_Condo,'k',ls=lsT1K,label='Present theory {} K'.format(T1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,Sw_P2_Condo,'k',ls=lsT2K,label='Present theory {} K'.format(T2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,Sw_P3_Condo,'k',ls=lsT3K,label='Present theory {} K'.format(T3),lw=linewidth)
			plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
			plt.ylabel('Swelling ratio $S_W$ $(cm^3/cm^3)$',fontsize=axis_size)
			plt.legend(loc=2,fontsize=size,numpoints=1)
			plt.axis(swelling_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

		if Entropy:
			#Plotting the phi's of the PS+CO2 mixture.
			figS_Condo = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')

			plt.plot(T0,S_1_P1_Condo,'r',ls=lsT1K,label='S_1_{} MPa'.format(P1),lw=linewidth)
			if P2!=0.0:
				plt.plot(T0,S_1_P2_Condo,'r',ls=lsT2K,label='S_1_{} MPa'.format(P2),lw=linewidth)
			if P3!=0.0:
				plt.plot(T0,S_1_P3_Condo,'r',ls=lsT3K,label='S_1_{} MPa'.format(P3),lw=linewidth)

			if Plot_S2:
					
				plt.plot(T0,S_2_P1_Condo,'m',ls=lsT1K,label='S_2_{} MPa'.format(P1),lw=linewidth)
				if P2!=0.0:
					plt.plot(T0,S_2_P2_Condo,'m',ls=lsT2K,label='S_2_{} MPa'.format(P2),lw=linewidth)
				if P3!=0.0:
					plt.plot(T0,S_2_P3_Condo,'m',ls=lsT3K,label='S_2_{} MPa'.format(P3),lw=linewidth)
			
			# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
			plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')
			# S_max=npy.max(S_1_P1)
			# print 'S_max is:', S_max
			Tg_line=xS_infty#0.310707*0.8708171#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361
			plt.axhline(y=Tg_line,lw=0.5,color='k', linestyle='-.')

			plt.xlabel('Temperature T (K)',fontsize=axis_size)
			plt.ylabel('Entropy',fontsize=axis_size)
			plt.legend(loc=4,fontsize=size,numpoints=1)
			# plt.axis(TD_axes)
			# figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

#Show plot windows.
plt.show()
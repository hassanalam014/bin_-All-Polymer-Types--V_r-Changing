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

def discard_zeros(x,y):
	
	for i in range(len(x)):
		if x[i]==0:
			y[i]=0

	for i in range(len(y)):
		if y[i]==0:
			x[i]=0

	x = npy.delete(x, npy.argwhere( (x >= 0) & (x <= 0) ))
	y = npy.delete(y, npy.argwhere( (y >= 0) & (y <= 0) ))
				
	return x,y

Polymer_Type='PS'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}

Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
# P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

forward=True
backward=False

Kier=True
Hassan=False  
Condo=False 
Find_Tg_at_P=True
Find_Pg_at_T=True

Condo_Original=False 
CondoFind_Tg_at_P=False
CondoFind_Pg_at_T=False

# For zeta=1.124 and delta=0.97
# P=[0.101325,0.75673333,1.41214167,2.06755,2.72295833]
# Tg_bisect_fwd=[372.0703125,362.6953125,252.1484375,270.5078125,292.7734375]
# Pg_bisect_fwd=[1.33789062,1.92382812,2.48242188,2.84179688,2.84960938,2.42382812,1.58398438]
# T=[250.0,266.66666667,283.33333333,300.0,316.66666667,333.33333333,350.0]


# For zeta=1.1050 and delta=0.97
# P = [0.101325 ,1.75673333 ,3.41214167 ,5.06755 ,6.72295833 ,8.37836667 ,10.033775 ,11.68918333 ,13.34459167 ,15.0]
# Tg_bisect_fwd = [373.2421875 ,365.8203125 ,359.5703125 ,353.3203125 ,347.8515625 ,343.1640625 ,339.2578125 ,336.1328125 ,334.9609375 ,333.7890625]
# Pg_bisect_fwd = [1.54882812]
# T = [366.66666667]

# For zeta=1.100 and delta=0.97
# P1 = [2.84570312,3.54301184,4.11662632,4.69024079,5.26385526,5.83746974,6.41108421,6.98469868,7.55831316,8.13192763,8.70554211,9.27915658,9.85277105,10.42638553,11.0000000]								
# Tg1 = [265.78947368,265.4296875,266.6015625,265.8203125,265.0390625,264.2578125,263.4765625,262.6953125,261.9140625,261.1328125,260.7421875,259.9609375,259.1796875,258.3984375,257.6171875]								
# P = [0.101325,0.67493947,1.24855395,1.82216842,2.39578289,2.96939737]							
# Tg_bisect_fwd = [372.8515625,366.6015625,361.1328125,355.6640625,349.8046875,343.5546875]					
# Pg_bisect_fwd = [1.87304688,2.84570312,3.32421875,3.77734375,4.16796875,4.44140625,4.58203125,4.55859375,4.36328125,4.00390625,3.49609375,2.86523438,2.12695312,1.32226562]								
# T = [250.000000,265.78947368,273.68421053,281.57894737,289.47368421,297.36842105,305.26315789,313.15789474,321.05263158,328.94736842,336.84210526,344.73684211,352.63157895,360.52631579]								

# # For zeta=1.088 and delta=0.97 #Full Data
# P1 = [0.101325,3.42621458,6.75110417,10.07599375,13.40088333,16.72577292,20.0506625,23.37555208,26.70044167,30.02533125,33.35022083,36.67511042,40.00000]
# P2 = [0.101325,1.09211364,2.08290227,3.07369091,4.06447955,5.05526818,6.04605682,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
# P=P1+P2

# Tg_bisect_fwd1 = [372.8515625,345.8984375,294.7265625,291.6015625,288.8671875,286.1328125,283.7890625,281.4453125,279.4921875,277.5390625,275.1953125,273.2421875,271.2890625]
# Tg_bisect_fwd2 = [372.8515625,364.6484375,356.8359375,348.6328125,270.5078125,329.1015625,290.4296875,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	
# Tg_bisect_fwd=Tg_bisect_fwd1+Tg_bisect_fwd2

# Pg_bisect_fwd1 = [1.87304688,3.95703125,5.91015625,5.34765625,2.93359375]
# Pg_bisect_fwd2 = [2.79882812,5.53515625,5.72265625,3.46484375]		
# Pg_bisect_fwd=Pg_bisect_fwd1+Pg_bisect_fwd2

# T1 = [250.0,275.0,300.0,325.0,350.0]					
# T2 = [263.63636364,290.90909091,318.18181818,345.45454545]									
# T=T1+T2


# For zeta=1.088 and delta=0.97 #Data Upto 11 MPa without retrograde
# P = [0.101325,1.09211364,2.08290227,3.07369091,3.42621458,5.05526818,5.34765625,5.72265625,5.91015625,6.04605682,6.75110417,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
# Tg_bisect_fwd = [372.8515625,364.6484375,356.8359375,348.6328125,345.8984375,329.1015625,325.0,318.18181818,300.0,290.4296875,294.7265625,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	


# For zeta=1.088 and delta=0.97 #Data Upto 11 MPa with retrograde
P = [0.101325,1.09211364,2.08290227,3.07369091,3.42621458,5.05526818,5.34765625,5.72265625,5.91015625]
Tg_bisect_fwd = [372.8515625,364.6484375,356.8359375,348.6328125,345.8984375,329.1015625,325.0,318.18181818,300.0]

Pg_bisect_fwd = [1.87304688,2.79882812,3.95703125,5.53515625]
T = [250.0,263.63636364,275.0,290.90909091]

indep_list = Tg_bisect_fwd+T
dep_list = P+Pg_bisect_fwd
indep_list, dep_list = zip(*sorted(zip(indep_list, dep_list)))
indep_list=list(indep_list)
dep_list=list(dep_list)
indep_list=indep_list[::-1]
dep_list=dep_list[::-1]

print indep_list
print dep_list

P1 = [1.87304688,2.79882812,3.95703125,5.53515625,5.91015625,5.72265625,6.75110417,7.03684545,8.02763409,9.01842273,10.00921136,11.00000]	
Tg1 = [250.0,263.63636364,275.0,290.90909091,300.0,318.18181818,294.7265625,294.3359375,293.5546875,292.3828125,291.6015625,290.8203125]	

indep_list = indep_list+Tg1
dep_list = dep_list+P1


#Setting font size
axis_size = 20				#Size of x and y axis variables names
title_size = 20				#Size of title
size = 14					#Size of legender
label_size = 16				#Size of x and y axis variable values
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Linestyles
lsT1K = '-'
lsT2K = '--'
lsT3K = ':'


#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'Plots'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)


#General line properties.
linewidth = 1
markersize = 8


#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

# plt.plot(P_exp,Tg_exp,color='b',marker='o',ls='',label='Tg_exp_condo',ms=markersize)

# Condo et al.			
Tg_exp_condo = [373.18,328.40,308.21,305.84,304.49]
P_exp_condo = [0.100,6.072,7.139,8.419,10.440]

# Pham el al. 90 nm				
Tg_exp_pham90nm = [383.11,348.03,322.87,308.10,298.06]
P_exp_pham90nm = [0.089,3.990,4.792,5.209,4.884]

# Pham et al. 17 nm				
Tg_exp_pham17nm = [362.47,347.58,322.65,307.98,297.61]
P_exp_pham17nm = [0.089,1.278,3.218,4.071,3.675]

# Wissinger Bulk				
Tg_exp_wissinger = [373.86,338.67,324.00,308.89]	
P_exp_wissinger = [0.100,3.645,4.864,6.072]	

plt.plot(P_exp_condo,Tg_exp_condo,color='b',marker='o',ls='',label='Condo',ms=markersize)
plt.plot(P_exp_pham90nm,Tg_exp_pham90nm,color='r',marker='P',ls='',label='Pham',ms=markersize)
# plt.plot(P_exp_pham17nm,Tg_exp_pham17nm,color='b',marker='^',ls='',label='Pham et al. - h~17nm',ms=markersize)
# plt.plot(P_exp_wissinger,Tg_exp_wissinger,color='b',marker='s',ls='',label='Wissinger',ms=markersize)

# plt.plot(P_exp_Condo,Tg_exp_Condo,color='k',marker='o',ls='',label='Tg_exp_condo',ms=markersize)
plt.plot(dep_list,indep_list,color='k',marker='x',lw=linewidth,ls='-.',label='Present Theory')


# if Find_Tg_at_P:
# 	if forward:
# 		P,Tg_bisect_fwd=discard_zeros(P,Tg_bisect_fwd)
# 		plt.plot(P,Tg_bisect_fwd,color='k',marker='x',lw=linewidth,ls='-.',label='Tg_bisect_fwd')
# 	if backward:
# 		P,Tg_bisect_bwd=discard_zeros(P,Tg_bisect_bwd)
# 		plt.plot(P,Tg_bisect_bwd,color='b',lw=linewidth,ls='-',label='Tg_bisect_bwd')

# if Find_Pg_at_T:
# 	if forward:
# 		Pg_bisect_fwd,T=discard_zeros(Pg_bisect_fwd,T)
# 		plt.plot(Pg_bisect_fwd,T,color='k',marker='x',lw=linewidth,ls='-.',label='Pg_bisect_fwd')
# 	if backward:
# 		Pg_bisect_bwd,T=discard_zeros(Pg_bisect_bwd,T)
# 		plt.plot(Pg_bisect_bwd,T,color='b',lw=linewidth,ls='-',label='Pg_bisect_bwd')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature Tg (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=1,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()
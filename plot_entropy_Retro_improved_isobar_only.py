# Author: Hassan Alam
# Date: 2019
#
# Description: The purpose of this file is to plot Polystyrene (PS) Thermodynamics Properties
#
from __future__ import division
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from sympy import *
# from p_params import *
# from s_params import *
# from calculateBinaryVariablesCHV import *
from plot_phi_Retro_improved import *
from All_Functions import*

##################################################################################################
##################################################################################################

if Condo_Original or Condo:

	if Polymer_Type=='PS':
		cepsilon_s=0
		cepsilon_p=7151.0
		cz=5.0
	if Polymer_Type=='PMMA':
		cepsilon_s=0
		cepsilon_p=7443.0
		cz=5.0

if Kier or Hassan:
		
	if Polymer_Type=='PS':
		kg=1.67 
		kepsilon_p=8013 
		kx=0.311

		hg=1.67 
		hepsilon_p=8013 
		hx=0.311

	if Polymer_Type=='PMMA':	
		kg=1.66 
		kepsilon_p=8094 
		kx=0.323

		hg=1.66 
		hepsilon_p=8094 
		hx=0.323

# phi_p=npy.zeros(len(T))
# phi_s=npy.zeros(len(T))
# phi_s0=npy.zeros(len(T))
# Fp=npy.zeros(len(T))
# Fs=npy.zeros(len(T))
# S_term1=npy.zeros(len(T))		
# S_term2=npy.zeros(len(T))		
# S_term3=npy.zeros(len(T))		
# S_term4=npy.zeros(len(T))
# S_term5=npy.zeros(len(T))
# S_term6=npy.zeros(len(T))
# S_term7=npy.zeros(len(T))

cS_1_bisect_fwd=npy.zeros(len(T))
cS_1_bisect_bwd=npy.zeros(len(T))
cS_2_bisect_fwd=npy.zeros(len(T))
cS_2_bisect_bwd=npy.zeros(len(T))

S_1_bisect_fwd=npy.zeros(len(T))
S_1_bisect_bwd=npy.zeros(len(T))
S_2_bisect_fwd=npy.zeros(len(T))
S_2_bisect_bwd=npy.zeros(len(T))

if Kier or Condo or Hassan:		
	if Kier:
		kwargs = {'g':kg,'epsilon_p':kepsilon_p,'Ppstar':kPpstar,'Tpstar':kTpstar,'Rpstar':kRpstar,'Psstar':kPsstar,'Tsstar':kTsstar,'Rsstar':kRsstar,'zeta':kzeta,'delta':kdelta}
	if Hassan:	
		kwargs = {'g':hg,'epsilon_p':hepsilon_p,'Ppstar':hPpstar,'Tpstar':hTpstar,'Rpstar':hRpstar,'Psstar':hPsstar,'Tsstar':hTsstar,'Rsstar':hRsstar,'zeta':hzeta,'delta':hdelta}
	if Condo:
		kwargs = {'z':cz,'epsilon_p':cepsilon_p,'epsilon_s':cepsilon_s,'Ppstar':cPpstar,'Tpstar':cTpstar,'Rpstar':cRpstar,'Psstar':cPsstar,'Tsstar':cTsstar,'Rsstar':cRsstar,'zeta':czeta}

	if forward:
		phi_s,phi_p,Rtilde,Rstilde,Rptilde=phi_sm_bisect_fwd,phi_pm_bisect_fwd,sum_phi_pm_sm_bisect_fwd,phi_s0_bisect_fwd,phi_p0_bisect_fwd

		if Kier:
			for i in range(0,len(T)):	
				S_1_bisect_fwd[i],S_2_bisect_fwd[i]= ThermodynamicVariables(P,T[i],phi_p[i],phi_s[i],Rptilde[i],Mp,Ms,**kwargs)
		if Hassan:	
			for i in range(0,len(T)):	
				S_1_bisect_fwd[i],S_2_bisect_fwd[i] = ThermodynamicVariables(P,T[i],phi_p[i],phi_s[i],Rptilde[i],Mp,Ms,**kwargs)
		if Condo:
			for i in range(0,len(T)):	
				S_1_bisect_fwd[i],S_2_bisect_fwd[i]= CondoThermodynamicVariables_Modified(P,T[i],phi_p[i],phi_s[i],Mp,Ms,**kwargs)

	if backward:
		phi_s,phi_p,Rtilde,Rstilde,Rptilde=phi_sm_bisect_bwd,phi_pm_bisect_bwd,sum_phi_pm_sm_bisect_bwd,phi_s0_bisect_bwd,phi_p0_bisect_bwd

		if Kier:
			for i in range(0,len(T)):	
				S_1_bisect_bwd[i],S_2_bisect_bwd[i]= ThermodynamicVariables(P,T[i],phi_p[i],phi_s[i],Rptilde[i],Mp,Ms,**kwargs)
		if Hassan:	
			for i in range(0,len(T)):	
				S_1_bisect_bwd[i],S_2_bisect_bwd[i] = ThermodynamicVariables(P,T[i],phi_p[i],phi_s[i],Rptilde[i],Mp,Ms,**kwargs)
		if Condo:
			for i in range(0,len(T)):	
				S_1_bisect_bwd[i],S_2_bisect_bwd[i]= CondoThermodynamicVariables_Modified(P,T[i],phi_p[i],phi_s[i],Mp,Ms,**kwargs)

if Condo_Original:
	
	kwargs = {'z':cz,'epsilon_p':cepsilon_p,'epsilon_s':cepsilon_s,'Ppstar':cPpstar,'Tpstar':cTpstar,'Rpstar':cRpstar,'Psstar':cPsstar,'Tsstar':cTsstar,'Rsstar':cRsstar,'zeta':czeta}

	if forward:
		phi_s,phi_p,Rtilde,Rstilde=cphi_sm_bisect_fwd,cphi_pm_bisect_fwd,Rtilde_bisect_fwd,Rstilde_bisect_fwd

		for i in range(0,len(T)):	
			cS_1_bisect_fwd[i],cS_2_bisect_fwd[i] = CondoThermodynamicVariables_Original(P,T[i],Rtilde[i],phi_s[i],Mp,Ms,**kwargs)
		
	if backward:
		phi_s,phi_p,Rtilde,Rstilde=cphi_sm_bisect_bwd,cphi_pm_bisect_bwd,Rtilde_bisect_bwd,Rstilde_bisect_bwd

		for i in range(0,len(T)):	
			cS_1_bisect_bwd[i],cS_2_bisect_bwd[i] = CondoThermodynamicVariables_Original(P,T[i],Rtilde[i],phi_s[i],Mp,Ms,**kwargs)

#Setting font size
axis_size = 20
title_size = 20
size = 14
label_size = 20
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Setting saved image properties
img_extension = '.png'
img_dpi = None
output_folder = 'Plots'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#General line properties.
linewidth = 1
markersize = 6

arrow_ls = 'dashdot'
show_arrows = True

#==================================================================================
#Plots.
figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
ax = plt.axes()

if Condo_Original:
	if forward:
		# plt.plot(T,cphi_pm_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='cphi_pm_bisect_fwd')
		# plt.plot(T,cphi_sm_bisect_fwd,'k',color='m',lw=linewidth,ls='-',label='cphi_sm_bisect_fwd')
		# plt.plot(T,Rstilde_bisect_fwd,'k',color='b',lw=linewidth,ls='-',label='Rstilde_bisect_fwd')
		# plt.plot(T,Rtilde_bisect_fwd,'k',color='y',lw=linewidth,ls='-',label='Rtilde_bisect_fwd')

		# plt.plot(T,kphi_pm_bisect_fwd,'k',color='r',lw=linewidth,ls='-',label='kphi_pm_bisect_fwd')
		# plt.plot(T,kphi_sm_bisect_fwd,'k',color='k',lw=linewidth,ls='-',label='kphi_sm_bisect_fwd')

		plt.plot(T,cS_1_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='cS_1_bisect_fwd')
		plt.plot(T,cS_2_bisect_fwd,'k',color='m',lw=linewidth,ls='-',label='cS_2_bisect_fwd')

	if backward:
		# plt.plot(T,cphi_pm_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='cphi_pm_bisect_bwd')
		# plt.plot(T,cphi_sm_bisect_bwd,'k',color='m',lw=linewidth,ls='-.',label='cphi_sm_bisect_bwd')
		# plt.plot(T,Rstilde_bisect_bwd,'k',color='b',lw=linewidth,ls='-.',label='Rstilde_bisect_bwd')
		# plt.plot(T,Rtilde_bisect_bwd,'k',color='y',lw=linewidth,ls='-.',label='Rtilde_bisect_bwd')

		# plt.plot(T,kphi_pm_bisect_bwd,'k',color='r',lw=linewidth,ls='-.',label='kphi_pm_bisect_bwd')
		# plt.plot(T,kphi_sm_bisect_bwd,'k',color='k',lw=linewidth,ls='-.',label='kphi_sm_bisect_bwd')

		plt.plot(T,cS_1_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='cS_1_bisect_bwd')
		plt.plot(T,cS_2_bisect_bwd,'k',color='m',lw=linewidth,ls='-.',label='cS_2_bisect_bwd')

if Kier or Hassan or Condo:
	if forward:
		# plt.plot(T,phi_pm_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='phi_pm_bisect_fwd')
		# plt.plot(T,phi_sm_bisect_fwd,'k',color='m',lw=linewidth,ls='-',label='phi_sm_bisect_fwd')
		# plt.plot(T,phi_s0_bisect_fwd,'k',color='b',lw=linewidth,ls='-',label='phi_s0_bisect_fwd')
		# plt.plot(T,sum_phi_pm_sm_bisect_fwd,'k',color='y',lw=linewidth,ls='-',label='sum_phi_pm_sm_bisect_fwd')

		plt.plot(T,S_1_bisect_fwd,'k',color='k',lw=linewidth,ls='-',label='S_1_bisect_fwd')
		plt.plot(T,S_2_bisect_fwd,'k',color='b',lw=linewidth,ls='-',label='S_2_bisect_fwd')

	if backward:
		# plt.plot(T,phi_pm_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='phi_pm_bisect_bwd')
		# plt.plot(T,phi_sm_bisect_bwd,'k',color='m',lw=linewidth,ls='-.',label='phi_sm_bisect_bwd')
		# plt.plot(T,phi_s0_bisect_bwd,'k',color='b',lw=linewidth,ls='-.',label='phi_s0_bisect_bwd')
		# plt.plot(T,sum_phi_pm_sm_bisect_bwd,'k',color='y',lw=linewidth,ls='-.',label='sum_phi_pm_sm_bisect_bwd')

		plt.plot(T,S_1_bisect_bwd,'k',color='k',lw=linewidth,ls='-',label='S_1_bisect_bwd')
		plt.plot(T,S_2_bisect_bwd,'k',color='b',lw=linewidth,ls='-',label='S_2_bisect_bwd')

# plt.plot(T,Fp,'k',color='r',lw=linewidth,ls='-',label='Fp')
# plt.plot(T,Fs,'k',color='g',lw=linewidth,ls='-',label='Fs')

# plt.plot(T,R,'k',color='y',lw=linewidth,ls='-',label='R')
# plt.plot(T,E,'k',color='g',lw=linewidth,ls='-',label='E')
# plt.plot(T,H,'k',color='m',lw=linewidth,ls='-',label='H')

# plt.plot(T,phi_p,'k',color='g',lw=linewidth,ls='-',label='phi_p')
# plt.plot(T,phi_s,'k',color='m',lw=linewidth,ls='-',label='phi_s')

# plt.plot(T,S,'k',color='k',lw=linewidth,ls='-',label='S')
# plt.plot(T,S_term1,'k',color='b',lw=linewidth,ls='-',label='S_term1')
# plt.plot(T,S_term2,'k',color='g',lw=linewidth,ls='-',label='S_term2')
# plt.plot(T,S_term3,'k',color='r',lw=linewidth,ls='-',label='S_term3')
# plt.plot(T,S_term4,'k',color='m',lw=linewidth,ls='-',label='S_term4')
# plt.plot(T,S_term5,'k',color='y',lw=linewidth,ls='-',label='S_term5')

# plt.plot(T,Cp,'k',color='c',lw=linewidth,ls='-',label='Cp')
# plt.plot(T,Cv,'k',color='r',lw=linewidth,ls='-',label='Cv')

# plt.axvline(x=378,lw=0.5,color='k', linestyle='-.')
plt.axhline(y=0.0,lw=0.5,color='k', linestyle='-.')

S_max=npy.max(S_1_bisect_fwd)
print 'S_max is:', S_max
Tg_line=0.323*1.1266666	#0.271*0.8708#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361				#0.310707*0.8708171
plt.axhline(y=Tg_line,lw=0.5,color='k', linestyle='-.')

Tg_line_2=0.323*1.1466666	#0.310707*0.90#0.2361#2.38*0.271##0.2558 #S_max*kx #0.2361				#0.310707*0.8708171
plt.axhline(y=Tg_line_2,lw=0.5,color='k', linestyle='-.')


plt.xlabel('Temperature T (K)',fontsize=axis_size)
plt.ylabel(r'TD Property',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()


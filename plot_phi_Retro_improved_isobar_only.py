import os,math,numpy as npy
from loadPhysicalConstants import *
from checkResults import *
from sympy import *
import warnings
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
import os,sys,math,matplotlib.pyplot as plt,numpy as npy
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
# from p_params import *
# from s_params import *
import cmath
# import time
from Parameters_of_Different_Polymers import *
from All_Functions import *
from Parameters_for_Mixtures_and_Tg import *

#####################################################################################################
#####################################################################################################

# X=binarySolubilitySwellingCHV(P,T[i],phis,phip,Mp,Ms,**kwargs)

Polymer_Type='PMMA'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}

Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

# equation='Kier'	#	'Kier', 'Condo', 'Hassan'
Calculate_Phi=True
Plot_Phi=True
Plot_Solubility=True

Kier=True
Hassan=False
Condo=False
Condo_Original=True

forward=True
backward=False 
direction=['fwd','bwd']

kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Kier':Kier,'Hassan':Hassan,'Condo':Condo,'Condo_Original':True}
cepsilon_s,cepsilon_p,cz,czeta,epsilon_p,g,x,delta,zeta=Parameters_for_Mixtures_and_Tg(**kwargs)

''''
if Kier:

	kPpstar,kTpstar,kRpstar,Mp,kPsstar,kTsstar,kRsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)

	if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Kier':
		kzeta=1.02127006
		kdelta=0.88102354

	if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Condo':
		kzeta=1.08820786	
		kdelta=0.97423316

	if Polymer_Type=='PMMA' and Solvent=='CO2' and Parameters_Paper=='Condo':
		kzeta=1.135		#1.12272978	#1.07643522	#	#1.1684#1.135			#1.10752004		#1.08621732		#1.135		
		kdelta=1.00		#0.73223348	#1.0		#	#0.5#0.74264552			#0.90758208		#1.05473203		#1.00	

if Hassan:

	hPpstar,hTpstar,hRpstar,Mp,hPsstar,hTsstar,hRsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)

	if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Kier':
		hzeta=0.92161759		#Hassan Values
		hdelta=0.88399214

	if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Condo':
		hzeta=0.98508210	
		hdelta=0.80548596

	if Polymer_Type=='PMMA' and Solvent=='CO2' and Parameters_Paper=='Condo':
		kzeta=1.02948288
		kdelta=0.90455520	

if Condo or Condo_Original:

	cPpstar,cTpstar,cRpstar,Mp,cPsstar,cTsstar,cRsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)

	if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Condo':
		#Condo Zeta:
		czeta=1.110			#Condo Values. Correlation value: 1.110, iteration value: 1.1240
		cdelta=100.00000	#Random value will never be used by the code

	if Polymer_Type=='PMMA' and Solvent=='CO2' and Parameters_Paper=='Condo':
		#Condo Zeta:
		czeta=1.1350
		cdelta=100.00000	#Random value will never be used by the code
'''

if Calculate_Phi:

	P=4.0#1.101325#7.39

	#Initializing the array of densities.
	T=npy.linspace(200,400,10)				#Tmax=1344479

	for i in range(0,len(direction)):

		exec "phi_pm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "phi_sm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "phi_s0_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "phi_p0_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "sum_phi_pm_sm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "X_bisect_%s=npy.zeros(len(T))" %(direction[i])


		#All Below Are Condo_Original:
		exec "cphi_pm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "cphi_sm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "Rtilde_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "Rstilde_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "Rptilde_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "kphi_pm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "kphi_sm_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "cX_bisect_%s=npy.zeros(len(T))" %(direction[i])
		exec "kX_bisect_%s=npy.zeros(len(T))" %(direction[i])

	if Condo_Original:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':czeta}

	if Condo_Original:
		if forward:
			
			for j in range(0,len(T)):
				print 'Iterating for T:', T[j], 'in forward direction for bisect method'
				cphi_sm_bisect_fwd[j],cphi_pm_bisect_fwd[j],Rtilde_bisect_fwd[j],Rstilde_bisect_fwd[j],Rptilde_bisect_fwd[j]=binaryPhaseEquilibriumCondo_Original_bisect('fwd',P,T[j],Mp,Ms,**kwargs)
				kphi_sm_bisect_fwd[j]=cphi_sm_bisect_fwd[j]*Rtilde_bisect_fwd[j]
				kphi_pm_bisect_fwd[j]=cphi_pm_bisect_fwd[j]*Rtilde_bisect_fwd[j]
				cX_bisect_fwd[j]=Solubility_Condo(P,T[i],cphi_sm_bisect_fwd[j],cphi_pm_bisect_fwd[j],Mp,Ms,**kwargs)
				kX_bisect_fwd[j]=Solubility_Condo(P,T[i],kphi_sm_bisect_fwd[j],kphi_pm_bisect_fwd[j],Mp,Ms,**kwargs)

		if backward:
				
			for j in range(0,len(T)):
				print 'Iterating for T:', T[j], 'in backward direction for bisect method'
				cphi_sm_bisect_bwd[j],cphi_pm_bisect_bwd[j],Rtilde_bisect_bwd[j],Rstilde_bisect_bwd[j],Rptilde_bisect_bwd[j]=binaryPhaseEquilibriumCondo_Original_bisect('bwd',P,T[j],Mp,Ms,**kwargs)
				kphi_sm_bisect_bwd[j]=cphi_sm_bisect_bwd[j]*Rtilde_bisect_bwd[j]
				kphi_pm_bisect_bwd[j]=cphi_pm_bisect_bwd[j]*Rtilde_bisect_bwd[j]
				cX_bisect_bwd[j]=Solubility_Condo(P,T[i],cphi_sm_bisect_bwd[j],cphi_pm_bisect_bwd[j],Mp,Ms,**kwargs)
				kX_bisect_bwd[j]=Solubility_Condo(P,T[i],kphi_sm_bisect_bwd[j],kphi_pm_bisect_bwd[j],Mp,Ms,**kwargs)

	if Kier or Hassan:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta,'method':'disparate','Kier':Kier,'Hassan':Hassan,'Condo':Condo}
	if Condo:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':czeta,'method':'disparate','Kier':Kier,'Hassan':Hassan,'Condo':Condo}

	if Kier or Hassan or Condo:
		if forward:
			
			for j in range(0,len(T)):
				print 'Iterating for T:', T[j], 'in forward direction for bisect method'
				_,_,phi_pm_bisect_fwd[j],phi_sm_bisect_fwd[j],phi_p0_bisect_fwd[j],phi_s0_bisect_fwd[j]=binaryPhaseEquilibriumCHV_bisect('fwd',P,T[j],Mp,Ms,**kwargs)
				phi_s0_bisect_fwd[j]=phi_pure_k_bisect('fwd',P,T[j],Ms,Psstar,Tsstar,Rsstar)
				sum_phi_pm_sm_bisect_fwd[j]=phi_pm_bisect_fwd[j]+phi_sm_bisect_fwd[j]
				X_bisect_fwd[j]=Solubility(P,T[i],phi_sm_bisect_fwd[j],phi_pm_bisect_fwd[j],Mp,Ms,**kwargs)

		if backward:
			
			for j in range(0,len(T)):
				print 'Iterating for T:', T[j], 'in backward direction for bisect method'
				_,_,phi_pm_bisect_bwd[j],phi_sm_bisect_bwd[j],phi_p0_bisect_bwd[j],phi_s0_bisect_bwd[j]=binaryPhaseEquilibriumCHV_bisect('bwd',P,T[j],Mp,Ms,**kwargs)
				phi_s0_bisect_bwd[j]=phi_pure_k_bisect('bwd',P,T[j],Ms,Psstar,Tsstar,Rsstar)
				sum_phi_pm_sm_bisect_bwd[j]=phi_pm_bisect_bwd[j]+phi_sm_bisect_bwd[j]
				X_bisect_bwd[j]=Solubility(P,T[i],phi_sm_bisect_bwd[j],phi_pm_bisect_bwd[j],Mp,Ms,**kwargs)

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
if Plot_Phi:

	figPUREPS=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
	ax = plt.axes()

	if Condo_Original:

		if forward:
			# plt.plot(T,cphi_pm_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='cphi_pm_bisect_fwd')
			# plt.plot(T,cphi_sm_bisect_fwd,'k',color='m',lw=linewidth,ls='-',label='cphi_sm_bisect_fwd')
			plt.plot(T,Rstilde_bisect_fwd,'k',color='b',lw=linewidth,ls='-',label='Rstilde_bisect_fwd')
			plt.plot(T,Rptilde_bisect_fwd,'k',color='c',lw=linewidth,ls='-',label='Rptilde_bisect_fwd')
			plt.plot(T,Rtilde_bisect_fwd,'k',color='y',lw=linewidth,ls='-',label='Rtilde_bisect_fwd')

			plt.plot(T,kphi_pm_bisect_fwd,'k',color='r',lw=linewidth,ls='-',label='kphi_pm_bisect_fwd')
			plt.plot(T,kphi_sm_bisect_fwd,'k',color='k',lw=linewidth,ls='-',label='kphi_sm_bisect_fwd')

		if backward:
			# plt.plot(T,cphi_pm_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='cphi_pm_bisect_bwd')
			# plt.plot(T,cphi_sm_bisect_bwd,'k',color='m',lw=linewidth,ls='-.',label='cphi_sm_bisect_bwd')
			plt.plot(T,Rstilde_bisect_bwd,'k',color='b',lw=linewidth,ls='-.',label='Rstilde_bisect_bwd')
			plt.plot(T,Rptilde_bisect_bwd,'k',color='c',lw=linewidth,ls='-.',label='Rptilde_bisect_bwd')
			plt.plot(T,Rtilde_bisect_bwd,'k',color='y',lw=linewidth,ls='-.',label='Rtilde_bisect_bwd')

			plt.plot(T,kphi_pm_bisect_bwd,'k',color='r',lw=linewidth,ls='-.',label='kphi_pm_bisect_bwd')
			plt.plot(T,kphi_sm_bisect_bwd,'k',color='k',lw=linewidth,ls='-.',label='kphi_sm_bisect_bwd')

		# print 'Rtilde_bisect_fwd=',Rtilde_bisect_fwd
		# print 'Rtilde_bisect_bwd=',Rtilde_bisect_bwd
		# print 'Rptilde_bisect_fwd=',Rptilde_bisect_fwd
		# print 'Rptilde_bisect_bwd=',Rptilde_bisect_bwd
		# print 'Rstilde_bisect_fwd=',Rstilde_bisect_fwd
		# print 'Rstilde_bisect_bwd=',Rstilde_bisect_bwd
		# print 'cphi_pm_bisect_fwd=',cphi_pm_bisect_fwd
		# print 'cphi_pm_bisect_bwd=',cphi_pm_bisect_bwd
		# print 'cphi_sm_bisect_fwd=',cphi_sm_bisect_fwd
		# print 'cphi_sm_bisect_bwd=',cphi_sm_bisect_bwd
		# print 'kphi_pm_bisect_fwd=',kphi_pm_bisect_fwd
		# print 'kphi_pm_bisect_bwd=',kphi_pm_bisect_bwd
		# print 'kphi_sm_bisect_fwd=',kphi_sm_bisect_fwd
		# print 'kphi_sm_bisect_bwd=',kphi_sm_bisect_bwd

	if Kier or Hassan or Condo:

		if forward:
			plt.plot(T,phi_pm_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='phi_pm_bisect_fwd')
			plt.plot(T,phi_sm_bisect_fwd,'k',color='m',lw=linewidth,ls='-',label='phi_sm_bisect_fwd')
			plt.plot(T,phi_s0_bisect_fwd,'k',color='b',lw=linewidth,ls='-',label='phi_s0_bisect_fwd')
			plt.plot(T,phi_p0_bisect_fwd,'k',color='c',lw=linewidth,ls='-',label='phi_p0_bisect_fwd')
			plt.plot(T,sum_phi_pm_sm_bisect_fwd,'k',color='y',lw=linewidth,ls='-',label='sum_phi_pm_sm_bisect_fwd')

		if backward:
			plt.plot(T,phi_pm_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='phi_pm_bisect_bwd')
			plt.plot(T,phi_sm_bisect_bwd,'k',color='m',lw=linewidth,ls='-.',label='phi_sm_bisect_bwd')
			plt.plot(T,phi_s0_bisect_bwd,'k',color='b',lw=linewidth,ls='-.',label='phi_s0_bisect_bwd')
			plt.plot(T,phi_p0_bisect_bwd,'k',color='c',lw=linewidth,ls='-.',label='phi_p0_bisect_bwd')
			plt.plot(T,sum_phi_pm_sm_bisect_bwd,'k',color='y',lw=linewidth,ls='-.',label='sum_phi_pm_sm_bisect_bwd')

	plt.xlabel('Temperature T (K)',fontsize=axis_size)
	plt.ylabel(r'phi',fontsize=axis_size)
	#plt.axis([300,500,0,1.5])
	plt.legend(loc=1,fontsize=size,numpoints=1)
	plt.subplots_adjust(bottom=0.3)

	# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

if Plot_Solubility:

	figX=plt.figure(num=None, figsize=(10,6), dpi=img_dpi, facecolor='w', edgecolor='k')
	ax = plt.axes()

	if Condo_Original:
		if forward:
			plt.plot(T,cX_bisect_fwd,'k',color='r',lw=linewidth,ls='-',label='cX_bisect_fwd')
			plt.plot(T,kX_bisect_fwd,'k',color='k',lw=linewidth,ls='-',label='kX_bisect_fwd')

		if backward:
			plt.plot(T,cX_bisect_bwd,'k',color='r',lw=linewidth,ls='-.',label='cX_bisect_bwd')
			plt.plot(T,kX_bisect_bwd,'k',color='k',lw=linewidth,ls='-.',label='kX_bisect_bwd')

	if Kier or Hassan or Condo:
		if forward:
			plt.plot(T,X_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='X_bisect_fwd')

		if backward:
			plt.plot(T,X_bisect_bwd,'k',color='g',lw=linewidth,ls='-.',label='X_bisect_bwd')

	plt.xlabel('Temperature T (K)',fontsize=axis_size)
	plt.ylabel(r'X(g/g)',fontsize=axis_size)
	#plt.axis([300,500,0,1.5])
	plt.legend(loc=1,fontsize=size,numpoints=1)
	plt.subplots_adjust(bottom=0.3)

	# figX.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()

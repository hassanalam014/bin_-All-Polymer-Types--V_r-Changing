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
# from plot_entropy_Retro_improved import *
from loadExperimentalData import *
from Parameters_of_Different_Polymers import *
from All_Functions import *
from Parameters_for_Mixtures_and_Tg import *
from To_get_colored_print import *

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

#Kier or Hassan:

def entropy_infty(P_random,Mp,Ms,**kwargs):

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	#Initializing the array of densities.
	T=npy.linspace(300,400,2)
	phip=npy.zeros(len(T))
	phis=npy.zeros(len(T))	
	phip0=npy.zeros(len(T))
	phis0=npy.zeros(len(T))
	Rstilde=npy.zeros(len(T))
	Rptilde=npy.zeros(len(T))
	Rtilde=npy.zeros(len(T))
	S_1=npy.zeros(len(T))
	S_2=npy.zeros(len(T))
	Xs=npy.zeros(len(T))
	Sw=npy.zeros(len(T))

	for i in range(0,len(T)):	

		result = calculateBinarySolubilitySwelling('CHV',P_random,T[i],Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol,forward=forward,backward=backward)
		Xs[i] = result[2]
		Sw[i] = result[3]
		phip[i] = result[4]
		phis[i] = result[5]
		Rtilde[i] = result[6]
		phip0[i] = result[7]
		phis0[i] = result[8]

		properties=calculateThermodynamicVariables(P_random,T[i],phip[i],phis[i],phip0[i],phis0[i],Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
		print properties
		S_1[i] = properties[0]
		S_2[i] = properties[1]

	S=S_1

	S_max=npy.max(S)
	S_infty=S[-1]
	print S_infty
	# S_infty=0.8708171	
	# S_infty=0.98268		#For PC Pure
	alpha_0=vh_pure_p/v_r
	Vratio=vh_pure_p/v_r
	S_infty=(Ppstar/(Rpstar*Tpstar))*((1/alpha_0)*(1+((alpha_0/Vratio)*(ln(1+g)))))
	return S_max,S_infty

def Find_Tg_Bisect_xS_infty(Tg,P,S_infty,Mp,Ms,g,epsilon_p,x,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol,forward,backward):

	result = calculateBinarySolubilitySwelling('CHV',P,Tg,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol,forward=forward,backward=backward)
	Xs = result[2]
	Sw = result[3]
	phip = result[4]
	phis = result[5]
	Rtilde = result[6]
	phip0 = result[7]
	phis0 = result[8]

	properties=calculateThermodynamicVariables(P,Tg,phip,phis,phip0,phis0,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
	print properties
	S_1 = properties[0]
	S_2 = properties[1]

	criterion=S_1-x*S_infty

	return criterion

def Find_Pg_Bisect_xS_infty(Pg,T,S_infty,Mp,Ms,g,epsilon_p,x,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol,forward,backward):

	result = calculateBinarySolubilitySwelling('CHV',Pg,T,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r,Kier=Kier,Hassan=Hassan,Condo=Condo,Hassan_Var_Vol=Hassan_Var_Vol,forward=forward,backward=backward)
	Xs = result[2]
	Sw = result[3]
	phip = result[4]
	phis = result[5]
	Rtilde = result[6]
	phip0 = result[7]
	phis0 = result[8]

	properties=calculateThermodynamicVariables(Pg,T,phip,phis,phip0,phis0,Mp,Ms,g=g,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,vh_pure_p=vh_pure_p,vh_pure_s=vh_pure_s,v_r=v_r)
	S_1 = properties[0]
	S_2 = properties[1]
	
	criterion=S_1-x*S_infty

	return criterion

def GlassTemperature(direction,P,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	max_Tg=450
	min_Tg=250
	step_Tg=25

	if direction=='fwd':
		start=min_Tg
		end=max_Tg
		step=step_Tg
		# print 'forward'
		
	elif direction=='bwd':
		start=max_Tg
		end=min_Tg
		step=-1*step_Tg
		# print 'backward'

	for i in range(start,end,step):
		Tg=0.0
		try:
			Tg = bisect(Find_Tg_Bisect_xS_infty,i,i+step,args=(P,S_infty,Mp,Ms,g,epsilon_p,x,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol,forward,backward),xtol=1e-3, rtol=1e-3)
		except:
			# print 'No value found'
			pass
		if Tg!=0.0:
			prRed('Hurry! Tg is:{} for direction {}'.format(Tg,direction))
			break
	if Tg==0.0:
		print 'Program Failed to get value of Tg in given bisect range in direction', direction

	return Tg

def GlassPressure(direction,T,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	max_Pg=7
	min_Pg=1
	step_Pg=1
	if direction=='fwd':
		start=min_Pg
		end=max_Pg
		step=step_Pg
		# print 'forward'
		
	elif direction=='bwd':
		start=max_Pg
		end=min_Pg
		step=-1*step_Pg
		# print 'backward'

	for i in range(start,end,step):
		Pg=0.0
		try:
			Pg = bisect(Find_Pg_Bisect_xS_infty,i,i+step,args=(T,S_infty,Mp,Ms,g,epsilon_p,x,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol,forward,backward),xtol=1e-3, rtol=1e-3)
		except:
			# print 'No value found'
			pass
		if Pg!=0.0:
			prRed('Hurry! Pg is:{} for direction {}'.format(Pg,direction))
			break
	if Pg==0.0:
		print 'Program Failed to get value of Pg in given bisect range in direction', direction

	return Pg

#Condo Original:

def Find_Tg_Bisect_CondoEntropy_Original(Tg,P,Mp,Ms,z,epsilon_s,epsilon_p,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,forward,backward):

	result = calculateBinarySolubilitySwelling('Condo_Original',P,Tg,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
	Xs = result[2]
	Sw = result[3]
	phip = result[4]
	phis = result[5]
	Rtilde = result[6]
	phip0 = result[7]
	phis0 = result[8]
	
	properties=calculateCondoThermodynamicVariables_Original(P,Tg,phip,phis,phip0,phis0,Mp,Ms,z=z,epsilon_s=epsilon_s,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
	S_1 = properties[0]
	S_2 = properties[1]

	return S_1

def Find_Pg_Bisect_CondoEntropy_Original(Pg,T,Mp,Ms,z,epsilon_s,epsilon_p,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,forward,backward):

	result = calculateBinarySolubilitySwelling('Condo_Original',Pg,T,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,forward=forward,backward=backward)
	Xs = result[2]
	Sw = result[3]
	phip = result[4]
	phis = result[5]
	Rtilde = result[6]
	phip0 = result[7]
	phis0 = result[8]
	
	properties=calculateCondoThermodynamicVariables_Original(Pg,T,phip,phis,phip0,phis0,Mp,Ms,z=z,epsilon_s=epsilon_s,epsilon_p=epsilon_p,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar)
	S_1 = properties[0]
	S_2 = properties[1]

	return S_1

def CondoGlassTemperature_Original(direction,P,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	max_Tg=600
	min_Tg=300
	step_Tg=25

	if direction=='fwd':
		start=min_Tg
		end=max_Tg
		step=step_Tg
		# print 'forward'
		
	elif direction=='bwd':
		start=max_Tg
		end=min_Tg
		step=-1*step_Tg
		# print 'backward'

	for i in range(start,end,step):
		Tg=0.0
		try:
			Tg = bisect(Find_Tg_Bisect_CondoEntropy_Original,i,i+step,args=(P,Mp,Ms,z,epsilon_s,epsilon_p,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,forward,backward),xtol=1e-3, rtol=1e-3)
		except:
			# print 'No value found'
			pass
		if Tg!=0.0:
			prRed('Hurry! Tg is:{} for direction {}'.format(Tg,direction))
			break
	if Tg==0.0:
		print 'Program Failed to get value of Tg in given bisect range in direction', direction

	return Tg

def CondoGlassPressure_Original(direction,T,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	max_Pg=13
	min_Pg=1
	step_Pg=1

	if direction=='fwd':
		start=min_Pg
		end=max_Pg
		step=step_Pg
		# print 'forward'
		
	elif direction=='bwd':
		start=max_Pg
		end=min_Pg
		step=-1*step_Pg
		# print 'backward'

	for i in range(start,end,step):
		Pg=0.0
		try:
			Pg = bisect(Find_Pg_Bisect_CondoEntropy_Original,i,i+step,args=(T,Mp,Ms,z,epsilon_s,epsilon_p,zeta,delta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,forward,backward),xtol=1e-3, rtol=1e-3)
		except:
			# print 'No value found'
			pass
		if Pg!=0.0:
			prRed('Hurry! Pg is:{} for direction {}'.format(Pg,direction))
			break
	if Pg==0.0:
		print 'Program Failed to get value of Pg in given bisect range in direction', direction

	return Pg

####################################################################################################

Polymer_Type='PC'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}

Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
# P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

forward=True
backward=False

Kier=False
Hassan=False  
Hassan_Var_Vol=True  
Condo=False 
Find_Tg_at_P=True
Find_Pg_at_T=True

Condo_Original=False 
CondoFind_Tg_at_P=False
CondoFind_Pg_at_T=False

Rescaling=True
if Rescaling:
	vh_pure_p=kB*Tpstar/Ppstar
	vh_pure_s=kB*Tsstar/Psstar
	v_r=vh_pure_s
	Tpstar,Tsstar=rescaling(Tpstar,Tsstar,vh_pure_p,vh_pure_s,v_r)

P=npy.linspace(0.101325,13,7)
T=npy.linspace(200,500,7)

kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Kier':Kier,'Hassan':Hassan,'Hassan_Var_Vol':Hassan_Var_Vol,'Condo':Condo,'Condo_Original':True}
cepsilon_s,cepsilon_p,cz,czeta,epsilon_p,g,x,delta,zeta=Parameters_for_Mixtures_and_Tg(**kwargs)
cdelta=100.0

zeta=1.12596087	#1.08820786		#PS zeta= 1.124 and delta=0.97423316 at 2.0 bar gives retrograde  #1.100
delta=1.00		#0.97423316	

if Find_Tg_at_P:
	#For Kier or Hassan or Condo:
	Tg_bisect_fwd=npy.zeros(len(P))
	Tg_bisect_bwd=npy.zeros(len(P))

if Find_Pg_at_T:
	#For Kier or Hassan or Condo:
	Pg_bisect_fwd=npy.zeros(len(T))
	Pg_bisect_bwd=npy.zeros(len(T))

if CondoFind_Tg_at_P:	
	#For Condo_Original:
	cTg_bisect_fwd=npy.zeros(len(P))
	cTg_bisect_bwd=npy.zeros(len(P))

if CondoFind_Pg_at_T:
	#For Condo_Original:
	cPg_bisect_fwd=npy.zeros(len(T))
	cPg_bisect_bwd=npy.zeros(len(T))

if Condo_Original:
	kwargs = {'z':cz,'epsilon_s':cepsilon_s,'epsilon_p':cepsilon_p,'zeta':czeta,'delta':cdelta,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'forward':forward,'backward':backward}

	if CondoFind_Tg_at_P:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			if forward:
				cTg_bisect_fwd[i]=CondoGlassTemperature_Original('fwd',P[i],Mp,Ms,**kwargs)
			if backward:
				cTg_bisect_bwd[i]=CondoGlassTemperature_Original('bwd',P[i],Mp,Ms,**kwargs)

	if CondoFind_Pg_at_T:
		for i in range(0,len(T)):
			print 'Iterating for T:', T[i], 'for bisect method'
			if forward:
				cPg_bisect_fwd[i]=CondoGlassPressure_Original('fwd',T[i],Mp,Ms,**kwargs)
			if backward:
				cPg_bisect_bwd[i]=CondoGlassPressure_Original('bwd',T[i],Mp,Ms,**kwargs)


if Kier or Hassan or Condo or Hassan_Var_Vol:
	kwargs = {'g':g,'epsilon_p':epsilon_p,'x':x,'zeta':zeta,'delta':delta,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'vh_pure_p':vh_pure_p,'vh_pure_s':vh_pure_s,'v_r':v_r,'Kier':Kier,'Hassan':Hassan,'Hassan_Var_Vol':Hassan_Var_Vol,'Condo':Condo,'forward':forward,'backward':backward}
	P_random=150.0
	if forward:
		Smax,S_infty=entropy_infty(P_random,Mp,Ms,**kwargs)
	if backward:
		Smax,S_infty=entropy_infty(P_random,Mp,Ms,**kwargs)

	kwargs = {'S_infty':S_infty,'g':g,'epsilon_p':epsilon_p,'x':x,'zeta':zeta,'delta':delta,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'vh_pure_p':vh_pure_p,'vh_pure_s':vh_pure_s,'v_r':v_r,'Kier':Kier,'Hassan':Hassan,'Hassan_Var_Vol':Hassan_Var_Vol,'Condo':Condo,'forward':forward,'backward':backward}

	if Find_Tg_at_P:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			if forward:
				Tg_bisect_fwd[i] = GlassTemperature('fwd',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_bisect_bwd[i] = GlassTemperature('bwd',P[i],Mp,Ms,**kwargs)

	if Find_Pg_at_T:
		for i in range(0,len(T)):
			print 'Iterating for T:', T[i], 'for bisect method'
			if forward:
				Pg_bisect_fwd[i] = GlassPressure('fwd',T[i],Mp,Ms,**kwargs)
			if backward:
				Pg_bisect_bwd[i] = GlassPressure('bwd',T[i],Mp,Ms,**kwargs)


#Prints:
if Condo_Original:
	if CondoFind_Tg_at_P:
		if forward:
			P,cTg_bisect_fwd=discard_zeros(P,cTg_bisect_fwd)
			print 'P is:', P
			print 'cTg_bisect_fwd is:', cTg_bisect_fwd
		if backward:
			P,cTg_bisect_bwd=discard_zeros(P,cTg_bisect_bwd)
			print 'P is:', P
			print 'cTg_bisect_bwd is:', cTg_bisect_bwd

	if CondoFind_Pg_at_T:
		if forward:
			cPg_bisect_fwd,T=discard_zeros(cPg_bisect_fwd,T)
			print 'cPg_bisect_fwd is:', cPg_bisect_fwd
			print 'T is:', T
		if backward:
			cPg_bisect_bwd,T=discard_zeros(cPg_bisect_bwd,T)
			print 'cPg_bisect_bwd is:', cPg_bisect_bwd
			print 'T is:', T

if Kier or Hassan or Condo or Hassan_Var_Vol:
	if Find_Tg_at_P:
		if forward:
			P,Tg_bisect_fwd=discard_zeros(P,Tg_bisect_fwd)
			print 'P is:', P
			print 'Tg_bisect_fwd is:', Tg_bisect_fwd
		if backward:
			P,Tg_bisect_bwd=discard_zeros(P,Tg_bisect_bwd)
			print 'P is:', P
			print 'Tg_bisect_bwd is:', Tg_bisect_bwd

	if Find_Pg_at_T:
		if forward:
			Pg_bisect_fwd,T=discard_zeros(Pg_bisect_fwd,T)
			print 'Pg_bisect_fwd is:', Pg_bisect_fwd
			print 'T is:', T
		if backward:
			Pg_bisect_bwd,T=discard_zeros(Pg_bisect_bwd,T)
			print 'Pg_bisect_bwd is:', Pg_bisect_bwd
			print 'T is:', T

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

plt.plot(P_exp,Tg_exp,color='b',marker='o',ls='',label='Tg_exp_condo',ms=markersize)
# plt.plot(P_exp_Condo,Tg_exp_Condo,color='k',marker='o',ls='',label='Tg_exp_condo',ms=markersize)

if Condo_Original:
	if CondoFind_Tg_at_P:
		if forward:
			plt.plot(P,cTg_bisect_fwd,color='g',lw=linewidth,ls='-',label='cTg_bisect_fwd')
		if backward:
			plt.plot(P,cTg_bisect_bwd,color='m',lw=linewidth,ls='-',label='cTg_bisect_bwd')

	if CondoFind_Pg_at_T:
		if forward:
			plt.plot(cPg_bisect_fwd,T,color='g',lw=linewidth,ls='-',label='cPg_bisect_fwd')
		if backward:
			plt.plot(cPg_bisect_bwd,T,color='m',lw=linewidth,ls='-',label='cPg_bisect_bwd')

if Kier or Hassan or Condo or Hassan_Var_Vol:
	if Find_Tg_at_P:
		if forward:
			plt.plot(P,Tg_bisect_fwd,color='k',marker='x',lw=linewidth,ls='-.',label='Tg_bisect_fwd')
		if backward:
			plt.plot(P,Tg_bisect_bwd,color='b',lw=linewidth,ls='-',label='Tg_bisect_bwd')

	if Find_Pg_at_T:
		if forward:
			plt.plot(Pg_bisect_fwd,T,color='k',marker='x',lw=linewidth,ls='-.',label='Pg_bisect_fwd')
		if backward:
			plt.plot(Pg_bisect_bwd,T,color='b',lw=linewidth,ls='-',label='Pg_bisect_bwd')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature Tg (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=1,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()


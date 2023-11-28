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
from plot_phi_Retro_improved import *

#Kier and Hassan:
def entropy_Kier_HassanOLD(P,T,phi_p,phi_s,Mp,Ms,**kwargs):
	
	# print 'Temperature is', T

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	vhm = delta*vhs
	alpha_p = alpha_p0*vhp/vhm
	alpha_s = alpha_s0*vhs/vhm
	Vratio=1.0 #hp/vhm
	#Following are My Theory Equations with, in general, v!=v_0:
	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))

	S=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/phi_p)*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((Vratio/alpha_s)*(phi_s/phi_p)*(ln(phi_s)))-((Vratio/alpha_p)*(ln(phi_p)))+((F*Vratio*epsilon_p)/(kB*T))-(ln(1-F)))

	return S

def entropy_Kier_Hassan(P,T,phi_p,phi_s,Rptilde,Mp,Ms,**kwargs):

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	# vhm = delta*vhs
	Vratio=1.0

	##########################################################################
	# Vratio=vhp/vhm
	# F_not valid=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))
	# S_2_wrong=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/phi_p)*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((Vratio/alpha_s)*(phi_s/phi_p)*(ln(phi_s)))-((Vratio/alpha_p)*(ln(phi_p)))+((F*Vratio*epsilon_p)/(kB*T))-(ln(1-F)))
	# S_1_wrong=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/phi_p)*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((Vratio/alpha_s)*(phi_s/phi_p)*(ln(phi_s)))-((Vratio/alpha_p)*(ln(phi_p)))+((F*Vratio*epsilon_p)/(kB*T))-(ln(1-F)))
	# S_term1_wrong=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/phi_p)*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p))))
	# S_term2_wrong=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/alpha_s)*(phi_s/phi_p)*(ln(phi_s))))
	# S_term3_wrong=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/alpha_p)*(ln(phi_p))))
	# S_term4_wrong=(Ppstar/(Rpstar*Tpstar))*(((F*Vratio*epsilon_p)/(kB*T)))
	# S_term5_wrong=(Ppstar/(Rpstar*Tpstar))*(-(ln(1-F)))
	# S_2_pure_wrong=(Ppstar/(Rpstar*Tpstar))*(-((1-Rptilde)*(ln(1-Rptilde))/Rptilde)-((ln(Rptilde))/alpha_p0)+((Tpstar/T)*g*(exp(-((Vratio*epsilon_p))/(kB*T)))/(1+g*exp(-((Vratio*epsilon_p))/(kB*T))))+((1/Vratio)*ln(1+g*exp(-(Vratio*epsilon_p)/(kB*T)))))
	###########################################################################
	
	# if P==0.101325:
	# 	print phi_s,phi_p,T,P

	#Following are My Theory Equations with, in general, v!=v_0:
	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))
	# S_1=(Ppstar/(Rpstar*Tpstar))*(-((vhp/(vhm*phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	# S_1=(Ppstar/(Rpstar*Tpstar))*(-((1/(delta_p*phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(delta_p*alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(delta_p*alpha_p))*(ln(phi_p)))+((F*epsilon_p)/(delta_p*kB*T))-(ln(1-F)))
	# S_1=(Ppstar/(Rpstar*Tpstar))*(-((vhp/(vhm*phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((vhp/(alpha_s*vhm))*(phi_s/phi_p)*(ln(phi_s)))-((vhp/(vhm*alpha_p))*(ln(phi_p)))+((F*Vratio*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	# In below equation v_r is taken as vhp:
	# Seems wrong: S_1=(Ppstar/(Rpstar*Tpstar))*(-((1/(phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((vhp/(alpha_s0))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p0))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))

	#Below assuming system to be in pure vhp state but phi_s present some how:
	alpha_s = alpha_s0*vhs/vhp
	S_1=(Ppstar/(Rpstar*Tpstar))*(-((1/(phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p0))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	S_2_slightly_approx_by_alpha_s0=(Ppstar/(Rpstar*Tpstar))*(-((1/(phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s0))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p0))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	S_2=S_2_slightly_approx_by_alpha_s0

	# For Pure:
	Vratio=1.0
	S_2_pure=(Ppstar/(Rpstar*Tpstar))*(-((1-Rptilde)*(ln(1-Rptilde))/Rptilde)-((ln(Rptilde))/alpha_p0)+((epsilon_p/(kB*T))*g*(exp(-((Vratio*epsilon_p))/(kB*T)))/(1+g*exp(-((Vratio*epsilon_p))/(kB*T))))+((1/Vratio)*ln(1+g*exp(-(Vratio*epsilon_p)/(kB*T)))))

	# S_2=S_2_pure

	# Successful Entropy:
	vhm = delta*vhs
	alpha_p = alpha_p0*vhp/vhm
	alpha_s = alpha_s0*vhs/vhm
	Vratio=1.0
	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))
	S_worked=(Ppstar/(Rpstar*Tpstar))*(-((Vratio/phi_p)*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((Vratio/alpha_s)*(phi_s/phi_p)*(ln(phi_s)))-((Vratio/alpha_p)*(ln(phi_p)))+((F*Vratio*epsilon_p)/(kB*T))-(ln(1-F)))
	# S_1=S_worked


	return S_1

def entropy_max_Kier_Hassan(direction,method,P,Mp,Ms,**kwargs):
		
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	#Initializing the array of densities.
	T=npy.linspace(300,3000,20)
	phi_p=npy.zeros(len(T))
	phi_s=npy.zeros(len(T))
	Rstilde=npy.zeros(len(T))
	Rptilde=npy.zeros(len(T))
	Rtilde=npy.zeros(len(T))
	S=npy.zeros(len(T))
	
	if Kier:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta,'method':'disparate','Kier':Kier}
	if Hassan:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta,'method':'disparate','Hassan':Hassan}
	
	for i in range(0,len(T)):	

		if method=='bisect':
			_,_,phi_p[i],phi_s[i],Rptilde[i],Rstilde[i]=binaryPhaseEquilibriumCHV_bisect(direction,P,T[i],Mp,Ms,**kwargs)
			Rstilde[i]=phi_pure_k_bisect(direction,P,T[i],Ms,Psstar,Tsstar,Rsstar)
			Rtilde[i]=phi_p[i]+phi_s[i]
		elif method=='nsolve':
			_,_,phi_p[i],phi_s[i],Rptilde[i],Rstilde[i]=binaryPhaseEquilibriumCHV_nsolve_first_time(direction,P,T[i],Mp,Ms,**kwargs)
			Rtilde[i]=phi_p[i]+phi_s[i]

	for i in range(0,len(T)):

		kwargs = {'g':g,'epsilon_p':epsilon_p,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta}
		S[i]=entropy_Kier_Hassan(P,T[i],phi_p[i],phi_s[i],Rptilde[i],Mp,Ms,**kwargs)

	S_max=npy.max(S)
	S_infty=S[-1]

	return S_max,S_infty

def GlassTemperatureCriterion_Kier_Hassan(T,direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,delta,g,epsilon_p,x,Smax,S_infty):
	
	if Kier:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta,'method':'disparate','Kier':Kier}
	if Hassan:
		kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta,'method':'disparate','Hassan':Hassan}

	if method=='bisect':
		_,_,phi_p,phi_s,Rptilde,Rstilde=binaryPhaseEquilibriumCHV_bisect(direction,P,T,Mp,Ms,**kwargs)
		Rstilde=phi_pure_k_bisect(direction,P,T,Ms,Psstar,Tsstar,Rsstar)
		Rtilde=phi_p+phi_s
	elif method=='nsolve':
		_,_,phi_p,phi_s,Rptilde,Rstilde=binaryPhaseEquilibriumCHV_nsolve_first_time(direction,P,T,Mp,Ms,**kwargs)
		Rtilde=phi_p+phi_s

	kwargs = {'g':g,'epsilon_p':epsilon_p,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':delta}
	S=entropy_Kier_Hassan(P,T,phi_p,phi_s,Rptilde,Mp,Ms,**kwargs)
	
	criterion=S-0.310707*0.8708171			#0.310707*0.8708171
	# criterion=S-x*Smax
	# criterion=S-x*S_infty

	return criterion

def GlassTemperature_Kier_Hassan(direction,method,P,Mp,Ms,**kwargs):
		
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	max_Tg=500
	min_Tg=200
	step_Tg=100

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
			Tg = bisect(GlassTemperatureCriterion_Kier_Hassan,i,i+step,args=(direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,delta,g,epsilon_p,x,Smax,S_infty))
		except:
			# print 'No value found'
			pass
		if Tg!=0.0:
			print 'Hurry! Tg is:', Tg, 'for direction', direction
			break
	if Tg==0.0:
		print 'Program Failed to get value of Tg in given bisect range in direction', direction

	return Tg

#Condo Modified:
def entropy_Condo_Modified(P,T,phi_p,phi_s,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	alpha_p0=(Mp*Ppstar)/(kB*Rpstar*Tpstar)
	alpha_s0=(Ms*Psstar)/(kB*Rsstar*Tsstar)

	# Tps_star=zeta*math.sqrt(Tpstar*Tsstar)

	Rtilde=phi_s+phi_p
	vhm=(phi_s/Rtilde)*vhs+(phi_p/Rtilde)*vhp
	alpha_pm = alpha_p0
	alpha_sm = alpha_s0					

	r=Rtilde/((phi_p/alpha_pm)+(phi_s/alpha_sm))					

	Fp=(((z-2)*exp(-epsilon_p/(kB*T)))/(1+(z-2)*exp(-epsilon_p/(kB*T))))
	Fs=(((z-2)*exp(-epsilon_s/(kB*T)))/(1+(z-2)*exp(-epsilon_s/(kB*T))))

	S_1=(Ppstar/(Rpstar*Tpstar))*(-((1-Rtilde)*ln(1-Rtilde)/phi_p)-(ln(phi_p)/alpha_pm)+(ln(alpha_pm)/alpha_pm)-(phi_s*ln(phi_s)/phi_p)-(1+((ln(2/z)-1)/r))-(((alpha_pm-2)/alpha_pm)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T))))-((phi_s/phi_p)*(-ln(1-Fs)+(Fs*epsilon_s/(kB*T))+1+((ln(2/z)-1)/r))))

	# vtilde=1/Rtilde
	# S_2=(Ppstar/(Rpstar*Tpstar))*(1/(phi_p/Rtilde))*(-1*((vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/r+((phi_s/Rtilde)/alpha_sm)*ln((phi_s/Rtilde)/alpha_sm)+((phi_p/Rtilde)/alpha_pm)*ln((phi_p/Rtilde)/alpha_pm)+1+(ln(2/z)-1)/r+((phi_s/Rtilde)/alpha_sm)*(alpha_sm-2)*(ln(1-Fs)-(Fs*epsilon_s/(kB*T)))+((phi_p/Rtilde)/alpha_pm)*(alpha_pm-2)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T)))))

	return S_1

def GlassTemperatureCriterion_Condo_Modified(T,direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,z,epsilon_p,epsilon_s):
	
	kwargs = {'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta,'delta':0.0,'method':'disparate','Condo':Condo}

	if method=='bisect':
		_,_,phi_p,phi_s,Rptilde,Rstilde=binaryPhaseEquilibriumCHV_bisect(direction,P,T,Mp,Ms,**kwargs)
		Rstilde=phi_pure_k_bisect(direction,P,T,Ms,Psstar,Tsstar,Rsstar)
		Rtilde=phi_p+phi_s
	elif method=='nsolve':
		_,_,phi_p,phi_s,Rptilde,Rstilde=binaryPhaseEquilibriumCHV_nsolve_first_time(direction,P,T,Mp,Ms,**kwargs)
		Rtilde=phi_p+phi_s

	kwargs = {'z':z,'epsilon_p':epsilon_p,'epsilon_s':epsilon_s,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta}
	S=entropy_Condo_Modified(P,T,phi_p,phi_s,Mp,Ms,**kwargs)
	criterion=S

	return criterion

def GlassTemperature_Condo_Modified(direction,method,P,Mp,Ms,**kwargs):
		
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	max_Tg=400
	min_Tg=200
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
			Tg = bisect(GlassTemperatureCriterion_Condo_Modified,i,i+step,args=(direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,z,epsilon_p,epsilon_s))
		except:
			# print 'No value found'
			pass
		if Tg!=0.0:
			print 'Hurry! Tg is:', Tg, 'for direction', direction
			break
	if Tg==0.0:
		print 'Program Failed to get value of Tg in given bisect range in direction', direction

	return Tg

#Condo Original:
def entropy_Original(T,direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,z,epsilon_p,epsilon_s):
			
	# print 'Temperature is', T
	
	kwargs = {'z':z,'epsilon_p':epsilon_p,'epsilon_s':epsilon_s,'Ppstar':Ppstar,'Tpstar':Tpstar,'Rpstar':Rpstar,'Psstar':Psstar,'Tsstar':Tsstar,'Rsstar':Rsstar,'zeta':zeta}

	if method=='bisect':
		phi_s,phi_p,Rtilde,Rstild,Rptilde=binaryPhaseEquilibriumCondo_Original_bisect(direction,P,T,Mp,Ms,**kwargs)
	elif method=='nsolve':
		phi_s,phi_p,Rtilde,Rstild=binaryPhaseEquilibriumCondo_Original_nsolve_first_time(direction,P,T,Mp,Ms,**kwargs)

	rp = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	rs = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	r=1/(phi_s/rs+phi_p/rp)

	vtilde=1/Rtilde

	Fp=((z-2)*exp(-epsilon_p/(kB*T)))/(1+(z-2)*exp(-epsilon_p/(kB*T)))
	Fs=((z-2)*exp(-epsilon_s/(kB*T)))/(1+(z-2)*exp(-epsilon_s/(kB*T)))

	S=-1*((vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/r+(phi_s/rs)*ln(phi_s/rs)+(phi_p/rp)*ln(phi_p/rp)+1+(ln(2/z)-1)/r+(phi_s/rs)*(rs-2)*(ln(1-Fs)-(Fs*epsilon_s/(kB*T)))+(phi_p/rp)*(rp-2)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T))))

	return S

def GlassTemperature_Original(direction,method,P,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	max_Tg=400
	min_Tg=200
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
			Tg = bisect(entropy_Original,i,i+step,args=(direction,method,P,Mp,Ms,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,zeta,z,epsilon_p,epsilon_s))
		except:
			# print 'No value found'
			pass
		if Tg!=0.0:
			print 'Hurry! Tg is:', Tg, 'for direction', direction
			break
	if Tg==0.0:
		print 'Program Failed to get value of Tg in given bisect range in direction', direction

	return Tg

####################################################################################################

# Tg_exp=[373.0,328.3,308.2,305.0,304.5]	
# P_exp=[0.101325,5.90,7.04,8.36,10.38]	

# Tg_exp_Condo=[374.1,339.2,323.9,308.8]
# P_exp_Condo=[0.008250,3.64,4.88,6.06]

if Condo_Original or Condo:
	cepsilon_s=0
	cepsilon_p=7151.0
	cz=5.0

if Kier or Hassan:
		
	kg=1.6674132169170892
	kepsilon_p=8012.605528610325
	kx=0.3107070707070707

	hg=1.6674132169170892
	hepsilon_p=8012.605528610325
	hx=0.3107070707070707

forward=True
backward=False
bisect_method=True
nsolve_method=False 

#Initializing the array of densities.
P=npy.linspace(0.101325,11,1)

#For Kier or Hassan or Condo:
Tg_bisect_fwd=npy.zeros(len(P))
Tg_bisect_bwd=npy.zeros(len(P))
Tg_nsolve_fwd=npy.zeros(len(P))
Tg_nsolve_bwd=npy.zeros(len(P))

#For Condo_Original:
cTg_bisect_fwd=npy.zeros(len(P))
cTg_bisect_bwd=npy.zeros(len(P))
cTg_nsolve_fwd=npy.zeros(len(P))
cTg_nsolve_bwd=npy.zeros(len(P))

if Condo_Original:
	kwargs = {'z':cz,'epsilon_p':cepsilon_p,'epsilon_s':cepsilon_s,'Ppstar':cPpstar,'Tpstar':cTpstar,'Rpstar':cRpstar,'Psstar':cPsstar,'Tsstar':cTsstar,'Rsstar':cRsstar,'zeta':czeta}

	if bisect_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			if forward:
				cTg_bisect_fwd[i]=GlassTemperature_Original('fwd','bisect',P[i],Mp,Ms,**kwargs)
			if backward:
				cTg_bisect_bwd[i]=GlassTemperature_Original('bwd','bisect',P[i],Mp,Ms,**kwargs)
		
	if nsolve_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for nsolve method'
			if forward:
				cTg_nsolve_fwd[i]=GlassTemperature_Original('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			if backward:
				cTg_nsolve_bwd[i]=GlassTemperature_Original('bwd','nsolve',P[i],Mp,Ms,**kwargs)

if Condo:
	kwargs = {'z':cz,'epsilon_p':cepsilon_p,'epsilon_s':cepsilon_s,'Ppstar':cPpstar,'Tpstar':cTpstar,'Rpstar':cRpstar,'Psstar':cPsstar,'Tsstar':cTsstar,'Rsstar':cRsstar,'zeta':czeta,'Condo':Condo}

	if bisect_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			if forward:
				Tg_bisect_fwd[i] = GlassTemperature_Condo_Modified('fwd','bisect',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_bisect_bwd[i] = GlassTemperature_Condo_Modified('bwd','bisect',P[i],Mp,Ms,**kwargs)
			
	if nsolve_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for nsolve method'
			if forward:
				Tg_nsolve_fwd[i]=GlassTemperature_Condo_Modified('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_nsolve_bwd[i]=GlassTemperature_Condo_Modified('bwd','nsolve',P[i],Mp,Ms,**kwargs)

if Kier:
	kwargs = {'g':kg,'epsilon_p':kepsilon_p,'x':kx,'Ppstar':kPpstar,'Tpstar':kTpstar,'Rpstar':kRpstar,'Psstar':kPsstar,'Tsstar':kTsstar,'Rsstar':kRsstar,'zeta':kzeta,'delta':kdelta,'Kier':Kier}

	if bisect_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			Smax,S_infty=entropy_max_Kier_Hassan('fwd','bisect',P[i],Mp,Ms,**kwargs)
			kwargs = {'Smax':Smax,'S_infty':S_infty,'g':kg,'epsilon_p':kepsilon_p,'x':kx,'Ppstar':kPpstar,'Tpstar':kTpstar,'Rpstar':kRpstar,'Psstar':kPsstar,'Tsstar':kTsstar,'Rsstar':kRsstar,'zeta':kzeta,'delta':kdelta,'Kier':Kier}
			if forward:
				Tg_bisect_fwd[i] = GlassTemperature_Kier_Hassan('fwd','bisect',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_bisect_bwd[i] = GlassTemperature_Kier_Hassan('bwd','bisect',P[i],Mp,Ms,**kwargs)
			
	if nsolve_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for nsolve method'
			Smax,S_infty=entropy_max_Kier_Hassan('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			kwargs = {'Smax':Smax,'S_infty':S_infty,'g':kg,'epsilon_p':kepsilon_p,'x':kx,'Ppstar':kPpstar,'Tpstar':kTpstar,'Rpstar':kRpstar,'Psstar':kPsstar,'Tsstar':kTsstar,'Rsstar':kRsstar,'zeta':kzeta,'delta':kdelta,'Kier':Kier}
			if forward:
				Tg_nsolve_fwd[i]=GlassTemperature_Kier_Hassan('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_nsolve_bwd[i]=GlassTemperature_Kier_Hassan('bwd','nsolve',P[i],Mp,Ms,**kwargs)

if Hassan:
	kwargs = {'g':hg,'epsilon_p':hepsilon_p,'x':hx,'Ppstar':hPpstar,'Tpstar':hTpstar,'Rpstar':hRpstar,'Psstar':hPsstar,'Tsstar':hTsstar,'Rsstar':hRsstar,'zeta':hzeta,'delta':hdelta,'Hassan':Hassan}

	if bisect_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for bisect method'
			Smax,S_infty=entropy_max_Kier_Hassan('fwd','bisect',P[i],Mp,Ms,**kwargs)
			kwargs = {'Smax':Smax,'S_infty':S_infty,'g':hg,'epsilon_p':hepsilon_p,'x':hx,'Ppstar':hPpstar,'Tpstar':hTpstar,'Rpstar':hRpstar,'Psstar':hPsstar,'Tsstar':hTsstar,'Rsstar':hRsstar,'zeta':hzeta,'delta':hdelta,'Hassan':Hassan}
			if forward:
				Tg_bisect_fwd[i] = GlassTemperature_Kier_Hassan('fwd','bisect',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_bisect_bwd[i] = GlassTemperature_Kier_Hassan('bwd','bisect',P[i],Mp,Ms,**kwargs)
			
	if nsolve_method:
		for i in range(0,len(P)):
			print 'Iterating for P:', P[i], 'for nsolve method'
			Smax,S_infty=entropy_max_Kier_Hassan('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			kwargs = {'Smax':Smax,'S_infty':S_infty,'g':hg,'epsilon_p':hepsilon_p,'x':hx,'Ppstar':hPpstar,'Tpstar':hTpstar,'Rpstar':hRpstar,'Psstar':hPsstar,'Tsstar':hTsstar,'Rsstar':hRsstar,'zeta':hzeta,'delta':hdelta,'Hassan':Hassan}
			if forward:
				Tg_nsolve_fwd[i]=GlassTemperature_Kier_Hassan('fwd','nsolve',P[i],Mp,Ms,**kwargs)
			if backward:
				Tg_nsolve_bwd[i]=GlassTemperature_Kier_Hassan('bwd','nsolve',P[i],Mp,Ms,**kwargs)

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

plt.plot(P_exp,Tg_exp,'b',marker='o',ls='',label='Tg_exp_condo',ms=markersize)
# plt.plot(P_exp_Condo,Tg_exp_Condo,'k',marker='o',ls='',label='Tg_exp_condo',ms=markersize)

if Condo_Original:
	if bisect_method:
		if forward:
			plt.plot(P,cTg_bisect_fwd,'k',color='g',lw=linewidth,ls='-',label='cTg_bisect_fwd')
		if backward:
			plt.plot(P,cTg_bisect_bwd,'k',color='m',lw=linewidth,ls='-',label='cTg_bisect_bwd')
	if nsolve_method:
		if forward:
			plt.plot(P,cTg_nsolve_fwd,'k',color='g',lw=linewidth,ls='-.',label='cTg_nsolve_fwd')
		if backward:
			plt.plot(P,cTg_nsolve_bwd,'k',color='m',lw=linewidth,ls='-.',label='cTg_nsolve_bwd')

if Kier or Hassan or Condo:
	if bisect_method:
		if forward:
			plt.plot(P,Tg_bisect_fwd,'k',color='k',lw=linewidth,ls='-',label='Tg_bisect_fwd')
		if backward:
			plt.plot(P,Tg_bisect_bwd,'k',color='b',lw=linewidth,ls='-',label='Tg_bisect_bwd')
	if nsolve_method:
		if forward:
			plt.plot(P,Tg_nsolve_fwd,'k',color='k',lw=linewidth,ls='-.',label='Tg_nsolve_fwd')
		if backward:
			plt.plot(P,Tg_nsolve_bwd,'k',color='b',lw=linewidth,ls='-.',label='Tg_nsolve_bwd')

plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel(r'Glass Temperature Tg (K)',fontsize=axis_size)
#plt.axis([300,500,0,1.5])
plt.legend(loc=1,fontsize=size,numpoints=1)
plt.subplots_adjust(bottom=0.3)

# figPUREPS.savefig('./'+output_folder+r'\pure_PMMA_TD_Property'+img_extension,dpi=img_dpi)

plt.show()


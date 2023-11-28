# Author: Kier von Konigslow
# Date: May 2017
#
# Description: The purpose of this file is to plot the saturated PS/CO2 binary mixture information
#			   based on both experimental and theoretical data.
#

import os,sys,math,csv,matplotlib.pyplot as plt,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from calculatePureVariables import calculateNewMolecularParameters
from wrapperFunctions import calculateBinarySolubilitySwelling
from calculateBinaryResidual import calculateBinarySSQ
from calculateBinaryVariablesCHV import *
from Parameters_of_Different_Polymers import *
from Parameters_for_Mixtures_and_Tg import *

Polymer_Type='PC'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}
Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper,'Kier':True,'Hassan':False,'Condo':False,'Condo_Original':False}
cepsilon_s,cepsilon_p,cz,czeta,epsilon_p,g,x,delta,zeta=Parameters_for_Mixtures_and_Tg(**kwargs)

zeta=1.17643522	#1.12272978	#1.02948288		#1.15666407#1.1684#1.135			#1.10752004			#1.08621732			#1.135
delta=0.97	#0.73223348	#0.90455520		#0.58402055#0.5#0.74264552			#0.90758208			#1.05473203			#1.00

P0 = npy.linspace(min(P0_X),max(P0_X),10)
# P0 = npy.linspace(4.5,20,10)

# print P0

T_scanner=T0_X[0]
T1=T0_X[0]

j=0
for i in range(0,len(T0_X)):
	if T_scanner==T0_X[i]:
		# P0_X_T1=P0_X
		# T0_X_T1=T0_X
		# X0_X_T1=X0_X
		# P0_S_T1=P0_S
		# T0_S_T1=T0_S
		# X0_S_T1=S0_S
		pass
	elif j==0:
			T_scanner=T0_X[i]
			T2=T0_X[i]
			index1=i
			P0_X_T1=P0_X[:index1]
			T0_X_T1=T0_X[:index1]
			X0_X_T1=X0_X[:index1]
			P0_S_T1=P0_S[:index1]
			T0_S_T1=T0_S[:index1]
			X0_S_T1=S0_S[:index1]

			P0_X_T2=P0_X[index1:]
			T0_X_T2=T0_X[index1:]
			X0_X_T2=X0_X[index1:]
			P0_S_T2=P0_S[index1:]
			T0_S_T2=T0_S[index1:]
			X0_S_T2=S0_S[index1:]

			j=j+1

	elif j==1:
			T_scanner=T0_X[i]
			T3=T0_X[i]
			index2=i
			P0_X_T2=P0_X[index1:index2]
			T0_X_T2=T0_X[index1:index2]
			X0_X_T2=X0_X[index1:index2]
			P0_S_T2=P0_S[index1:index2]
			T0_S_T2=T0_S[index1:index2]		
			X0_S_T2=S0_S[index1:index2]

			P0_X_T3=P0_X[index2:]
			T0_X_T3=T0_X[index2:]
			X0_X_T3=X0_X[index2:]
			P0_S_T3=P0_S[index2:]
			T0_S_T3=T0_S[index2:]
			X0_S_T3=S0_S[index2:]
			j=j+1

# print T0_X_T1
# print T0_X_T2
# print T0_X_T3
# print T0_X

#===================================================================================
#Setting global plotting parameters.
#===================================================================================

verbose = True
DHV = True
OSL = False
TSL = False

#===================================================================================
#Determining theoretical solubility/swelling.
#===================================================================================

#Set = [ zps , delta ]

if DHV:
	#Fit to all data, X only.
	# solu = [1.02127006,0.88102354] #Kiers values
	# solu = [0.92161759,0.88399214] #Hassan: My values

	# zeta = solu[0]
	# delta = solu[1]
	print('DHV mixture parameters zeta = {} and delta = {}.'.format(zeta,delta))
	
	gammas,vhs,epsilons = calculateNewMolecularParameters(Psstar,Tsstar,Rsstar,Ms)
	vh = delta*vhs/NA
	print('The hole volume is vh = {}.'.format(vh))
	'''
	P=[6.71321,8.25575333,9.79829667,11.34084,12.88338333,14.42592667,15.96847,17.51101333,19.05355667,20.5961,6.71321,8.25575333,9.79829667,11.34084,12.88338333,14.42592667,15.96847,17.51101333,19.05355667,20.5961,6.71321,8.25575333,9.79829667,11.34084,12.88338333,14.42592667,15.96847,17.51101333,19.05355667,20.5961]
	T=[403.0,403.0,403.0,403.0,403.0,403.0,403.0,403.0,403.0,403.0,423.0,423.0,423.0,423.0,423.0,423.0,423.0,423.0,423.0,423.0,463.0,463.0,463.0,463.0,463.0,463.0,463.0,463.0,463.0,463.0,]
	phip_Kier=  [0.867053639,0.860683885,0.854515717,0.848575914,0.842889614,0.837477997,0.832355469,0.827526852,0.822985508,0.818713407,0.856714353,0.850712756,0.844850162,0.83914095,0.833598238,0.828232945,0.82305273,0.818060983,0.813256049,0.808630919,0.833303198,0.827796979,0.822353196,0.816975093,0.81166505,0.806424339,0.801252878,0.796148999,0.791109247,0.786128226]
	phis_Kier=  [0.025733468,0.031316214,0.036732952,0.041962729,0.046985968,0.051786309,0.056352834,0.060682294,0.064780568,0.068662589,0.023757947,0.028979671,0.034089007,0.039075046,0.043927906,0.048639458,0.053204139,0.057619734,0.061887965,0.06601471,0.021026671,0.025744681,0.030414016,0.035032446,0.039598404,0.044111156,0.048570987,0.052979367,0.057339093,0.061654404]
	phip_Hassan=[0.10280427,0.148365207,0.197750274,0.242414664,0.280927024,0.313817777,0.341948986,0.366136186,0.387082972,0.405384283,0.098290273,0.128295347,0.16618883,0.204986622,0.241101058,0.273536785,0.302334619,0.327861034,0.350550515,0.370819651,0.099765756,0.118432669,0.141271225,0.167235291,0.19447166,0.221472119,0.247361768,0.271747975,0.294522T1,0.315724819]
	phis_Hassan=[0.334007861,0.392458527,0.410975285,0.411528213,0.404957371,0.395794208,0.385983326,0.376375798,0.367324678,0.3589422,0.268207313,0.333458125,0.369285271,0.383932502,0.387020871,0.384268508,0.378669427,0.37177747,0.364407482,0.356985605,0.189709844,0.240274671,0.28109079,0.309669572,0.327443375,0.337209754,0.341506636,0.34219498,0.340536048,0.337361067]
	mu_pure_s=[-21850.46194,-21239.45021,-20747.9747,-20340.88719,-19996.42291,-19700.05705,-19441.49637,-19213.07213,-19008.83972,-18824.07081,-23044.0052,-22387.13312,-21855.14539,-21411.1679,-21032.41175,-20703.74098,-20414.52961,-20156.97977,-19925.16171,-19714.4407,-25468.09947,-24723.93824,-24115.94684,-23603.82913,-23162.73779,-22776.23935,-22432.87695,-22124.32766,-21844.34292,-21588.10481]
	phis0_Hassan=[0.0705918488002,	0.0893362473962,	0.109209792647,	0.130278686765,	0.152564496502,	0.176007900101,	0.200425357396,	0.225476094203,	0.250669733583,	0.275437424049,	0.0658640219935,	0.082815456286,	0.100513495885,	0.118965304993,	0.138150252433,	0.158008278575,	0.178428647412,	0.199243111461,	0.220228727396,	0.241124188507,	0.0583736442586,	0.0727767419612,	0.0875374019288,	0.102635237587,	0.118039210176,	0.133706298753,	0.149580808686,	0.165594667751,	0.181668994408,	0.197717042154]
	phis0_Kier= [0.070591849,	0.089336247,	0.109209793,	0.130278687,0.152564497,	0.1760079,	0.200425357,	0.225476094,	0.250669734,	0.275437424,	0.065864022,	0.082815456,	0.100513496,	0.118965305,	0.138150252,	0.158008279,	0.178428647,	0.199243111,	0.220228727,	0.241124189,	0.058373644,	0.072776742,	0.087537402,	0.102635238,	0.11803921,	0.133706299,	0.149580809,	0.165594668,	0.181668994,	0.197717042]
	# phip0_Kier=[0.901155258835725,	0.902179861076194,	0.903187474462627,	0.904178575746837,	0.905153621330203,	0.906113048439133,	0.907057276215268,	0.907986706727838,	0.908901725914829,	0.909802704458972,	0.889317795681950,	0.890464018198803,	0.891590480887656,	0.892697773960780,	0.893786460696637,	0.894857079102351,	0.895910143447433,	0.896946145680699,	0.897965556741033,	0.898968827771527,	0.864368878000839,	0.865783151478616,	0.867170837607098,	0.868532839072539,	0.869870011452431,	0.871183166536206,	0.872473075352867,	0.873740470936429,	0.874986050856320,	0.876210479536607]

	EOS_Kier=npy.zeros(len(P))
	EOS_Hassan=npy.zeros(len(P))
	mu_Kier=npy.zeros(len(P))
	mu_Hassan=npy.zeros(len(P))

	for i in range(0,len(P)):
		EOS_Kier[i]=EOS_mix(phis_Kier[i],P[i],T[i],Mp,Ms,phip_Kier[i],delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
		EOS_Hassan[i]=EOS_mix(phis_Hassan[i],P[i],T[i],Mp,Ms,phip_Hassan[i],delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
		mu_Kier[i]=chemicalEquilibrium(phip_Kier[i],P[i],T[i],phis0_Kier[i],mu_pure_s[i],Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
		mu_Hassan[i]=chemicalEquilibrium(phip_Hassan[i],P[i],T[i],phis0_Hassan[i],mu_pure_s[i],Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
	
	print EOS_Kier
	print EOS_Hassan
	print mu_Kier
	print mu_Hassan
	'''
	Pmin = min(P0_X)
	Pmax = max(P0_X)
	Tmin = min(T0_X)
	Tmax = max(T0_X)
	print('The pressure range is {}-{}MPa and the temperature range is {}-{}K.'.format(Pmin,Pmax,Tmin,Tmax))

	result = calculateBinarySolubilitySwelling('CHV',P0,T1,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T1_DHV = result[2]
	Sw_T1_DHV = result[3]
	result = calculateBinarySolubilitySwelling('CHV',P0,T2,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T2_DHV = result[2]
	Sw_T2_DHV = result[3]
	result = calculateBinarySolubilitySwelling('CHV',P0,T3,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T3_DHV = result[2]
	Sw_T3_DHV = result[3]

if OSL:
	#Fit to all data, X only.
	solu = 0.9

	zeta = solu
	print('OSL mixture parameters zeta = {}.'.format(zeta))
	
	result = calculateBinarySolubilitySwelling('OSL',P0,T1,Mp,Ms,zeta=zeta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T1_OSL = result[2]
	Sw_T1_OSL = result[3]
	result = calculateBinarySolubilitySwelling('OSL',P0,T2,Mp,Ms,zeta=zeta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T2_OSL = result[2]
	Sw_T2_OSL = result[3]
	result = calculateBinarySolubilitySwelling('OSL',P0,T3,Mp,Ms,zeta=zeta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T3_OSL = result[2]
	Sw_T3_OSL = result[3]

if TSL:
	#Fit to all data, X only.
	solu = [0.9,4.7]

	zeta = solu[0]
	delta = solu[1]
	print('TSL mixture parameters zeta = {} and delta = {}.'.format(zeta,delta))

	result = calculateBinarySolubilitySwelling('OSL',P0,T1,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T1_TSL = result[2]
	Sw_T1_TSL = result[3]
	result = calculateBinarySolubilitySwelling('OSL',P0,T2,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T2_TSL = result[2]
	Sw_T2_TSL = result[3]
	result = calculateBinarySolubilitySwelling('OSL',P0,T3,Mp,Ms,zeta=zeta,delta=delta,Ppstar=Ppstar,Tpstar=Tpstar,Rpstar=Rpstar,Psstar=Psstar,Tsstar=Tsstar,Rsstar=Rsstar,method='disparate',verbose=verbose)
	ms_T3_TSL = result[2]
	Sw_T3_TSL = result[3]

#===================================================================================
#Plotting the PS/CO2 mixture results.
#===================================================================================

#Setting font size
axis_size = 20
title_size = 20
size = 16
label_size = 16
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size


#Setting saved image properties
img_extension = '.pdf'
img_dpi = None
output_folder = 'plot_mixture'

#Checking for existence of output directory. If such a directory doesn't exist, one is created.
if not os.path.exists('./'+output_folder):
    os.makedirs('./'+output_folder)

#Defining axes
solubility_axes = [5,25,0.0,0.13]
swelling_axes = [5,25,1.0,1.17]

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
markT3K = '^'
mark473K = 's'

#Linestules
ls453K = '-'
lsT3K = '--'
ls473K = ':'

#General line properties.
linewidth = 3
markersize = 8

if DHV:
    #Plotting the solubility of the PS+CO2 mixture.
    figX = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
    plt.plot(P0_X_T1,X0_X_T1,'k',marker=mark453K,ls='',label='Saturated PS/CO2 mixture at T1K',ms=markersize)
    plt.plot(P0_X_T2,X0_X_T2,'k',marker=markT3K,ls='',label='T2K',ms=markersize)
    plt.plot(P0_X_T3,X0_X_T3,'k',marker=mark473K,ls='',label='T3K',ms=markersize)
    plt.plot(P0,ms_T1_DHV,'k',ls=ls453K,label='Present theory T1K',lw=linewidth)
    plt.plot(P0,ms_T2_DHV,'k',ls=lsT3K,label='Present theory T2K',lw=linewidth)
    plt.plot(P0,ms_T3_DHV,'k',ls=ls473K,label='Present theory T3K',lw=linewidth)
    plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
    plt.ylabel('Gas solubility $\chi_s$ (g/g)',fontsize=axis_size)
    plt.legend(loc=2,fontsize=size,numpoints=1)
    # plt.axis(solubility_axes)
    figX.savefig('./'+output_folder+r'\bin_PS_CO2_Solubility'+img_extension,dpi=img_dpi)

'''
if DHV:
    #Plotting the swelling of the PS+CO2 mixture.
    figS = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
    plt.plot(P0_S_T1,S0_S_T1,'k',marker=mark453K,ls='',label='Saturated PS/CO2 mixture at T1K',ms=markersize)
    plt.plot(P0_S_T2,S0_S_T2,'k',marker=markT3K,ls='',label='T2K',ms=markersize)
    plt.plot(P0_S_T3,S0_S_T3,'k',marker=mark473K,ls='',label='T3K',ms=markersize)
    plt.plot(P0,Sw_T1_DHV,'k',ls=ls453K,label='Present theory T1K',lw=linewidth)
    plt.plot(P0,Sw_T2_DHV,'k',ls=lsT3K,label='Present theory T2K',lw=linewidth)
    plt.plot(P0,Sw_T3_DHV,'k',ls=ls473K,label='Present theory T3K',lw=linewidth)
    plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
    plt.ylabel('Swelling ratio $S_W$ $(cm^3/cm^3)$',fontsize=axis_size)
    plt.legend(loc=2,fontsize=size,numpoints=1)
    plt.axis(swelling_axes)
    figS.savefig('./'+output_folder+r'\bin_PS_CO2_Swelling'+img_extension,dpi=img_dpi)

#Plotting solubility at T1K.
figT1X = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_X_T1,X0_X_T1,'ok',label='Saturated PS/CO2 mixture at T1K',ms=markersize)
if DHV:
	plt.plot(P0,ms_T1_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,ms_T1_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,ms_T1_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Solubility $\chi_s$ (g/g)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(solubility_axes)
figT1X.savefig('./'+output_folder+r'\bin_PS_CO2_T1K_X'+img_extension,dpi=img_dpi)

#Plotting swelling at T1K.
figT1S = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_S_T1,S0_S_T1,'ok',label='Saturated PS/CO2 mixture at T1K',ms=markersize)
if DHV:
	plt.plot(P0,Sw_T1_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,Sw_T1_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,Sw_T1_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Swelling Ratio $S_W$ ($cm^3$/$cm^3$)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(swelling_axes)
figT1S.savefig('./'+output_folder+r'\bin_PS_CO2_T1K_S'+img_extension,dpi=img_dpi)

#Plotting solubility at T2K.
figT2X = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_X_T2,X0_X_T2,'ok',label='Saturated PS/CO2 mixture at T2K',ms=markersize)
if DHV:
	plt.plot(P0,ms_T2_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,ms_T2_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,ms_T2_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Solubility $\chi_s$ (g/g)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(solubility_axes)
figT2X.savefig('./'+output_folder+r'\bin_PS_CO2_T2K_X'+img_extension,dpi=img_dpi)

#Plotting swelling at T2K.
figT2S = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_S_T2,S0_S_T2,'ok',label='Saturated PS/CO2 mixture at T2K',ms=markersize)
if DHV:
	plt.plot(P0,Sw_T2_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,Sw_T2_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,Sw_T2_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Swelling Ratio $S_W$ ($cm^3$/$cm^3$)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(swelling_axes)
figT2S.savefig('./'+output_folder+r'\bin_PS_CO2_T2K_S'+img_extension,dpi=img_dpi)

#Plotting solubility at T3K.
figT3X = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_X_T3,X0_X_T3,'ok',label='Saturated PS/CO2 mixture at T3K',ms=markersize)
if DHV:
	plt.plot(P0,ms_T3_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,ms_T3_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,ms_T3_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Solubility $\chi_s$ (g/g)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(solubility_axes)
figT3X.savefig('./'+output_folder+r'\bin_PS_CO2_T3K_X'+img_extension,dpi=img_dpi)

#Plotting swelling at T3K.
figT3S = plt.figure(num=None, figsize=(12, 10), dpi=img_dpi, facecolor='w', edgecolor='k')
plt.plot(P0_S_T3,S0_S_T3,'ok',label='Saturated PS/CO2 mixture at T3K',ms=markersize)
if DHV:
	plt.plot(P0,Sw_T3_DHV,label=dhv_label,c=dhv_color,ls=dhv_line,lw=linewidth)
if TSL:
	plt.plot(P0,Sw_T3_TSL,label=tsl_label,c=tsl_color,ls=tsl_line,lw=linewidth)
if OSL:
	plt.plot(P0,Sw_T3_OSL,label=osl_label,c=osl_color,ls=osl_line,lw=linewidth)
plt.xlabel('Pressure P (MPa)',fontsize=axis_size)
plt.ylabel('Swelling Ratio $S_W$ ($cm^3$/$cm^3$)',fontsize=axis_size)
plt.legend(loc=2,fontsize=size,numpoints=1)
plt.axis(swelling_axes)
figT3S.savefig('./'+output_folder+r'\bin_PS_CO2_T3K_S'+img_extension,dpi=img_dpi)
'''

#Show plot windows.
plt.show()
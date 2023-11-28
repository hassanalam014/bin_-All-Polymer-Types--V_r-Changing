# Author: Kier von Konigslow
# Date: May 2017
#
# Description	: The purpose of this file is to estimate the multicomponent fluid parameters
#				  for the PS/CO2 binary mixture.
#

import os,sys,math,csv,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
from lmfit import minimize, Parameters, report_fit
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from calculateBinaryResidual import calculateBinaryResidualCHV
from Parameters_of_Different_Polymers import *

Polymer_Type='PC'
Solvent='CO2'
Parameters_Paper='Condo'
kwargs = {'Polymer_Type':Polymer_Type,'Solvent':Solvent,'Parameters_Paper':Parameters_Paper}

Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp=Parameters_of_Different_Polymers(**kwargs)
P0_X,T0_X,X0_X,P0_S,T0_S,S0_S=loadExperimentSwXData(**kwargs)

#Initializing the parameters.
params = Parameters()
#The following code sets up the model's parameters. It includes both fitting parameters and parameters that will remain fixed
#for the fitting. The values given are the inital guesses of fitting parameters and values of fixed parameters.
#           	(Name,		Value,		Vary?,	Min,	Max,	Expr)
params.add_many(('zeta',	1.125,		True,	0.0,	None,	None),
				('delta',	0.97,		True,	0.7,	1.0,	None),
				('Ppstar',	Ppstar,		False,	None,	None,	None),
				('Tpstar',	Tpstar,		False,	None,	None,	None),
				('Rpstar',	Rpstar,		False,	None,	None,	None),
				('Mp',		Mp,			False,	None,	None,	None),
				('Psstar',	Psstar,		False,	None,	None,	None),
				('Tsstar',	Tsstar,		False,	None,	None,	None),
				('Rsstar',	Rsstar,		False,	None,	None,	None),
				('Ms',		Ms,			False,	None,	None,	None),
				('fs',		0.0,		False,	None,	None,	None),
				('verbose',	False,		False,	None,	None,	None))

#fs=Weight of swelling residual in simultaneous solubility-swelling fit
#Running the Levenberg-Marquart algorithm on the residuals in order to do least squares fitting. This will return the fitted value of the RESIDUALS.
#These need to be added to the experimental datapints to find the fitted pressures.
print('For POLYMER: {} and SOLVENT: {}.'.format(Polymer_Type,Solvent))
print('Using {} parameters.'.format(Parameters_Paper))

fit = minimize(calculateBinaryResidualCHV,params,args=(P0_X,T0_X,X0_X,P0_S,T0_S,S0_S,'X','disparate'))

#Reporting the values of the parameters. NEED TO FIGURE OUT HOW TO PRINT THIS TO FILE.
report_fit(fit.params)
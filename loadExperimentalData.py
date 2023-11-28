# Author: Kier von Konigslow
# Date: May 2017
#
# Description	: The purpose of this module is to load the experimental solubility
#				  and swelling data for the PS/CO2 binary mixture.
#

import os,sys,math,numpy as npy
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadData import loadBinaryData

def loadExperimentSwXData(**kwargs):
    print 'who are you'
    for key,value in kwargs.items():
		exec "%s='%s'" % (key,value)

    print 'Polymer Type', Polymer_Type, 'Solvent', Solvent, 'Parameters from', Parameters_Paper

    if Polymer_Type=='PS' and Solvent=='CO2':
        #======================================================
        #Solubility Data
        #======================================================

        P0_X_403,T0_X_403,X0_X_403 = loadBinaryData('Data/403K_PS_CO2_X.csv')
        P0_X = P0_X_403
        T0_X = T0_X_403
        X0_X = X0_X_403

        P0_X_423,T0_X_423,X0_X_423 = loadBinaryData('Data/423K_PS_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_423),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_423),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_423),axis=0)

        P0_X_463,T0_X_463,X0_X_463 = loadBinaryData('Data/463K_PS_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_463),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_463),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_463),axis=0)

        #======================================================
        #Swelling Data
        #======================================================

        P0_S_403,T0_S_403,S0_S_403 = loadBinaryData('Data/403K_PS_CO2_Sw.csv')
        P0_S = P0_S_403
        T0_S = T0_S_403
        S0_S = S0_S_403

        P0_S_423,T0_S_423,S0_S_423 = loadBinaryData('Data/423K_PS_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_423),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_423),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_423),axis=0)

        P0_S_463,T0_S_463,S0_S_463 = loadBinaryData('Data/463K_PS_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_463),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_463),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_463),axis=0)

    if Polymer_Type=='PMMA' and Solvent=='CO2' and True:

        #======================================================
        #Solubility Data
        #======================================================

        P0_X_306,T0_X_306,X0_X_306 = loadBinaryData('Data/306K_PMMA_CO2_X.csv')
        P0_X = P0_X_306
        T0_X = T0_X_306
        X0_X = X0_X_306

        P0_X_315,T0_X_315,X0_X_315 = loadBinaryData('Data/315K_PMMA_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_315),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_315),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_315),axis=0)

        P0_X_332,T0_X_332,X0_X_332 = loadBinaryData('Data/332K_PMMA_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_332),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_332),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_332),axis=0)

        #======================================================
        #Swelling Data
        #======================================================

        P0_S_306,T0_S_306,S0_S_306 = loadBinaryData('Data/306K_PMMA_CO2_Sw.csv')
        P0_S = P0_S_306
        T0_S = T0_S_306
        S0_S = S0_S_306

        P0_S_315,T0_S_315,S0_S_315 = loadBinaryData('Data/315K_PMMA_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_315),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_315),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_315),axis=0)

        P0_S_332,T0_S_332,S0_S_332 = loadBinaryData('Data/332K_PMMA_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_332),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_332),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_332),axis=0)

    #New Paper of PMMA CO2, Paper Title: Simultaneous Measurement of Swelling and Sorption in a Supercritical CO2-Poly(methyl methacrylate) System

    if Polymer_Type=='PMMA' and Solvent=='CO2'and True:

        #======================================================
        #Solubility Data
        #======================================================

        P0_X_323,T0_X_323,X0_X_323 = loadBinaryData('Data/323K_PMMA_CO2_X.csv')
        
        P0_X = npy.concatenate((P0_X,P0_X_323),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_323),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_323),axis=0)

        # P0_X = P0_X_323
        # T0_X = T0_X_323
        # X0_X = X0_X_323

        P0_X_338,T0_X_338,X0_X_338 = loadBinaryData('Data/338K_PMMA_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_338),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_338),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_338),axis=0)

        P0_X_353,T0_X_353,X0_X_353 = loadBinaryData('Data/353K_PMMA_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_353),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_353),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_353),axis=0)

        #======================================================
        #Swelling Data
        #======================================================

        P0_S_323,T0_S_323,S0_S_323 = loadBinaryData('Data/323K_PMMA_CO2_Sw.csv')
        P0_S = P0_S_323
        T0_S = T0_S_323
        S0_S = S0_S_323

        P0_S_338,T0_S_338,S0_S_338 = loadBinaryData('Data/338K_PMMA_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_338),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_338),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_338),axis=0)

        P0_S_353,T0_S_353,S0_S_353 = loadBinaryData('Data/353K_PMMA_CO2_Sw.csv')
        P0_S = npy.concatenate((P0_S,P0_S_353),axis=0)
        T0_S = npy.concatenate((T0_S,T0_S_353),axis=0)
        S0_S = npy.concatenate((S0_S,S0_S_353),axis=0)


    if Polymer_Type=='PC' and Solvent=='CO2':
        
        #======================================================
        #Solubility Data
        #======================================================

        P0_X_313,T0_X_313,X0_X_313 = loadBinaryData('Data/313K_PC_CO2_X.csv')
        P0_X = P0_X_313
        T0_X = T0_X_313
        X0_X = X0_X_313

        P0_X_323,T0_X_323,X0_X_323 = loadBinaryData('Data/323K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_323),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_323),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_323),axis=0)

        P0_X_333,T0_X_333,X0_X_333 = loadBinaryData('Data/333K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_333),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_333),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_333),axis=0)

        P0_X_348,T0_X_348,X0_X_348 = loadBinaryData('Data/348K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_348),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_348),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_348),axis=0)

        P0_X_373,T0_X_373,X0_X_373 = loadBinaryData('Data/373K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_373),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_373),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_373),axis=0)

        # P0_X_373,T0_X_373,X0_X_373 = loadBinaryData('Data/373K_PC_CO2_X.csv')
        # P0_X = P0_X_373
        # T0_X = T0_X_373
        # X0_X = X0_X_373

        P0_X_398,T0_X_398,X0_X_398 = loadBinaryData('Data/398K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_398),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_398),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_398),axis=0)

        P0_X_423,T0_X_423,X0_X_423 = loadBinaryData('Data/423K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_423),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_423),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_423),axis=0)

        P0_X_448,T0_X_448,X0_X_448 = loadBinaryData('Data/448K_PC_CO2_X.csv')
        P0_X = npy.concatenate((P0_X,P0_X_448),axis=0)
        T0_X = npy.concatenate((T0_X,T0_X_448),axis=0)
        X0_X = npy.concatenate((X0_X,X0_X_448),axis=0)

        #======================================================
        #Swelling Data
        #======================================================

        P0_S_308,T0_S_308,S0_S_308 = loadBinaryData('Data/308K_PC_CO2_Sw.csv')
        P0_S = P0_S_308
        T0_S = T0_S_308
        S0_S = S0_S_308


    return P0_X,T0_X,X0_X,P0_S,T0_S,S0_S


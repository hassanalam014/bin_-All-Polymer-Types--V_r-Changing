#Sanchez-Lacombe parameters for the pure POLYMERS.

#======================================================
#Characteristic Parameters
#======================================================
# Pstar = [MPa]
# Tstar = [K]
# Rstar = [g/cm3]

#############################################################################################
# Polymer_Type='PMMA' #PS or PMMA or DME or LPP or BPP or PLA or LDPE
# Solvent='CO2' 
#############################################################################################

P_atm = 0.101325
M_infinity = 9E9

def Parameters_of_Different_Polymers(**kwargs):
	
    for key,value in kwargs.items():
		exec "%s='%s'" % (key,value)
    
    print 'Polymer Type', Polymer_Type, 'Solvent', Solvent, 'Parameters taken from', Parameters_Paper

    #PURE POLYMERS DATA:
    if Polymer_Type=='PS' and Parameters_Paper=='Condo':
        #Sanchez-Lacombe parameters for the pure POLYMER.
        # Condo PS Parameters:
        Mp = 196000.0
        Ppstar = 357.0
        Tpstar = 735.0
        Rpstar = 1.105
        # # Iterated Values From Pure PS Flexible System:
        # g=1.67 
        # epsilon_p=8013 
        # x=0.311

    if Polymer_Type=='PS' and Parameters_Paper=='Kier':

        # Iterated Values From Solubility of PS/CO2 System:
        Ppstar = 421.762455
        Tpstar = 687.788143
        Rpstar = 1.11768655


    if Polymer_Type=='PMMA' and Parameters_Paper=='Condo':
        #Sanchez-Lacombe parameters for the pure POLYMER.
        # Condo PMMA Parameters:
        Mp = 196000.0
        Ppstar = 503.0
        Tpstar = 696.0
        Rpstar = 1.269
        # Iterated Values From Pure PS Flexible System:
        # g=1.66 
        # epsilon_p=8094 
        # x=0.323

    if Polymer_Type=='PC':
        # For PC Ref: Zoller Paper, A Studey of PVT Relationships of Four Related Amorphous Polymers
        Mp = M_infinity
        Ppstar = 574.4             #From huge List of SL EOS Parameters
        Tpstar = 728.0             #From huge List of SL EOS Parameters
        Rpstar = 1.2925            #From huge List of SL EOS Parameters

    #PURE SOLVENT DATA:
    if Solvent=='CO2' and Parameters_Paper=='Kier':
       #Sanchez-Lacombe parameters for the pure SOLVENT.
        Ms = 44.01
        Psstar = 419.915520
        Tsstar = 341.772507
        Rsstar = 1.39745988
        Vsstar = 0.71558405

    if Solvent=='CO2' and Parameters_Paper=='Condo':
        #Sanchez-Lacombe parameters for the pure SOLVENT.
        # Condo CO2 Parameters:
        Ms = 44.01
        Psstar = 574.0
        Tsstar = 308.64
        Rsstar = 1.505

    #Experimental Data of Binary Mixtures:
    if Polymer_Type=='PS' and Solvent=='CO2':

        #Experiment Data PS;
        Tg_exp=[373.0,328.3,308.2,305.0,304.5]	
        P_exp=[0.101325,5.90,7.04,8.36,10.38]
        
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

        # Tg_exp_Condo=[374.1,339.2,323.9,308.8]
        # P_exp_Condo=[0.008250,3.64,4.88,6.06]

    if Polymer_Type=='PMMA' and Solvent=='CO2':

        #Experiment Data PMMA;
        Tg_exp=	[277.8,297.7,317.7,337.8,347.9]
        P_exp=	[3.75,5.14,5.70,5.09,3.73]

    if Polymer_Type=='PC' and Solvent=='CO2':

        #Experiment Data PC;
        Tg_exp = [423.02,413.201,401.982,394.405,383.186,417.975,414.327,407.875,399.741,396.934,368.8913,363.2814,348.1364,334.1149,268.80497,253.3957,242.7511]
        P_exp =	[0.0,1.506763545,2.82165807,4.375730258,5.586016853,0.0,0.699796046,1.641171158,2.642292555,3.24000873,5.675426033,6.377709608,8.006367128,9.231548498,9.992580308,8.976108173,7.98959784]

    return (Ppstar,Tpstar,Rpstar,Mp,Psstar,Tsstar,Rsstar,Ms,P_exp,Tg_exp)

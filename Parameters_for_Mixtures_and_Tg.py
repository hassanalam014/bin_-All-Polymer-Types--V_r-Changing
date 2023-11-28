#############################################################################################
# Polymer_Type='PMMA' #PS or PMMA or DME or LPP or BPP or PLA or LDPE
# Solvent='CO2' 
#############################################################################################

def Parameters_for_Mixtures_and_Tg(**kwargs):
	
    for key,value in kwargs.items():
		exec "%s='%s'" % (key,value)
    
    print 'Polymer Type', Polymer_Type, 'Solvent', Solvent, 'Parameters taken from', Parameters_Paper

    if Condo_Original or Condo:

        if Polymer_Type=='PS':
            cepsilon_s=0
            cepsilon_p=7151.0
            cz=5.0

        if Polymer_Type=='PMMA':
            cepsilon_s=0
            cepsilon_p=7443.0
            cz=5.0

        if Polymer_Type=='PC':
            cepsilon_s=0
            cepsilon_p=6247.0
            cz=4.0

    if Kier or Hassan or Hassan_Var_Vol:
            
        if Polymer_Type=='PS':
            g=1.67 
            epsilon_p=8013 
            x=0.311

        if Polymer_Type=='PMMA':	
            g=1.66 
            epsilon_p=8094 
            x=0.323

        if Polymer_Type=='PC':	
            g=0.84 
            epsilon_p=8273.0
            x=0.317

    if Kier:

        if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Kier':
            zeta=1.02127006
            delta=0.88102354

        if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Condo':
            zeta=1.08820786     #1.08820786(Original)   #1.11000815(from Xs condo formula)
            delta=0.97423316    #0.97423316(Original)   #1.00003292(from Xs condo formula)

        if Polymer_Type=='PMMA' and Solvent=='CO2' and Parameters_Paper=='Condo':
            zeta=1.13500002       #1.135      #1.13500002(from Xs condo formula)		#1.12272978	#1.07643522	#	#1.1684#1.135			#1.10752004		#1.08621732		#1.135		
            delta=1.00003613      #1.00       #1.00003613(from Xs condo formula)		#0.73223348	#1.0		#	#0.5#0.74264552			#0.90758208		#1.05473203		#1.00	

        if Polymer_Type=='PC' and Solvent=='CO2':
            zeta=1.06719283	#From Kier Solubility alpha
            # zeta=1.10694978 #From Condo Solubility alpha	
            delta=1.000

            #SuperCritical Solubility Fit by Condo Solubility alpha
            zeta=1.12596087
            delta=0.99999371

    if Condo or Condo_Original:

        if Polymer_Type=='PS' and Solvent=='CO2' and Parameters_Paper=='Condo':
            #Condo Zeta:
            czeta=1.1240			#Condo Values. Correlation value: 1.110, iteration value: 1.1240

        if Polymer_Type=='PMMA' and Solvent=='CO2' and Parameters_Paper=='Condo':
            #Condo Zeta:
            czeta=1.1350

        if Polymer_Type=='PC' and Solvent=='CO2':
            #Condo Zeta:
            czeta=1.000

    return cepsilon_s,cepsilon_p,cz,czeta,epsilon_p,g,x,delta,zeta
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
# import cmath
# import time
# from Parameters_of_Different_Polymers import *

def EOS_pure_k(phi_pure_k,P,T,Mk,Pkstar,Tkstar,Rkstar,vh_pure_k,v_r):

	alpha_0=vh_pure_k/v_r
	alpha_pure_k = (Pkstar*Mk)/(kB*Tkstar*Rkstar)

	#Equation of state in Pure solvent phase.
	EOS_pure_k = v_r*P/(kB*T)+((1.0/alpha_0)-(1.0/alpha_pure_k))*phi_pure_k+(1.0/alpha_0)*ln(1.0-phi_pure_k)+(Tkstar/T)*(phi_pure_k**2)

	return EOS_pure_k

def phi_pure_k_bisect(direction,P,T,Mk,Pkstar,Tkstar,Rkstar,vh_pure_k,v_r):
	
	step_size_of_phi=0.2
	range_max=1/step_size_of_phi

	if direction=='fwd':
		start=0
		end=int(range_max)+1
		step=1
		
	elif direction=='bwd':
		start=int(range_max)
		end=-1
		step=-1

	for i in range(start,end,step):
		k=i*step_size_of_phi

		if k+step*0.00001<0.0:
			guess1=0.00001
		elif k+step*0.00001>1.0:
			guess1=0.99999
		else:
			guess1=k+step*0.00001

		if k+step*0.00001+step*step_size_of_phi<0.0:
			guess2=0.00001
		elif k+step*0.00001+step*step_size_of_phi>1.0:
			guess2=0.99999
		else:
			guess2=k+step*0.00001+step*step_size_of_phi

		phi_pure_k=0.0
		try:
			phi_pure_k = bisect(EOS_pure_k,abs(guess1),abs(guess2),args=(P,T,Mk,Pkstar,Tkstar,Rkstar,vh_pure_k,v_r))
		except:
			# print 'No value found between', k, 'to', k+step_size_of_phi
			pass
		if phi_pure_k!=0.0:
			# print 'Hurry! phi_pure_k is:', phi_pure_k
			break
	
	if phi_pure_k==0.0:
		print 'Program Failed to get value of phi_pure_k in', direction, 'direction'

	return phi_pure_k

def EOS_mix(phi_mix_s,P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol):
	
	# print 'are you running by any chance'
	# vh_pure_p = kB*Tpstar/Ppstar
	# vh_pure_s = kB*Tsstar/Psstar

	alpha_pure_p = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_pure_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	Tps_star=zeta*math.sqrt(Tpstar*Tsstar)#*((T/Tpstar)**0.333333)

	if Kier or Hassan:	
		vhm = delta*v_r									#Hassan: This is definition of delta
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm							
		alpha_mix_p = alpha_pure_p				#Hassan: v_h is volume of hole in mixture or pure solvent
		alpha_mix_s = alpha_pure_s							

	elif Condo:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		alpha_mix_p = alpha_pure_p				#Hassan: v_h is volume of hole in mixture or pure solvent
		alpha_mix_s = alpha_pure_s							

	elif Hassan_Var_Vol:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm	
		alpha_mix_p = alpha_pure_p				#Hassan: v_h is volume of hole in mixture or pure solvent
		alpha_mix_s = alpha_pure_s	

		# pfactor=Mp/alpha_pure_p
		# sfactor=Ms/alpha_pure_s
		# vhm=(phi_mix_s*sfactor*vh_pure_s+phi_mix_p*pfactor*vh_pure_p)/(phi_mix_s*sfactor+phi_mix_p*pfactor)						

	alpha_0=vhm/v_r
	#Equation of state in mixture phase.
	EOS_mixture = v_r*P/(kB*T)+(1.0/alpha_0-1.0/alpha_mix_p)*phi_mix_p+(1.0/alpha_0-1.0/alpha_mix_s)*phi_mix_s+(1.0/alpha_0)*ln(1.0-phi_mix_p-phi_mix_s)+(Tsstar/T)*(phi_mix_s**2)+(Tpstar/T)*(phi_mix_p**2)+2*(Tps_star/T)*phi_mix_p*phi_mix_s

	return EOS_mixture

def phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol):
	
	step_size_of_phi=0.2
	range_max=1/step_size_of_phi

	if direction=='fwd':
		start=0
		end=int(range_max)+1
		step=1
		
	elif direction=='bwd':
		start=int(range_max)
		end=-1
		step=-1

	for i in range(start,end,step):
		k=i*step_size_of_phi

		
		if k+step*0.00001<0.0:
			guess1=0.00001
		elif k+step*0.00001>1.0:
			guess1=0.99999
		else:
			guess1=k+step*0.00001

		if k+step*0.00001+step*step_size_of_phi<0.0:
			guess2=0.00001
		elif k+step*0.00001+step*step_size_of_phi>1.0:
			guess2=0.99999
		else:
			guess2=k+step*0.00001+step*step_size_of_phi


		phi_mix_s=0.0
		try:
			phi_mix_s = bisect(EOS_mix,abs(guess1),abs(guess2),args=(P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol))
		except:
			# print 'No value found between', k, 'to', k+step_size_of_phi
			pass
		if phi_mix_s!=0.0 and phi_mix_s>0.0 :
			# print 'Hurry! phi_mix_s is:', phi_mix_s
			break
	if phi_mix_s==0.0 or phi_mix_s<0.0:
		phi_mix_s=0.0
		# print 'Program Failed to get value of phi_mix_s at assumed phi_mix_p in', direction, 'direction'
		pass

	return phi_mix_s

def chemicalEquilibrium(phi_mix_p,direction,P,T,phi_pure_s,mu_pure_s,Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol):

	# vh_pure_p = kB*Tpstar/Ppstar
	# vh_pure_s = kB*Tsstar/Psstar

	alpha_pure_p = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_pure_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	Tps_star=zeta*math.sqrt(Tpstar*Tsstar)#*((T/Tpstar)**0.333333)

	phi_mix_s=phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol)

	# ####Falto Addition Starts
	# phi_mix_p_falto=0.81545345
	# phi_mix_s_falto=phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p_falto,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
	# print 'phi_mix_s_falto is:', phi_mix_s_falto, 'at phi_mix_p_falto', phi_mix_p_falto, 'at temperature', T
	# ### Falto Addition Ends

	if Kier or Hassan:	
		vhm = delta*v_r									#Hassan: This is definition of delta
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm							
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s	
	elif Condo:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s							

	elif Hassan_Var_Vol:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm	
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s

		# pfactor=Mp/alpha_pure_p
		# sfactor=Ms/alpha_pure_s
		# vhm=(phi_mix_s*sfactor*vh_pure_s+phi_mix_p*pfactor*vh_pure_p)/(phi_mix_s*sfactor+phi_mix_p*pfactor)						

	#Equilibrium Condition is:
	mu_mix_s=mu_pure_s
	# print 'why are not you runnign'

	alpha_0=vhm/v_r

	#Solvent chemical potential in mixture phase.
	if Kier:
		residual_mu= -mu_mix_s+(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(0+ln(phi_mix_s))-((1.0/alpha_0))*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Kier
		# residual_mu= -mu_mix_s+alpha_mix_s*(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-(1.0/alpha_0)*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Including all factors
	elif Hassan or Hassan_Var_Vol:
		residual_mu= -mu_mix_s+alpha_mix_s*(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-(1.0/alpha_0)*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p-(((vh_pure_s-vh_pure_p)/(alpha_0*vhm))*(((1-phi_mix_p-phi_mix_s)*(phi_mix_p))/((phi_mix_p+phi_mix_s)**2)))) #Including all factors
	elif Condo:
		residual_mu= -mu_mix_s+alpha_mix_s*(-1+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Including all factors

	##################Falto Trying Start
	
	# phi_mix_p_falto=0.81545345
	# phi_mix_s_falto=phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p_falto,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar)
	# Rtilde_falto=phi_mix_p_falto+phi_mix_s_falto
	# residual_mu_falto= -mu_mix_s+alpha_mix_s*(-1+(1.0/alpha_mix_s)*(1+ln(phi_mix_s_falto))-ln(1-phi_mix_p_falto-phi_mix_s_falto)-2*(Tsstar/T)*phi_mix_s_falto-2*(Tps_star/T)*phi_mix_p_falto) #Including all factors
	# print 'residual is', residual_mu_falto, 'at phi_mix_p=', phi_mix_p_falto, 'at phi_mix_s=', phi_mix_s_falto, 'at Rtilde=', Rtilde_falto,'at T', T, 'and P',P
	
	######################Falto Trying Ends

	# print 'residual_mu is', residual_mu
	if phi_mix_s==0.0:
		# print 'Unfortunately program may pick phi_mix_s=', phi_mix_s
		residual_mu=-9000000000000000
	
	# print 'residual is:', residual_mu
	# print 'phi_mix_s for above residual is:', phi_mix_s, 'and phi_mix_p:', phi_mix_p
	return residual_mu

def chemicalEquilibrium_bisect(direction,P,T,phi_pure_s,mu_pure_s,Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol):
	
	step_size_refinement=[0.2,0.1,0.01]

	for m in range(0,len(step_size_refinement)):	

		step_size_of_phi=step_size_refinement[m]
		range_max=1/step_size_of_phi

		if direction=='fwd':
			start=0
			end=int(range_max)+1
			step=1
			# print 'forward'
			
		elif direction=='bwd':
			start=int(range_max)
			end=-1
			step=-1
			# print 'backward'

		for i in range(start,end,step):
			k=i*step_size_of_phi
					
			if k+step*0.001<0.0:
				guess1=0.001
			elif k+step*0.001>1.0:
				guess1=0.999
			else:
				guess1=k+step*0.001

			if k+step*0.001+step*step_size_of_phi<0.0:
				guess2=0.001
			elif k+step*0.001+step*step_size_of_phi>1.0:
				guess2=0.999
			else:
				guess2=k+step*0.001+step*step_size_of_phi

			phi_mix_p=0.0
			try:
				phi_mix_p = bisect(chemicalEquilibrium,abs(guess1),abs(guess2),args=(direction,P,T,phi_pure_s,mu_pure_s,Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol))
			except:
				# print 'No value found between', k+step*0.001, 'to', k+step*0.001+step*step_size_of_phi
				pass
			if phi_mix_p!=0.0 and phi_mix_p>0.0:
				# print 'Hurry! phi_mix_p is:', phi_mix_p
				break
		if phi_mix_p!=0.0 and phi_mix_p>0.0:
			# print 'Hurry! phi_mix_p is:', phi_mix_p
			break
		elif phi_mix_p==0.0 or phi_mix_p<0.0:
			# print 'Program Failed to get value of phi_mix_p in', direction, 'direction'
			pass
	if phi_mix_p==0.0 or phi_mix_p<0.0:
		# print 'this is failing'
		print 'Program Failed to get value of phi_mix_p in', direction, 'direction'

	phi_mix_s=phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol)
	
	#########################################################
	# vh_pure_p = kB*Tpstar/Ppstar
	# vh_pure_s = kB*Tsstar/Psstar
	alpha_pure_p = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_pure_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	Tps_star=zeta*math.sqrt(Tpstar*Tsstar)#*((T/Tpstar)**0.333333)

	phi_mix_s=phi_mix_s_bisect(direction,P,T,Mp,Ms,phi_mix_p,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol)

	if Kier or Hassan:	
		vhm = delta*v_r									#Hassan: This is definition of delta
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm							
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s	
	elif Condo:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s							

	elif Hassan_Var_Vol:
		Rtilde=phi_mix_s+phi_mix_p
		vhm=(phi_mix_s/Rtilde)*vh_pure_s+(phi_mix_p/Rtilde)*vh_pure_p
		# alpha_mix_p = alpha_pure_p*vh_pure_p/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_mix_s = alpha_pure_s*vh_pure_s/vhm	
		alpha_mix_p = alpha_pure_p
		alpha_mix_s = alpha_pure_s

		# pfactor=Mp/alpha_pure_p
		# sfactor=Ms/alpha_pure_s
		# vhm=(phi_mix_s*sfactor*alpha_pure_s+phi_mix_p*pfactor*alpha_pure_p)/(phi_mix_s*sfactor+phi_mix_p*pfactor)						

	#Equilibrium Condition is:
	mu_mix_s=mu_pure_s
	# print 'why are not you runnign'

	alpha_0=vhm/v_r
	
	#Solvent chemical potential in mixture phase.
	if Kier:
		residual_mu= -mu_mix_s+(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(0+ln(phi_mix_s))-((1.0/alpha_0))*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Kier
		# residual_mu= -mu_mix_s+alpha_mix_s*(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-(1.0/alpha_0)*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Including all factors
	elif Hassan or Hassan_Var_Vol:
		residual_mu= -mu_mix_s+alpha_mix_s*(-(1.0/alpha_0)+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-(1.0/alpha_0)*ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Including all factors
	elif Condo:
		residual_mu= -mu_mix_s+alpha_mix_s*(-1+(1.0/alpha_mix_s)*(1+ln(phi_mix_s))-ln(1-phi_mix_p-phi_mix_s)-2*(Tsstar/T)*phi_mix_s-2*(Tps_star/T)*phi_mix_p) #Including all factors

	print '#####################################################'
	print 'residual is:', residual_mu
	# print 'phi_mix_s for above residual is:', phi_mix_s
	print '#####################################################'
	##########################################################

	if phi_mix_s==0.0:
		print 'Program Failed to get value of phi_mix_s in', direction, 'direction'

	return phi_mix_p,phi_mix_s

def binaryPhaseEquilibriumCHV_bisect(direction,P,T,Mp,Ms,**kwargs):

	#Reference:
	# -p --> polymer
	# -s --> solvent

	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))
	
	#PURE FLUID PARAMETERS.
	alpha_pure_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	
	phi_pure_s=phi_pure_k_bisect(direction,P,T,Ms,Psstar,Tsstar,Rsstar,vh_pure_s,v_r)
	phi_pure_p=phi_pure_k_bisect(direction,P,T,Mp,Ppstar,Tpstar,Rpstar,vh_pure_p,v_r)

	# print phi_pure_s
	# print phi_pure_p
	# print Kier
	# print Hassan 
	# print Condo
	# print delta
	# print zeta

	#CHECKING IF PURE VOLUME FRACTION RESULTS ARE VALID.
	checkVolumeFraction(phi_pure_s,'phi_s')
	checkVolumeFraction(phi_pure_p,'phi_p')
	
	# print phi_pure_s
	#Solvent Chemical Potential in pure solvent phase.
	alpha_0=vh_pure_s/v_r
	if Kier:
		mu_pure_s = (-(1.0/alpha_0)+(1.0/alpha_pure_s)*(0+ln(phi_pure_s))-(1.0/alpha_0)*ln(1-phi_pure_s)-2*(Tsstar/T)*phi_pure_s) #Kier
		# mu_pure_s = alpha_pure_s*(-(1.0/alpha_0)+(1.0/alpha_pure_s)*(1+ln(phi_pure_s))-(1.0/alpha_0)*ln(1-phi_pure_s)-2*(Tsstar/T)*phi_pure_s) #Including 1 and prefactor
		# print mu_pure_s
	elif Hassan or Hassan_Var_Vol:
		mu_pure_s = alpha_pure_s*(-(1.0/alpha_0)+(1.0/alpha_pure_s)*(1+ln(phi_pure_s))-(1.0/alpha_0)*ln(1-phi_pure_s)-2*(Tsstar/T)*phi_pure_s) #Including 1 and prefactor
	elif Condo:
		mu_pure_s = alpha_pure_s*(-1+(1.0/alpha_pure_s)*(1+ln(phi_pure_s))-ln(1-phi_pure_s)-2*(Tsstar/T)*phi_pure_s) #Including 1 and prefactor
		# print 'Kier mu_pure_s is:', mu_pure_s, 'at temperature=', T, 'and phi_pure_s:', phi_pure_s
	phi_mix_p,phi_mix_s=chemicalEquilibrium_bisect(direction,P,T,phi_pure_s,mu_pure_s,Mp,Ms,delta,zeta,Ppstar,Tpstar,Rpstar,Psstar,Tsstar,Rsstar,vh_pure_p,vh_pure_s,v_r,Kier,Hassan,Condo,Hassan_Var_Vol)

	#CHECKING IF MIXTURE VOLUME FRACTION RESULTS ARE VALID.
	checkVolumeFraction([phi_mix_p,phi_mix_s],['phi_p','phi_s'])

	return [P,T,phi_mix_p,phi_mix_s,phi_pure_p,phi_pure_s]

#Bisect Method Original Condo

def EOS_Original(Rtilde,Ptilde,Ttilde,r):

	#Equation of state in Pure solvent phase.
	EOS=Rtilde**2+Ptilde+Ttilde*(ln(1-Rtilde)+(1-1/r)*Rtilde)

	return EOS

def EOS_Original_bisect(direction,Ptilde,Ttilde,r):
	
	step_size=0.2
	range_max=1/step_size

	if direction=='fwd':
		start=0
		end=int(range_max)+1
		step=1
		
	elif direction=='bwd':
		start=int(range_max)
		end=-1
		step=-1

	for i in range(start,end,step):
		k=i*step_size

		
		if k+step*0.00001<0.0:
			guess1=0.00001
		elif k+step*0.00001>1.0:
			guess1=0.99999
		else:
			guess1=k+step*0.00001

		if k+step*0.00001+step*step_size<0.0:
			guess2=0.00001
		elif k+step*0.00001+step*step_size>1.0:
			guess2=0.99999
		else:
			guess2=k+step*0.00001+step*step_size

		Rtilde=0.0
		try:
			Rtilde = bisect(EOS_Original,abs(guess1),abs(guess2),args=(Ptilde,Ttilde,r))
		except:
			# print 'No value found between', k, 'to', k+step_size
			pass
		if Rtilde!=0.0:
			# print 'Hurry! Rtilde is:', Rtilde
			break
	
	if Rtilde==0.0:
		print 'Program Failed to get value of Rtilde in', direction, 'direction'

	return Rtilde

def chemicalEquilibrium_Original(phi_s,direction,mus0,P,T,Psstar,Tsstar,Rsstar,Ppstar,Tpstar,Rpstar,rs,rp,Xsp):

	phi_p=1-phi_s
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	vhm=phi_s*vhs+phi_p*vhp
	Tstar=phi_s*Tsstar+phi_p*Tpstar-phi_s*phi_p*T*Xsp
	Pstar=kB*Tstar/vhm

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	r=1/(phi_s/rs+phi_p/rp)
	Rtilde=EOS_Original_bisect(direction,Ptilde,Ttilde,r)
	vtilde=1/Rtilde
	Pstilde=P/Psstar
	Tstilde=T/Tsstar

	#Equilibrium Condition is:
	musm=mus0

	#Solvent chemical potential in mixture phase.
	residual_mu= -musm+(ln(phi_s)+(1-rs/rp)*phi_p+rs*Rtilde*Xsp*phi_p**2+rs*((-1*Rtilde+Pstilde*vtilde)/Tstilde+(vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/rs))

	# ##################Falto Trying Start
	# phi_s_falto=phi_s
	# phi_p_falto=1-phi_s_falto
	# vhm_falto=phi_s_falto*vhs+phi_p_falto*vhp
	# Tstar_falto=phi_s_falto*Tsstar+phi_p_falto*Tpstar-phi_s_falto*phi_p_falto*T*Xsp
	# Pstar_falto=kB*Tstar_falto/vhm_falto

	# Ptilde_falto=P/Pstar_falto
	# Ttilde_falto=T/Tstar_falto
	# r_falto=1/(phi_s_falto/rs+phi_p_falto/rp)
	# Rtilde_falto=EOS_Original_bisect(direction,Ptilde_falto,Ttilde_falto,r_falto)
	# vtilde_falto=1/Rtilde_falto
	# Pstilde=P/Psstar
	# Tstilde=T/Tsstar

	# #Equilibrium Condition is:
	# # musm=mus0
	# musm=0
	# # kphi_p_falto=0.81545345
	# # Rtilde_falto=0.909237527009
	# # phi_p_falto=kphi_p_falto/Rtilde_falto
	# # phi_s_falto=1.0-phi_p_falto
	# # # r_falto=1/(phi_s_falto/rs+phi_p_falto/rp)
	# # # Rtilde_falto=EOS_Original_bisect(direction,Ptilde,Ttilde,r_falto)
	# # vtilde_falto=1/Rtilde_falto
	# cesidual_mu_falto= -musm+(ln(phi_s_falto)+(1-rs/rp)*phi_p_falto+rs*Rtilde_falto*Xsp*phi_p_falto**2+rs*((-1*Rtilde_falto+Pstilde*vtilde_falto)/Tstilde+(vtilde_falto-1)*ln(1-Rtilde_falto)+ln(Rtilde_falto)/rs))
	# kphi_p_falto=phi_p_falto*Rtilde_falto
	# kphi_s_falto=phi_s_falto*Rtilde_falto
	# Tps_star=zeta*math.sqrt(Tpstar*Tsstar)
	# alpha_0=vhm_falto/vhs
	# print alpha_0
	# residual_mu_falto= -musm+rs*(-1/alpha_0+(1.0/rs)*(1+ln(kphi_s_falto))-(1/alpha_0)*ln(1-kphi_p_falto-kphi_s_falto)-2*(Tsstar/T)*kphi_s_falto-2*(Tps_star/T)*kphi_p_falto) #Including all factors
	# # residual_mu_falto= -musm+rs*(-1+(1.0/rs)*(0+ln(kphi_s_falto))-ln(1-kphi_p_falto-kphi_s_falto)-2*(Tsstar/T)*kphi_s_falto-2*(Tps_star/T)*kphi_p_falto) #Kier

	# print 'cesidual is', cesidual_mu_falto, 'at phi_p_falto=', kphi_p_falto, 'at phi_s_falto=', kphi_s_falto, 'at Rtilde=', Rtilde_falto,'at T', T, 'and P',P
	# print 'residual is', residual_mu_falto, 'at phi_p_falto=', kphi_p_falto, 'at phi_s_falto=', kphi_s_falto, 'at Rtilde=', Rtilde_falto,'at T', T, 'and P',P

	######################Falto Trying Ends



	if Rtilde==0.0:
		# print 'Unfortunately Rtilde can be be calulated at this guess phi_s:', phi_s
		residual_mu=-9000000000000000
	
	# print 'residual_mu is:', residual_mu
	# print 'phi_s for above residual is:', phi_s
	return residual_mu

def chemicalEquilibrium_Original_bisect(direction,mus0,P,T,Psstar,Tsstar,Rsstar,Ppstar,Tpstar,Rpstar,rs,rp,Xsp):
	
	step_size=0.2
	range_max=1/step_size

	if direction=='fwd':
		start=0
		end=int(range_max)+1
		step=1
		# print 'forward'
		
	elif direction=='bwd':
		start=int(range_max)
		end=-1
		step=-1
		# print 'backward'

	for i in range(start,end,step):
		k=i*step_size
				
		if k+step*0.001<0.0:
			guess1=0.001
		elif k+step*0.001>1.0:
			guess1=0.999
		else:
			guess1=k+step*0.001

		if k+step*0.001+step*step_size<0.0:
			guess2=0.001
		elif k+step*0.001+step*step_size>1.0:
			guess2=0.999
		else:
			guess2=k+step*0.001+step*step_size

		phi_s=0.0
		try:
			phi_s = bisect(chemicalEquilibrium_Original,abs(guess1),abs(guess2),args=(direction,mus0,P,T,Psstar,Tsstar,Rsstar,Ppstar,Tpstar,Rpstar,rs,rp,Xsp))
		except:
			# print 'No value found between', k, 'to', k+step_size
			pass
		if phi_s!=0.0 and phi_s>0.0:
			# print 'Hurry! phi_s is:', phi_s
			break

	if phi_s==0.0 or phi_s<0.0:
		print 'Program Failed to get value of phi_s in', direction, 'direction'
	# print phi_s
	return phi_s

def binaryPhaseEquilibriumCondo_Original_bisect(direction,P,T,Mp,Ms,**kwargs):

	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))
	
	#PURE FLUID PARAMETERS.
	rs = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	rp = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	
	Pstilde=P/Psstar
	Tstilde=T/Tsstar

	# print Psstar
	# print Tsstar
	# print Rsstar
	# print rs
	# print Ms

	Rstilde=EOS_Original_bisect(direction,Pstilde,Tstilde,rs)

	Pptilde=P/Ppstar
	Tptilde=T/Tpstar
	Rptilde=EOS_Original_bisect(direction,Pptilde,Tptilde,rp)
	# print Rstilde
	vstilde=1/Rstilde

	#Solvent chemical potential in mixture phase.
	mus0=rs*((-1*Rstilde+Pstilde*vstilde)/Tstilde+(vstilde-1)*ln(1-Rstilde)+ln(Rstilde)/rs)
	# print 'condo mu_pure_s is:', mus0, 'at temperature=', T, 'and phi_pure_s:', Rstilde

	Tspstar=zeta*math.sqrt(Tsstar*Tpstar)
	Xsp=(Tsstar+Tpstar-2*Tspstar)/T

	# R=Rstilde*Rsstar
	phi_s=chemicalEquilibrium_Original_bisect(direction,mus0,P,T,Psstar,Tsstar,Rsstar,Ppstar,Tpstar,Rpstar,rs,rp,Xsp)

	phi_p=1-phi_s
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	vhm=phi_s*vhs+phi_p*vhp
	Tstar=phi_s*Tsstar+phi_p*Tpstar-phi_s*phi_p*T*Xsp
	Pstar=kB*Tstar/vhm

	Ptilde=P/Pstar
	Ttilde=T/Tstar
	r=1/(phi_s/rs+phi_p/rp)
	Rtilde=EOS_Original_bisect(direction,Ptilde,Ttilde,r)

	return [phi_s,phi_p,Rtilde,Rstilde,Rptilde]

###########################################################################
###########################################################################
#Entropy Functions

def ThermodynamicVariables(P,T,phi_p,phi_s,Rptilde,Mp,Ms,**kwargs):

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)

	vhp = vh_pure_p
	vhs = vh_pure_s

	# vhp = kB*Tpstar/Ppstar
	# vhs = kB*Tsstar/Psstar

	#Removing Rescaling:
	alpha_0_s=vhs/v_r
	alpha_0_p=vhp/v_r
	# alpha_0=vhp/v_r
	Tpstar=Tpstar*alpha_0_p
	Tsstar=Tsstar*alpha_0_s
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	Vratio=1.0
	#Following are My Theory Equations with, in general, v!=v_0:
	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))

	#Below assuming system to be in pure vhp state but phi_s present some how:
	alpha_s = alpha_s0*vhs/vhp
	S_1=(Ppstar/(Rpstar*Tpstar))*(-((1/(phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p0))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	'''
	#To Reintroducing Rescaling:

	Tpstar=Tpstar/alpha_0_p
	Tsstar=Tsstar/alpha_0_s
	alpha_p = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	vhm = delta*v_r
	alpha_0=vhm/v_r
	Vratio=alpha_0
	epsilon_p=epsilon_p/alpha_0

	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))	
	S_2=(Ppstar/(Rpstar*Tpstar))*(-((1/(alpha_0*phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	'''
	#New: For "Hassan_Var_Vol" i.e. variable hole volume
	
	#To Reintroducing Rescaling:

	Tpstar=Tpstar/alpha_0_p
	Tsstar=Tsstar/alpha_0_s
	alpha_p = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	Rtilde=phi_s+phi_p
	vhm=(phi_s/Rtilde)*vhs+(phi_p/Rtilde)*vhp

	# pfactor=Mp/alpha_p
	# sfactor=Ms/alpha_s
	# vhm=(phi_s*sfactor*vhs+phi_p*pfactor*vhp)/(phi_s*sfactor+phi_p*pfactor)						

	alpha_0=vhm/v_r
	Vratio=alpha_0
	epsilon_p=epsilon_p/alpha_0

	F=(g*exp(-Vratio*epsilon_p/(kB*T)))/(1+g*exp(-Vratio*epsilon_p/(kB*T)))	
	S_2=(Ppstar/(Rpstar*Tpstar))*(-((1/(alpha_0*phi_p))*(1-phi_s-phi_p)*(ln(1-phi_s-phi_p)))-((1/(alpha_s))*(phi_s/phi_p)*(ln(phi_s)))-((1/(alpha_p))*(ln(phi_p)))+((F*epsilon_p)/(kB*T))-((1/Vratio)*ln(1-F)))
	
	return [S_1,S_2]

def CondoThermodynamicVariables_Modified(P,T,phi_p,phi_s,Mp,Ms,**kwargs):
	
	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	alpha_p0=(Mp*Ppstar)/(kB*Rpstar*Tpstar)
	alpha_s0=(Ms*Psstar)/(kB*Rsstar*Tsstar)

	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar

	# Tps_star=zeta*math.sqrt(Tpstar*Tsstar)

	Rtilde=phi_s+phi_p
	vhm=(phi_s/Rtilde)*vhs+(phi_p/Rtilde)*vhp
	alpha_pm = alpha_p0
	alpha_sm = alpha_s0					

	r=Rtilde/((phi_p/alpha_pm)+(phi_s/alpha_sm))					

	Fp=(((z-2)*exp(-epsilon_p/(kB*T)))/(1+(z-2)*exp(-epsilon_p/(kB*T))))
	Fs=(((z-2)*exp(-epsilon_s/(kB*T)))/(1+(z-2)*exp(-epsilon_s/(kB*T))))

	S_1=(Ppstar/(Rpstar*Tpstar))*(-((1-Rtilde)*ln(1-Rtilde)/phi_p)-(ln(phi_p)/alpha_pm)+(ln(alpha_pm)/alpha_pm)-(phi_s*ln(phi_s)/phi_p)-(1+((ln(2/z)-1)/r))-(((alpha_pm-2)/alpha_pm)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T))))-((phi_s/phi_p)*(-ln(1-Fs)+(Fs*epsilon_s/(kB*T))+1+((ln(2/z)-1)/r))))

	# -((1-Rtilde)*ln(1-Rtilde)/phi_p)-(ln(phi_p)/alpha_pm)+(ln(alpha_pm)/alpha_pm)-(phi_s*ln(phi_s)/phi_p)-(1+((ln(2/z)-1)/r))-(((alpha_pm-2)/alpha_pm)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T))))-((phi_s/phi_p)*(-ln(1-Fs)+(Fs*epsilon_s/(kB*T))+1+((ln(2/z)-1)/r)))

	S_term1=(Ppstar/(Rpstar*Tpstar))*(-((1-Rtilde)*ln(1-Rtilde)/phi_p))
	S_term2=(Ppstar/(Rpstar*Tpstar))*(-(phi_s*ln(phi_s)/phi_p))
	S_term3=(Ppstar/(Rpstar*Tpstar))*(-(ln(phi_p)/alpha_pm))
	S_term4=(Ppstar/(Rpstar*Tpstar))*(-(((alpha_pm-2)/alpha_pm)*(-(Fp*epsilon_p/(kB*T)))))
	S_term5=(Ppstar/(Rpstar*Tpstar))*(-(((alpha_pm-2)/alpha_pm)*(ln(1-Fp))))
	S_term6=(Ppstar/(Rpstar*Tpstar))*(-((phi_s/phi_p)*(-ln(1-Fs)+(Fs*epsilon_s/(kB*T))+1+((ln(2/z)-1)/r))))
	S_term7=(Ppstar/(Rpstar*Tpstar))*(-(((ln(2/z)-1)/r)))

	vtilde=1/Rtilde
	S_2=(Ppstar/(Rpstar*Tpstar))*(1/(phi_p/Rtilde))*(-1*((vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/r+((phi_s/Rtilde)/alpha_sm)*ln((phi_s/Rtilde)/alpha_sm)+((phi_p/Rtilde)/alpha_pm)*ln((phi_p/Rtilde)/alpha_pm)+1+(ln(2/z)-1)/r+((phi_s/Rtilde)/alpha_sm)*(alpha_sm-2)*(ln(1-Fs)-(Fs*epsilon_s/(kB*T)))+((phi_p/Rtilde)/alpha_pm)*(alpha_pm-2)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T)))))

	return S_1,S_2

def CondoThermodynamicVariables_Original(P,T,Rtilde,phi_s,Mp,Ms,**kwargs):
			
	phi_p=1.0-phi_s
	# print 'Temperature is', T

	for key,value in kwargs.items():
		exec "%s=%s" % (key,value)
	
	#PURE FLUID PARAMETERS.
	# vhp = kB*Tpstar/Ppstar
	# vhs = kB*Tsstar/Psstar
	rp = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	rs = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	
	# vhm=phi_s*vhs+phi_p*vhp
	
	# Tspstar=zeta*math.sqrt(Tsstar*Tpstar)
	# Xsp=(Tsstar+Tpstar-2*Tspstar)/T
	# Tstar=phi_s*Tsstar+phi_p*Tpstar-phi_s*phi_p*T*Xsp
	# Pstar=kB*Tstar/vhm

	# Ptilde=P/Pstar
	# Ttilde=T/Tstar

	r=1/(phi_s/rs+phi_p/rp)

	vtilde=1/Rtilde

	# Pstilde=P/Psstar
	# Tstilde=T/Tsstar	
	# Pptilde=P/Ppstar
	# Tptilde=T/Tpstar

	Fp=((z-2)*exp(-epsilon_p/(kB*T)))/(1+(z-2)*exp(-epsilon_p/(kB*T)))
	Fs=((z-2)*exp(-epsilon_s/(kB*T)))/(1+(z-2)*exp(-epsilon_s/(kB*T)))

	S_1=-1*((vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/r+(phi_s/rs)*ln(phi_s/rs)+(phi_p/rp)*ln(phi_p/rp)+1+(ln(2/z)-1)/r+(phi_s/rs)*(rs-2)*(ln(1-Fs)-(Fs*epsilon_s/(kB*T)))+(phi_p/rp)*(rp-2)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T))))
	S_2=(Ppstar/(Rpstar*Tpstar))*(1/phi_p)*(-1*((vtilde-1)*ln(1-Rtilde)+ln(Rtilde)/r+(phi_s/rs)*ln(phi_s/rs)+(phi_p/rp)*ln(phi_p/rp)+1+(ln(2/z)-1)/r+(phi_s/rs)*(rs-2)*(ln(1-Fs)-(Fs*epsilon_s/(kB*T)))+(phi_p/rp)*(rp-2)*(ln(1-Fp)-(Fp*epsilon_p/(kB*T)))))

	return [S_1,S_2]

###########################################################################

def binarySolubilitySwellingCHV(P,T,Mp,Ms,**kwargs):

	# Xs:the total solubility m of solvent and co-solvent.
	# Sw:the volume swelling.

	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))

	#PURE FLUID PARAMETERS.
	# vhp = kB*Tpstar/Ppstar
	# vhs = kB*Tsstar/Psstar
	vhp = vh_pure_p
	vhs = vh_pure_s
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)

	# CALCULATION OF VOLUME FRACTIONS AT P, T.

	if forward:
		[Pd,Td,phip,phis,phip0,phis0] = binaryPhaseEquilibriumCHV_bisect('fwd',P,T,Mp,Ms,**kwargs)
	if backward:
		[Pd,Td,phip,phis,phip0,phis0] = binaryPhaseEquilibriumCHV_bisect('bwd',P,T,Mp,Ms,**kwargs)

	# [Pd,Td,phip,phis,phip0,phis0] = binaryPhaseEquilibriumCHV(P,T,Mp,Ms,**kwargs)

	#MIXTURE PARAMETERS.
	if Kier or Hassan:
		vhm = delta*v_r
		# alpha_p = alpha_p0*vhp/vhm
		# alpha_s = alpha_s0*vhs/vhm
		alpha_p = alpha_p0
		alpha_s = alpha_s0
	elif Hassan_Var_Vol:
		Rtilde=phis+phip
		vhm=(phis/Rtilde)*vhs+(phip/Rtilde)*vhp
		# alpha_p = alpha_p0*vhp/vhm				#Hassan: v_h is volume of hole in mixture or pure solvent
		# alpha_s = alpha_s0*vhs/vhm	
		alpha_p = alpha_p0
		alpha_s = alpha_s0

		# pfactor=Mp/alpha_p0
		# sfactor=Ms/alpha_s0
		# vhm=(phis*sfactor*vhs+phip*pfactor*vhp)/(phis*sfactor+phip*pfactor)


	#CALCULATION OF SOLVENT SOLUBILITY (MASS FRACTION) AT P, T.
	Xs = (Ms*phis/alpha_s)/(Mp*phip/alpha_p+Ms*phis/alpha_s)		#Kier Formula
	# Xs = (Ms*phis/alpha_s0)/(Mp*phip/alpha_p0+Ms*phis/alpha_s0)		#Condo Formula
	
	#CALCULATION OF VOLUME SWELLING (Vp/Vm) IN THE MIXTURE AT P, T.
	Sw = phip0/phip
	
	#PRINTING OF RESULTS OF SOLUBILITY AND SWELLING.
	print('At P = {}, T = {}, zeta = {}, delta = {};'.format(P,T,zeta,delta))
	print('Xs = {}, Sw = {};'.format(Xs,Sw))
	print('phip = {}, phis = {}, phip0 = {}, phis0 = {};'.format(phip,phis,phip0,phis0))

	Rtilde=phip+phis

	return [P,T,Xs,Sw,phip,phis,Rtilde,phip0,phis0]

def binarySolubilitySwellingCondo_Original(P,T,Mp,Ms,**kwargs):

	# Xs:the total solubility m of solvent and co-solvent.
	# Sw:the volume swelling.

	for key,value in kwargs.iteritems():
		exec('{} = value'.format(key))

	#PURE FLUID PARAMETERS.
	vhp = kB*Tpstar/Ppstar
	vhs = kB*Tsstar/Psstar
	alpha_p0 = (Ppstar*Mp)/(kB*Tpstar*Rpstar)
	alpha_s0 = (Psstar*Ms)/(kB*Tsstar*Rsstar)
	
	# CALCULATION OF VOLUME FRACTIONS AT P, T.

	if forward:
		[cphis,cphip,Rtilde,Rstilde,Rptilde]=binaryPhaseEquilibriumCondo_Original_bisect('fwd',P,T,Mp,Ms,**kwargs)
	if backward:
		[cphis,cphip,Rtilde,Rstilde,Rptilde]=binaryPhaseEquilibriumCondo_Original_bisect('bwd',P,T,Mp,Ms,**kwargs)

	kphis=cphis*Rtilde
	kphip=cphip*Rtilde
	kphip0=Rptilde
	kphis0=Rstilde


	#CALCULATION OF SOLVENT SOLUBILITY (MASS FRACTION) AT P, T.
	Xs = (Ms*kphis/alpha_s0)/(Mp*kphip/alpha_p0+Ms*kphis/alpha_s0)
	# Xs_wrong = (kphis*Rsstar)/(kphis*Rsstar+kphip*Rpstar)		#Not working
	# vhm=cphis*vhs+cphip*vhp
	# alpha_s=alpha_s0*vhs/vhm
	# alpha_p=alpha_p0*vhp/vhm
	# Xs_wrong = (Ms*kphis/alpha_s)/(Mp*kphip/alpha_p+Ms*kphis/alpha_s) #Not working

	#CALCULATION OF VOLUME SWELLING (Vp/Vm) IN THE MIXTURE AT P, T.
	Sw = kphip0/(kphip)
	
	#PRINTING OF RESULTS OF SOLUBILITY AND SWELLING.
	print('At P = {}, T = {}, zeta = {};'.format(P,T,zeta))
	print('Xs = {}, Sw = {};'.format(Xs,Sw))
	print('kphip = {}, kphis = {}, kphip0 = {}, kphis0 = {};'.format(kphip,kphis,kphip0,kphis0))

	return [P,T,Xs,Sw,kphip,kphis,Rtilde,kphip0,kphis0]

#####################################################################################################
#####################################################################################################

def calculateThermodynamicVariables(P0,T0,phi_p,phi_s,phi_p0,phi_s0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		TD = ThermodynamicVariables(P0,T0,phi_p,phi_s,phi_p0,Mp,Ms,**kwargs)
		result = TD
	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(4)]
		T = range(0,len(P0))
		S_1 = range(0,len(P0))
		S_2 = range(0,len(P0))
		
		for i in range(0,len(P0)):
			TD = ThermodynamicVariables(P0[i],T0,phi_p[i],phi_s[i],phi_p0[i],Mp,Ms,**kwargs)
			T[i] = T0
			S_1[i] = TD[0]
			S_2[i] = TD[1]

		result[0] = P0
		result[1] = T
		result[2] = S_1
		result[3] = S_2

	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		S_1 = range(0,len(T0))
		S_2 = range(0,len(T0))

		for i in range(0,len(T0)):
			TD = ThermodynamicVariables(P0,T0[i],phi_p[i],phi_s[i],phi_p0[i],Mp,Ms,**kwargs)
			P[i] = P0
			S_1[i] = TD[0]
			S_2[i] = TD[1]
	
		result[0] = P
		result[1] = T0
		result[2] = S_1
		result[3] = S_2
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		S_1 = range(0,len(T0))
		S_2 = range(0,len(T0))

		for i in range(0,len(T0)):
			TD = ThermodynamicVariables(P0[i],T0[i],phi_p[i],phi_s[i],phi_p0[i],Mp,Ms,**kwargs)
			S_1[i] = TD[0]
			S_2[i] = TD[1]

		result[0] = P0
		result[1] = T0
		result[2] = S_1
		result[3] = S_2
	
	else:
		raise ValueError('In calculateThermodynamicVariables: Unknown error involving P0 and T0.')

	return result

def calculateCondoThermodynamicVariables_Original(P0,T0,kphi_p,kphi_s,kphi_p0,kphi_s0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		Rtilde=kphi_p+kphi_s
		cphi_s=kphi_s/Rtilde
		TD = CondoThermodynamicVariables_Original(P0,T0,Rtilde,cphi_s,Mp,Ms,**kwargs)
		result = TD
	
	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(4)]
		T = range(0,len(P0))
		S_1 = range(0,len(P0))
		S_2 = range(0,len(P0))
		Rtilde = range(0,len(P0))
		cphi_s = range(0,len(P0))
	
		for i in range(0,len(P0)):

			Rtilde[i]=kphi_p[i]+kphi_s[i]
			cphi_s[i]=kphi_s[i]/Rtilde[i]
			TD = CondoThermodynamicVariables_Original(P0[i],T0,Rtilde[i],cphi_s[i],Mp,Ms,**kwargs)
			T[i] = T0
			S_1[i] = TD[0]
			S_2[i] = TD[1]

		result[0] = P0
		result[1] = T
		result[2] = S_1
		result[3] = S_2

	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		S_1 = range(0,len(T0))
		S_2 = range(0,len(T0))
		Rtilde = range(0,len(T0))
		cphi_s = range(0,len(T0))
		# print Rtilde
		# print T0

		for i in range(0,len(T0)):
			Rtilde[i]=kphi_p[i]+kphi_s[i]
			cphi_s[i]=kphi_s[i]/Rtilde[i]
			TD = CondoThermodynamicVariables_Original(P0,T0[i],Rtilde[i],cphi_s[i],Mp,Ms,**kwargs)
			P[i] = P0
			S_1[i] = TD[0]
			S_2[i] = TD[1]
	
		result[0] = P
		result[1] = T0
		result[2] = S_1
		result[3] = S_2
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(4)]
		P = range(0,len(T0))
		S_1 = range(0,len(T0))
		S_2 = range(0,len(T0))
		Rtilde = range(0,len(T0))
		cphi_s = range(0,len(T0))

		for i in range(0,len(T0)):
			Rtilde[i]=kphi_p[i]+kphi_s[i]
			cphi_s[i]=kphi_s[i]/Rtilde[i]
			TD = CondoThermodynamicVariables_Original(P0[i],T0[i],Rtilde[i],cphi_s[i],Mp,Ms,**kwargs)
			S_1[i] = TD[0]
			S_2[i] = TD[1]

		result[0] = P0
		result[1] = T0
		result[2] = S_1
		result[3] = S_2
	
	else:
		raise ValueError('In calculateCondoThermodynamicVariables_Original: Unknown error involving P0 and T0.')

	return result

def calculateBinarySolubilitySwelling(theory,P0,T0,Mp,Ms,**kwargs):

	if not isListOrNpyArray(P0) and not isListOrNpyArray(T0):
		exec "XSw = binarySolubilitySwelling%s(P0,T0,Mp,Ms,**kwargs)" % (theory)
		result = XSw
	
	elif not isListOrNpyArray(T0) and isListOrNpyArray(P0):
		result = [[range(0,len(P0))] for x in range(9)]
		T = range(0,len(P0))
		m_s = range(0,len(P0))
		Sw = range(0,len(P0))
		phip = range(0,len(P0))
		phis = range(0,len(P0))
		Rtilde = range(0,len(P0))
		phip0 = range(0,len(P0))
		phis0 = range(0,len(P0))
		
		for i in range(0,len(P0)):
			exec "XSw = binarySolubilitySwelling%s(P0[i],T0,Mp,Ms,**kwargs)" % (theory)
			T[i] = XSw[1]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
			phip[i]=XSw[4]
			phis[i]=XSw[5]
			Rtilde[i]=XSw[6]
			phip0[i]=XSw[7]
			phis0[i]=XSw[8]

		result[0] = P0
		result[1] = T
		result[2] = m_s
		result[3] = Sw
		result[4]=phip
		result[5]=phis
		result[6]=Rtilde
		result[7]=phip0
		result[8]=phis0

	elif not isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(9)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		Sw = range(0,len(T0))
		phip = range(0,len(T0))
		phis = range(0,len(T0))
		Rtilde = range(0,len(T0))
		phip0 = range(0,len(T0))
		phis0 = range(0,len(T0))

		for i in range(0,len(T0)):
			exec "XSw = binarySolubilitySwelling%s(P0,T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = XSw[0]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
			phip[i]=XSw[4]
			phis[i]=XSw[5]
			Rtilde[i]=XSw[6]
			phip0[i]=XSw[7]
			phis0[i]=XSw[8]
	
		result[0] = P
		result[1] = T0
		result[2] = m_s
		result[3] = Sw
		result[4]=phip
		result[5]=phis
		result[6]=Rtilde
		result[7]=phip0
		result[8]=phis0
	
	elif isListOrNpyArray(P0) and isListOrNpyArray(T0):
		result = [[range(0,len(T0))] for x in range(9)]
		P = range(0,len(T0))
		m_s = range(0,len(T0))
		Sw = range(0,len(T0))
		phip = range(0,len(T0))
		phis = range(0,len(T0))
		Rtilde = range(0,len(T0))
		phip0 = range(0,len(T0))
		phis0 = range(0,len(T0))

		for i in range(0,len(T0)):
			exec "XSw = binarySolubilitySwelling%s(P0[i],T0[i],Mp,Ms,**kwargs)" % (theory)
			P[i] = XSw[0]
			m_s[i] = XSw[2]
			Sw[i] = XSw[3]
			phip[i]=XSw[4]
			phis[i]=XSw[5]
			Rtilde[i]=XSw[6]
			phip0[i]=XSw[7]
			phis0[i]=XSw[8]

		result[0] = P0
		result[1] = T0
		result[2] = m_s
		result[3] = Sw
		result[4]=phip
		result[5]=phis
		result[6]=Rtilde
		result[7]=phip0
		result[8]=phis0
	
	else:
		raise ValueError('In calculateBinarySwelling: Unknown error involving P0 and T0.')
	
	return result

#####################################################################################################
#####################################################################################################

def rescaling(Tpstar,Tsstar,vh_pure_p,vh_pure_s,v_r):

	alpha_0_s=vh_pure_s/v_r
	Tsstar=Tsstar/alpha_0_s

	alpha_0_p=vh_pure_p/v_r
	Tpstar=Tpstar/alpha_0_p

	return Tpstar,Tsstar
# Author: Kier von Konigslow
# Date: November 2017
#
import os,sys,math,csv,matplotlib.pyplot as plt,numpy as npy
# from p_params import *
# from s_params import *
from loadExperimentalData import *
lib_path = os.path.abspath(os.path.join('..'))
sys.path.append(lib_path)
from loadPhysicalConstants import *
from calculatePureVariables import calculateNewMolecularParameters
# from wrapperFunctions import calculateBinarySolubilitySwelling
# from calculateBinaryResidual import calculateBinarySSQ
from checkResults import *
from sympy import *
import warnings
import cmath
from scipy.optimize import bisect,fsolve
from scipy.interpolate import interp1d
from Parameters_of_Different_Polymers import *
from All_Functions import *
from Parameters_for_Mixtures_and_Tg import *

###########################################################################
###########################################################################

def discard_zeros(x,y):
	
	for i in range(len(x)):
		if x[i]==0:
			y[i]=0

	for i in range(len(y)):
		if y[i]==0:
			x[i]=0

	# print x
	# print y

	x = npy.delete(x, npy.argwhere( (x >= 0) & (x <= 0) ))
	y = npy.delete(y, npy.argwhere( (y >= 0) & (y <= 0) ))

	# print x
	# print y
				
	return x,y

#Setting font size
axis_size = 20
title_size = 20
size = 16
label_size = 16
plt.rcParams['xtick.labelsize'] = label_size
plt.rcParams['ytick.labelsize'] = label_size

#Linestyles
lsT1K = '-'
lsT2K = '--'
lsT3K = ':'

#General line properties.
linewidth = 2
markersize = 8

x = npy.array([1,0,3,4,5,6,7,8,9,10])
y = npy.array([1,2,9,16,0,36,49,64,81,100])

print x,y
x,y=discard_zeros(x,y)
print x,y
plt.plot(x,y,'r',marker='x',ls=lsT1K)

P=[0.101325,0.75673333,1.41214167,2.06755,2.72295833]
Tg_bisect_fwd=[372.0703125,362.6953125,252.1484375,270.5078125,292.7734375]
Pg_bisect_fwd=[1.33789062,1.92382812,2.48242188,2.84179688,2.84960938,2.42382812,1.58398438]
T=[250.0,266.66666667,283.33333333,300.0,316.66666667,333.33333333,350.0,]


#Show plot windows.
plt.show()
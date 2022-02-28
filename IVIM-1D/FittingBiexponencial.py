# Author: João Antônio
# E-mail: joao.ferreira@estudante.ufscar.br

"""Biexponential fitting of a biexponential decay signal. The signal was 
simulated using parameters found on literature for healthy and pathologic 
regions of a human brain. The fitting process was made using the
function curve_fit from the library scipy.optmize.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#Simulation Parameters
bvalueMax=1000
Deltabvalue=2
bvalue= np.arange(0,bvalueMax,Deltabvalue)
SignalDecay=[]


#Select Parameters
Combination=6
if Combination==1:
    PerfusionFraction=9.0E-2
    DiffusionCoefficient=7.0E-4
    PseudoDiffusionCoefficient=1.0E-2
elif Combination==2:
    PerfusionFraction=3.0E-2
    DiffusionCoefficient=9.0E-4
    PseudoDiffusionCoefficient=9.0E-2
elif Combination==3:
    PerfusionFraction=1.25E-1
    DiffusionCoefficient=7.5E-4
    PseudoDiffusionCoefficient=1.0E-2
elif Combination==4:
    PerfusionFraction=5.0E-2
    DiffusionCoefficient=8.0E-4
    PseudoDiffusionCoefficient=9.0E-2
elif Combination==5:
    PerfusionFraction=8.0E-2
    DiffusionCoefficient=2.0E-3
    PseudoDiffusionCoefficient=6.0E-3
elif Combination==6:
    PerfusionFraction=1.8E-1
    DiffusionCoefficient=1.0E-3
    PseudoDiffusionCoefficient=3.5E-2
else:
    print('Revisar valor da combinação')

    
#Biexponential Function
def Biexponential(x,a,b,c):
    return a*np.exp(-1*b*x) + (1-a)*np.exp(-1*c*x)

    
#Signal Simulation    
for i in enumerate(bvalue):
    Signal=Biexponential(i[1],PerfusionFraction,DiffusionCoefficient,PseudoDiffusionCoefficient)
    SignalDecay.append(Signal)

#Signal Fitting
Parameters, Covariancie=curve_fit(Biexponential,bvalue,SignalDecay,p0=[0.0, 0.0, 0.0],bounds=(-np.inf, np.inf),method='dogbox')

SimulatedPerfusionFraction=Parameters[0]
SimulatedDiffusionCoefficient=Parameters[1]
SimulatedPseudoDiffusionCoefficient=Parameters[2]

BiexponentialFitting=Biexponential(bvalue,*Parameters)
Residue=SignalDecay-BiexponentialFitting

#Graphics Generation
fig1, axes =plt.subplots(nrows=2,ncols=2,figsize=(12,8))
axes[0,0].set_title('Simulated Signal')
axes[0,0].set_xlabel('bvalue[mm^2/s]')
axes[0,0].set_ylabel('Magnetization [u.a.]')
axes[0,0].plot(bvalue,SignalDecay,color="black")
axes[0,1].set_title('Estimated Signal')
axes[0,1].set_xlabel('bvalue[mm^2/s]')
axes[0,1].set_ylabel('Magnetization [u.a.]')
axes[0,1].plot(bvalue,BiexponentialFitting,color="red")
axes[1,0].set_title('Comparison')
axes[1,0].set_xlabel('bvalue[mm^2/s]')
axes[1,0].set_ylabel('Magnetization [u.a.]')
axes[1,0].plot(bvalue,SignalDecay,color="black")
axes[1,0].plot(bvalue,BiexponentialFitting,color="red")
plt.legend(['Simulated Signal','Estimated Signal'])
axes[1,1].set_title('Residue')
axes[1,1].set_xlabel('bvalue[mm^2/s]')
axes[1,1].set_ylabel('Magnetization [u.a.]')
axes[1,1].plot(bvalue,Residue,color="blue")
fig1.tight_layout()

#Printing Parameters
print(SimulatedPerfusionFraction)
print(SimulatedDiffusionCoefficient)
print(SimulatedPseudoDiffusionCoefficient)

################################################

'''plt.title("Double Exponential Decay Signal vs b-value")
plt.xlabel('bvalue[mm^2/s]')
plt.ylabel('Magnetização [u.a.]')
plt.plot(bvalue,SignalDecay,color="black")
plt.show()'''

'''plt.title("Ajuste do Sinal de Decaimento vs b-value")
plt.xlabel('bvalue[mm^2/s]')
plt.ylabel('Magnetização [u.a.]')
plt.plot(bvalue,FittingBiexponencial,color="red")
plt.show()'''

'''plt.title("Resíduo vs b-value")
plt.xlabel('bvalue[mm^2/s]')
plt.ylabel('Magnetização [u.a.]')
plt.plot(bvalue,Residuo,color="blue")
plt.show()'''
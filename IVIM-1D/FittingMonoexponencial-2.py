# Author: João Antônio
# E-mail: joao.ferreira@estudante.ufscar.br

"""Segmented fitting of a biexponential decay signal. While a section of 
the signal is fitted as a monoexponential curve, the remainig is fitted as a 
biexponenctial. The signal was simulated using parameters found on literature 
for healthy and pathologic regions of a human brain. The fitting process was 
made using the function curve_fit from the library scipy.optmize.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


#Simulation Parameters
bvalueMax=900
Deltabvalue=2
bvalue=np.arange(0,bvalueMax,Deltabvalue)
SignalDecay=[]

#Select Parameters
Combination=6
if Combination==1:
    PerfusionFractionComplement=9.1E-1
    PerfusionFraction=9.0E-2
    DiffusionCoefficient=7.0E-4
    PseudoDiffusionCoefficient=1.0E-2
    bvalueCut=400
elif Combination==2:
    PerfusionFractionComplement=9.7E-1
    PerfusionFraction=3.0E-2
    DiffusionCoefficient=9.0E-4
    PseudoDiffusionCoefficient=9.0E-2
    bvalueCut=200
elif Combination==3:
    PerfusionFractionComplement=8.75E-1
    PerfusionFraction=1.25E-1
    DiffusionCoefficient=7.5E-4
    PseudoDiffusionCoefficient=1.0E-2
    bvalueCut=400
elif Combination==4:
    PerfusionFractionComplement=9.5E-1
    PerfusionFraction=5.0E-2
    DiffusionCoefficient=8.0E-4
    PseudoDiffusionCoefficient=9.0E-2
    bvalueCut=200
elif Combination==5:
    PerfusionFractionComplement=9.2E-1
    PerfusionFraction=8.0E-2
    DiffusionCoefficient=2.0E-3
    PseudoDiffusionCoefficient=6.0E-3
    bvalueCut=425
elif Combination==6:
    PerfusionFractionComplement=8.2E-1
    PerfusionFraction=1.8E-1
    DiffusionCoefficient=1.0E-3
    PseudoDiffusionCoefficient=3.5E-2
    bvalueCut=200
else:
    print('Review tehe bvalueCut used')


#Biexponential and exponential function
def biexponential(x,a,b,c,d):
    return a*np.exp(-1*b*x) + c*np.exp(-1*d*x)

def exponential (x,a,b):
    return a*np.exp(-1*b*x)


#Signal Simulation
for i in enumerate(bvalue):
    Signal= biexponential(i[1],PerfusionFractionComplement,DiffusionCoefficient,PerfusionFraction,PseudoDiffusionCoefficient)
    SignalDecay.append(Signal)

#Signal Segmentation
bvalueEnd=bvalue[slice(bvalueCut,450,1)]
SignalDecayEnd=SignalDecay[slice(bvalueCut,450,1)]
bvalueBeginning=bvalue[slice(0,bvalueCut,1)]
SignalDecayBeginning=SignalDecay[slice(0,bvalueCut,1)]

#Signal fitting
Parameters, Covariancie= curve_fit(exponential,bvalueEnd,SignalDecayEnd,p0=[0.0,0.0],bounds=(-np.inf, np.inf),method='dogbox')

SimulatedPerfusionFractionComplement=Parameters[0]
SimulatedDiffusionCoefficient=Parameters[1]

ExponentialFitting=exponential(bvalue,*Parameters)


#Testing the bvalueCut
DeltaMin=1.0E-5
Delta=abs((ExponentialFitting[bvalueCut]-SignalDecay[bvalueCut])/ExponentialFitting[bvalueCut])
if Delta<DeltaMin:
    print("The bvalueCut used is valid")
else:
    print("The bvalueCut used isn't valid")


#Biexponential function
def biexponential2(x,c,d):
    
    '''In this  case, the biexponential function takes the values acquired 
    in the previous fittig as its parameters'''
    
    return (SimulatedPerfusionFractionComplement)*np.exp(-1*SimulatedDiffusionCoefficient*x) + c*np.exp(-1*d*x)


#Signal fitting
Parameters2, Covariancie2= curve_fit(biexponential2,bvalueBeginning,SignalDecayBeginning,p0=[0.0,0.0],bounds=(-np.inf, np.inf))

SimulatedPerfusionFraction=Parameters2[0]
SimulatedPseudoDiffusionCoefficient=Parameters2[1]

BiexponentialFitting=biexponential2(bvalue,*Parameters2)
Residue=SignalDecay-BiexponentialFitting


#Graphic Generation
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


#Printing Parametters
print(SimulatedPerfusionFraction)
print(SimulatedDiffusionCoefficient)
print(SimulatedPseudoDiffusionCoefficient)

##############################################################################
''''plt.title("Sinal de Decaimento Biexponencial Final vs b-value")
plt.xlabel('bvalue[mm^2/s]')
plt.ylabel('Magnetização [u.a.]')
plt.plot(bvalueFim,SinalDecaimentoFim,color="black")
plt.show()'''''''''

'''plt.title('Sinal Difusão - Estimado')
plt.xlabel('bvalue[mm^2/s]')
plt.ylabel('Magnetização [u.a.]')
plt.plot(bvalue,FittingExponencial,color="black")
plt.show()'''




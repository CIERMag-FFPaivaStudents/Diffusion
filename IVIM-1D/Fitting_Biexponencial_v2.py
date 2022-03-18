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


def SelectParameters(CombinationNumber):
    
    '''
    Function to select the coefficients that are going to be used on the simulation.
    
    Parameter
    ---------
    CombinationNumber: Float64
        Select the combination of coeffcients that are going to be used.
    
    
    Returns
    ----------
    PFC: Float 64
        Value of the Perfusion Fraction Complement
    PF: Float 64
        Value of the Perfusion Fraction
    DC: Float 64
        Value of the Diffusion Coefficient
    PDC: Float 64
        Value of the Pseudodiffusion Coefficient  
        
    '''
    
    if CombinationNumber==1:
        PF=9.0E-2
        DC=7.0E-4
        PDC=1.0E-2
    elif CombinationNumber==2:
        PF=3.0E-2
        DC=9.0E-4
        PDC=9.0E-2
    elif CombinationNumber==3:
        PF=1.25E-1
        DC=7.5E-4
        PDC=1.0E-2
    elif CombinationNumber==4:
        PF=5.0E-2
        DC=8.0E-4
        PDC=9.0E-2
    elif CombinationNumber==5:
        PF=8.0E-2
        DC=2.0E-3
        PDC=6.0E-3
    elif CombinationNumber==6:
        PF=1.0E-1
        DC=7.0E-4
        PDC=5.0E-2
    else:
        print('Review the bvalueCut used')
        
    return [PF,DC,PDC]


def biexponential(x,a,b,c):
    
    '''Function to simulate a biexponential decay.
    
    Parameter
    ---------
    x: int32
        Integer containig abscissa information
    a: Float64
        Float containing the extent of the first exponential
    b: Float64
        Float containing the coefficient of the first exponential
    c: Float64
        Float containing the extent of the second exponential
    
    Returns
    ----------
    biexponential: Float64  
        Float containing the value of the biexponential decay signal at
        a given point
    '''
    
    biexponential=a*np.exp(-1*b*x) + (1-a)*np.exp(-1*c*x)
    
    return biexponential


def SignalSimulation(X,PF,DC,PDC):
    
    ''' Function to simulate the biexponential decay signal for given values
    of the IVIM-MRI parameters.
    
    Parameters
    ----------
    X: Array of int32
        Array containing the abscissa's values that are going to be used on the
        simulation
    PFC: Float 64
        Value of the Perfusion Fraction Complement
    PF: Float 64
        Value of the Perfusion Fraction
    DC: Float 64
        Value of the Diffusion Coefficient
    PDC: Float 64
        Value of the Pseudodiffusion Coefficient
        
    Returns
    -------
    YSimulated: List
        List containing the signal's values that were simulated
    '''
    YSimulated=[]
    for i in enumerate(X):
        Y= biexponential(i[1],PF,DC,PDC)
        YSimulated.append(Y)
        
    return YSimulated

    

def GraphicGeneration(X,YSimulated,YEstimated,Rest):
    
    '''Function to generate a graphic of the simulated signal, the estimated 
    signal, a comparisson of both of them and their residue.
    
    Parameters
    ---------
    X: Array of int32
        Array containing the abscissa's values for all graphics
    YSimulated: List
        List containing the values of the estimated signal
    YEstimated: Array of Float64
        Array containing the values of the simulated signal
    Residue: Array of float64
        Array containing the difference between the YSimulated and the
        YEstimated
    '''
    
    fig1, axes =plt.subplots(nrows=2,ncols=2,figsize=(12,8))
    axes[0,0].set_title('Simulated Signal')
    axes[0,0].set_xlabel('bvalue[mm^2/s]')
    axes[0,0].set_ylabel('Magnetization [u.a.]')
    axes[0,0].plot(X,YSimulated,color="black")
    axes[0,1].set_title('Estimated Signal')
    axes[0,1].set_xlabel('bvalue[mm^2/s]')
    axes[0,1].set_ylabel('Magnetization [u.a.]')
    axes[0,1].plot(X,YEstimated,color="red")
    axes[1,0].set_title('Comparison')
    axes[1,0].set_xlabel('bvalue[mm^2/s]')
    axes[1,0].set_ylabel('Magnetization [u.a.]')
    axes[1,0].plot(X,YSimulated,color="black")
    axes[1,0].plot(X,YEstimated,color="red")
    plt.legend(['Simulated Signal','Estimated Signal'])
    axes[1,1].set_title('Residue')
    axes[1,1].set_xlabel('bvalue[mm^2/s]')
    axes[1,1].set_ylabel('Magnetization [u.a.]')
    axes[1,1].plot(X,Rest,color="blue")
    fig1.tight_layout()

if __name__=='__main__':

    SimulatedParameters=SelectParameters(6)
    bvalueMax=900
    Deltabvalue=2
    bvalue=np.arange(0,bvalueMax,Deltabvalue)
    SignalDecay=SignalSimulation(bvalue,*SimulatedParameters)
    
    EstimatedParameters, Covariancie= curve_fit(biexponential,bvalue,SignalDecay,p0=[0.0,0.0,0.0],bounds=(-np.inf, np.inf),method='dogbox')
    EstimatedPerfusionFraction=EstimatedParameters[0]
    EstimatedDiffusionCoefficient=EstimatedParameters[1]
    EstimatedPseudoDiffusionCoefficient=EstimatedParameters[2]
    
    BiexponentialFitting=biexponential(bvalue,*EstimatedParameters)
    Residue=SignalDecay-BiexponentialFitting
    
    GraphicGeneration(bvalue,SignalDecay,BiexponentialFitting,Residue)
    
    print(EstimatedPerfusionFraction)
    print(EstimatedDiffusionCoefficient)
    print( EstimatedPseudoDiffusionCoefficient)

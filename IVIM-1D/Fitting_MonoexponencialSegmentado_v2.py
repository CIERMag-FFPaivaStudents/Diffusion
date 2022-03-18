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
        PFC=9.1E-1
        PF=9.0E-2
        DC=7.0E-4
        PDC=1.0E-2
    elif CombinationNumber==2:
        PFC=9.7E-1
        PF=3.0E-2
        DC=9.0E-4
        PDC=9.0E-2
    elif CombinationNumber==3:
        PFC=8.75E-1
        PF=1.25E-1
        DC=7.5E-4
        PDC=1.0E-2
    elif CombinationNumber==4:
        PFC=9.5E-1
        PF=5.0E-2
        DC=8.0E-4
        PDC=9.0E-2
    elif CombinationNumber==5:
        PFC=9.2E-1
        PF=8.0E-2
        DC=2.0E-3
        PDC=6.0E-3
    elif CombinationNumber==6:
        PFC=9.0E-1
        PF=1.0E-1
        DC=7.0E-4
        PDC=5.0E-2    
    else:
        print('Review the bvalueCut used')
        
    return [PFC,PF,DC,PDC]


def biexponential(x,a,b,c,d):
    
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
    d: Float64
        Float containing the coefficient of the second exponential
    
    Returns
    ----------
    biexponential: Float64  
        Float containing the value of the biexponential decay signal at
        a given point
    '''
    
    biexponential=a*np.exp(-1*b*x) + c*np.exp(-1*d*x)
    
    return biexponential


def exponential (x,a,b):

    '''Function to simulate a exponential decay.
    
    Parameter
    ---------
    x: int32
        Integer containig abscissa information
    a: Float64
        Float containing the extent of the exponential
    b: Float64
        Float containing the coefficient of the exponential
   
    
    Returns
    ----------
    exponential: Float64  
        Float containing the value of the exponential decay signal at
        a given point
    '''
    exponential=a*np.exp(-1*b*x)
    
    return exponential

def biexponential2(x,c,d):
    
    '''Function to simulate a biexponential decay. In this  case, the biexponential
    function takes the values acquired in the previous fittig as its parameters.
    
    Parameter
    ---------
    x: int32
        Integer containig abscissa's information
    c: Float64
        Float containing the extent of the second exponential
    d: Float64
        Float containing the coefficient of the second exponential
    
    Returns
    ----------
    biexponential: Float64  
        Float containing the value of the biexponential decay signal at
        a given point
    '''
    biexponential2=(EstimatedPerfusionFractionComplement)*np.exp(-1*EstimatedDiffusionCoefficient*x) + c*np.exp(-1*d*x)
    return biexponential2


def SignalSimulation(X,PFC,PF,DC,PDC):
    
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
        Y= biexponential(i[1],PFC,DC,PF,PDC)
        YSimulated.append(Y)
        
    return YSimulated


def SignalSegmentation(X,YSimulated,Cut):
    
    """ Function to split the X and the YEstimated in two slices
    according to the Cut point.
    
    Parameters
    ----------
    X: Array of int32
        Array containing the abscissa's values that are going to be segmented
    YEstimatedy: List
        List containing the signal's values that are going to be segmented
    Cut: int
        Integer that informs the point in X where the signal will be splited
        
    Returns
    -------
    XEnd: Array of int32
        Array containing the points of X that are superior than Cut
    YEnd: List
        List containing the points of the YEstimated which correspond to a 
        point of X superior than Cut
        
    """
    XEnd=X[slice(Cut,450,1)]
    YEnd=YSimulated[slice(Cut,450,1)]
    return [XEnd,YEnd]

def SegmentationVerification(XEnd,YEnd):
    
    '''Function that verifies if the function generated in the segmentation
    process is a exponential.
    
    Parameters
    ----------
    XEnd: Array of int32
        Array containing the points of X that are superior than Cut
    YEnd: List
        List containing the points of the YEstimated which correspond to a 
        point of X superior than Cut
        
    Returns
    -------
    Inform if the Cut used is valid or not
    '''
    
    LinearizedSignal=np.log(YEnd)
    Residual=abs(np.polyfit(XEnd,LinearizedSignal,1,full=True)[1])
    if Residual<1.0E-5:
        print("The bvalueCut used is valid")
    else:
        print("The bvalueCut used is not valid")
    return Residual

def GraphicGeneration(X,YSimulated,YEstimated,Residue):
    
    '''Function to generate a graphic of the simulated signal, the estimated 
    signal, a comparisson of both of them and their residue.
    
    Parameters
    ---------
    bvalue: Array of int32
        Array containing the abscissa's values for all graphics
    SignalDecay: List
        List containing the values of the estimated signal
    BiexponentialFitting: Array of Float64
        Array containing the values of the simulated signal
    Residue: Array of float64
        Array containing the difference between the SignalDecay and the
        BiexponentialFitting
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
    axes[1,1].plot(bvalue,Residue,color="blue")
    fig1.tight_layout()

if __name__=='__main__':

    SimulatedParameters=SelectParameters(6)
    bvalueMax=900
    Deltabvalue=2
    bvalue=np.arange(0,bvalueMax,Deltabvalue)
    SignalDecay=SignalSimulation(bvalue,*SimulatedParameters)

    bvalueCut=100
    bvalueEnd=SignalSegmentation(bvalue,SignalDecay,bvalueCut)[0]
    SignalDecayEnd=SignalSegmentation(bvalue,SignalDecay,bvalueCut)[1]
    
    
    SegmentationVerification(bvalueEnd,SignalDecayEnd)
    print(SegmentationVerification(bvalueEnd,SignalDecayEnd))
    
    EstimatedParameters, Covariancie= curve_fit(exponential,bvalueEnd,SignalDecayEnd,p0=[0.0,0.0],bounds=(-np.inf, np.inf),method='dogbox')
    EstimatedPerfusionFractionComplement=EstimatedParameters[0]
    EstimatedDiffusionCoefficient=EstimatedParameters[1]

    EstimatedParameters2, Covariancie2= curve_fit(biexponential2,bvalue,SignalDecay,p0=[0.0,0.0],bounds=(-np.inf, np.inf))
    EstimatedPerfusionFraction=EstimatedParameters2[0]
    EstimatedPseudoDiffusionCoefficient=EstimatedParameters2[1]

    BiexponentialFitting=biexponential2(bvalue,*EstimatedParameters2)
    Residue=SignalDecay-BiexponentialFitting

    GraphicGeneration(bvalue,SignalDecay,BiexponentialFitting,Residue)
    
    print(EstimatedPerfusionFraction)
    print(EstimatedDiffusionCoefficient)
    print( EstimatedPseudoDiffusionCoefficient)
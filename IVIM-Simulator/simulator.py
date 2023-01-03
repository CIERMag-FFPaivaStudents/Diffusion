#Author: Gustavo Solcia
#E-mail: gustavo.solcia@usp.br

"""IVIM simulator with BrainWeb data.

"""

import os
import sys
import numpy as np
from tools.importData import importImage
import matplotlib.pyplot as plt

def diffSpinEcho(TE, TR, field, b, WM_data, GM_data, CSF_data):
    if (field==1.5):
        T1_WM = 600 
        T1_GM = 900
        T1_CSF = 3500

        T2_WM = 80
        T2_GM = 100
        T2_CSF = 2000
    elif (field==3.0):
        T1_WM =830
        T1_GM = 1330
        T1_CSF = 4000

        T2_WM =80
        T2_GM =110
        T2_CSF =2000
    else:
        sys.exit('Wrong field value')

    f_WM = 0.068
    f_GM = 0.10

    Dstar_WM = 25.3e-3
    Dstar_GM = 21.9e-3
    
    D_WM = 0.74e-3
    D_GM = 0.81e-3
    D_CSF = 2e-3 
    
    spinEcho_data = (WM_data*(1-np.exp(-TR/T1_WM))*np.exp(-TE/T2_WM) +
                    GM_data*(1-np.exp(-TR/T1_GM))*np.exp(-TE/T2_GM) +
                    CSF_data*(1-np.exp(-TR/T1_CSF))*np.exp(-TE/T2_CSF))

    diffSpinEcho_data = spinEcho_data*(WM_data*(f_WM*np.exp(-b*Dstar_WM)+(1-f_WM)*np.exp(-b*D_WM))+
                        GM_data*(f_GM*np.exp(-b*Dstar_GM)+(1-f_GM)*np.exp(-b*D_GM))+
                        CSF_data*np.exp(-b*D_CSF))

    return diffSpinEcho_data

def addNoise(mag_data, noise_level):

    shape = np.shape(mag_data)
    num_voxels = shape[0]*shape[1]*shape[2]

    gaussian_noise1 = np.random.normal(0, noise_level*1e-2, num_voxels).reshape(shape)
    gaussian_noise2 = np.random.normal(0, noise_level*1e-2, num_voxels).reshape(shape)
    noise_data = np.sqrt((mag_data+gaussian_noise1)**2/np.sqrt(2) + 
            (mag_data+gaussian_noise2)**2/np.sqrt(2))
    
    return noise_data
if __name__=='__main__':

    TE = float(sys.argv[1]) 
    TR = float(sys.argv[2])
    field = float(sys.argv[3])
    noise = float(sys.argv[4])
    bvals_file = open(str(sys.argv[5]), 'r')

    dataPath = 'data/dataset/'
    WM_filename = os.path.join(dataPath, 't1_brainweb_1mm_BET_pve_2.nii.gz')
    GM_filename = os.path.join(dataPath, 't1_brainweb_1mm_BET_pve_1.nii.gz')
    CSF_filename = os.path.join(dataPath, 't1_brainweb_1mm_BET_pve_0.nii.gz') 

    WM_data = importImage(WM_filename)
    GM_data = importImage(GM_filename)
    CSF_data = importImage(CSF_filename)

    bvals = [float(b) for b in bvals_file.read().split('\n') if b.isdigit()]
    bvals_file.close()

    diffSpinEcho_data = []
    for b in bvals:
        tmp_data = diffSpinEcho(TE, TR, field, b, WM_data, GM_data, CSF_data)
        diffSpinEcho_data.append(addNoise(tmp_data,noise))

    diffSpinEcho_data = np.array(diffSpinEcho_data)
    plt.imshow(diffSpinEcho_data[0,:,:,90].T, cmap=plt.cm.gray)
    plt.xlim(reversed(plt.xlim()))
    plt.ylim(reversed(plt.ylim()))
    plt.axis('off')
    plt.show()

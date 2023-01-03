#Author: Gustavo Solcia
#E-mail: gustavo.solcia@usp.br

"""
Program created to import a gray and white matter segmentation from the BrainWeb dataset from 2006.
"""

import os
import numpy as np
import nibabel as nib

def importImage(filename):
    """Receives data path as string and returns numpy array data.

    """

    img = nib.load(filename)

    data_array = img.get_fdata()

    return data_array


if __name__=='__main__':
    import matplotlib.pyplot as plt

    path = '../data/dataset/'
    filename = os.path.join(path, 't1_brainweb_1mm.nii')

    data_array = importImage(filename)

    plt.imshow(data_array[:,:,90].T, cmap=plt.cm.gray)
    plt.xlim(reversed(plt.xlim()))
    plt.ylim(reversed(plt.ylim()))
    plt.axis('off')
    plt.show()

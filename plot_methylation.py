__author__ = 'ravi'
from matplotlib import pyplot as plt

import numpy as np
import sys

sample1,sample2 = sys.argv[1],sys.argv[2]


methyl_array1 = np.loadtxt(sample1,usecols=(3,4))
methyl_array2 = np.loadtxt(sample2,usecols=(3,4))

rep1 = methyl_array1[:,0]/(methyl_array1[:,1] + methyl_array1[:,0])

rep2 = methyl_array2[:,0]/(methyl_array2[:,1] + methyl_array2[:,0])


plt.scatter(rep1,rep2)

plt.show()
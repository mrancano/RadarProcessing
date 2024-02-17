#%%
import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
from PIL import Image
import radardataprocessor

tools = radardataprocessor.radartools("FinalProject")
filename = 'finalDataTranspose1'

number_of_samples1 = 1185
number_of_lines1 = 8357
header_length = 0
#create chirp
s = -1.037037*10**12 #Chirp Slope
tao = 2.7*10**-5 #interval
sampling_rate1 = 32*10**6


chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,number_of_samples1)
print("Chirp Created")
sample_array , sample_spectra = tools.createSampleArray(filename,number_of_lines1,number_of_samples1,header_length,number_of_samples1)
print("Sample Array Created")
image1 = Image.fromarray(8*np.abs(sample_array))

compressedarray = tools.compressSignal(chirp_vec,sample_array)

image2 = Image.fromarray(np.abs(compressedarray))

image1.show()
image2.show()
# %%
image2 = Image.fromarray(np.abs(compressedarray)/4)

image1.show()
image2.show()

# %%

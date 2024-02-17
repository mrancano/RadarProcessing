#%%
import numpy as np
import radardataprocessor
import matplotlib.pyplot as plt
import scipy as sp
from PIL import Image

#PROBLEM 1
#create chirp
s = 4.18966*10**11 #Chirp Slope
tao = 37.12*10**-6 #interval

sampling_rate1 = 18.96*10**6
number_of_samples1 = 4903
number_of_lines1 = 1024
header_length = 412


tools = radardataprocessor.radartools("Problem1")

chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,number_of_samples1)

sample_array , sample_spectra = tools.createSampleArray("ersdata",number_of_lines1,number_of_samples1,header_length,number_of_samples1)

#print(len(t_vec))
#print(len(chirp_vec))
#print(len(sample_array[1]))

raw_data = tools.rawdataArray("ersdata",number_of_lines1,number_of_samples1,header_length,number_of_samples1)

image1 = Image.fromarray(np.abs(sample_array))

compressedarray = tools.compressSignal(chirp_vec,sample_array)

image2 = Image.fromarray(np.abs(compressedarray))

#image1.show()
#image2.show()

image3 = Image.fromarray(raw_data)

image3.show()
# %%

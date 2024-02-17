#%%
import numpy as np
import radardataprocessor
import matplotlib.pyplot as plt
import scipy as sp
from PIL import Image
import copy 

#PRE-Processing
#create chirp
s = 4.189166*10**11 #Chirp Slope
tao = 37.12*10**-6 #interval

sampling_rate1 = 18.96*10**6
number_of_samples1 = 4903
number_of_lines1 = 10100
header_length = 412


tools = radardataprocessor.radartools("Problem1")

chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,number_of_samples1)
print("Chirp Created")
sample_array , sample_spectra = tools.createSampleArray("ersdata.hw4",number_of_lines1,number_of_samples1,header_length,number_of_samples1)
print("Sample Array Created")
image1 = Image.fromarray(8*np.abs(sample_array))

compressedarray = tools.compressSignal(chirp_vec,sample_array)

image2 = Image.fromarray(np.abs(compressedarray))

image1.show()
image2.show()

#%%
rawdata = tools.rawdataArray("ersdata.hw4",number_of_lines1,number_of_samples1,header_length)
image3 = Image.fromarray(rawdata)
image3.show()
# %%
#THIS SECTION WAS JUST FOR EXPERIMENTING-------
dopplerfreq = -300
compressedarray = compressedarray[:,:4200]
length = 10
velocity = 7125
wavelength = 0.0566
r_0 = 830000
r_dc = np.sqrt(r_0**2*(1+(wavelength*dopplerfreq/(2*velocity))**2))
az_samples = 10100
prf = 1679.9
chirp_rate = -velocity**2*2/(r_dc*wavelength)
sampling_rate_az = prf
tao_az = r_0*wavelength/(velocity*length)

az_chirp_vec , az_t_vec = tools.createchirp(chirp_rate,tao_az,prf,az_samples)

processed_array = tools.compressSignal(az_chirp_vec,compressedarray.T)

image3 = Image.fromarray(np.abs(processed_array.T)/64)

image3.show()
#-------------------------

# %%
#AZIMUTH COMPRESSSION
dr = velocity/(2*sampling_rate1)

processed_array = np.zeros(np.shape(compressedarray))

for i,bin in enumerate(compressedarray.T):
    if i % 50 == 0:
        print(i)
    '''
    if (i+1) % 800 == 0:
        break
    '''
    range1 = r_0+(i-1)*dr
    r_dc = np.sqrt(range1**2*(1+(wavelength*dopplerfreq/(2*velocity))**2))
    chirp_rate1 = -velocity**2*2/(r_dc*wavelength)
    tao_az1 = 0.8*r_dc*wavelength/(velocity*length)
    az_chirp_vec , az_t_vec = tools.createchirp(chirp_rate1,tao_az1,prf,az_samples)
    processed_bin = tools.compressSignal(az_chirp_vec,bin)
    processed_array[:,i] = processed_bin

image4 = Image.fromarray(np.abs(processed_array)/64)

image4.show()


# %%
#AVERAGING 5 PIXELS
nto1 = np.zeros((int(np.shape(processed_array)[0]/5)+1,int(np.shape(processed_array)[1])))
for i in range(np.shape(processed_array)[1]):
    az = processed_array[:,i]
    new = 0
    for j in range(len(az)):
        new+=np.abs(az[j])
        if j % 5 == 0:
            nto1[int(j/5),i] = new/5
            new = 0
#%%
image6 = Image.fromarray(nto1/64)
image6.show()
# %%
line = rawdata[:,412:]
print(np.sum(line)/np.size(line))
# %%

#%%
import numpy as np
import radardataprocessor
import matplotlib.pyplot as plt
import scipy as sp
from PIL import Image
import copy 
sampling_rate1 = 16*10**6
number_of_samples1 = 3444
number_of_lines1 = 20481
header_length = 112
tools = radardataprocessor.radartools("HW6")
sample_array , sample_spectra = tools.createSampleArray("data.hw6",number_of_lines1,number_of_samples1,header_length,number_of_samples1)
print("Sample Array Created")
# %%
image1 = Image.fromarray(8*np.abs(sample_array)) #CHECK HEADER LENGTH AND DATA ARRAY EXTRACTION
image1.show()
#%%
#RANGE PROCESSING
s = -5*10**11
tao = 1.1751132812499998e-05



chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,number_of_samples1)

compressedarray = tools.compressSignal(chirp_vec,sample_array)
print("Compressed Array Created")

image2 = Image.fromarray(np.abs(compressedarray))
image2.show()
# %%
slist = np.linspace(-4.8*10**11,-5.2*10**11,10)

taolist = np.linspace(10*10**-6,30*10**-6,10)
taolist = np.linspace(1.1751132812499998e-05,9.809277343749997e-06,5)


for i in range(len(taolist)):
    s = -5*10**11
    
    tao = 1.1751132812499998e-05
    chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,number_of_samples1)


    compressedarray = tools.compressSignal(chirp_vec,sample_array)

    image2 = Image.fromarray(np.abs(compressedarray))
    
    #plt.figure()
    #plt.title("tao = "+str(tao)+" s = "+str(round(s)))

    # Display the image within the matplotlib window
    #plt.imshow(image2)
    #plt.show()

    image2.show()
# %%
dopplerfreq = 100
length = 12
velocity = 7125
wavelength = 0.236057
r_0 = 696000/np.cos(35*np.pi/180)
r_dc = np.sqrt(r_0**2*(1+(wavelength*dopplerfreq/(2*velocity))**2))
az_samples = 20481
prf = 2159.827
chirp_rate = -velocity**2*2/(r_dc*wavelength)
sampling_rate_az = prf
tao_az = r_0*wavelength/(velocity*length)


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
#%%
#AVERAGING 4 PIXELS
nto1 = np.zeros((int(np.shape(processed_array)[0]/4)+1,int(np.shape(processed_array)[1])))
for i in range(np.shape(processed_array)[1]):
    az = processed_array[:,i]
    new = 0
    for j in range(len(az)):
        new+=np.abs(az[j])
        if j % 4 == 0:
            nto1[int(j/4),i] = new/4
            new = 0
image6 = Image.fromarray(nto1/64)
image6.show()
#%%
#FASTER AVERAGING 5 PIXELS
nto1 = np.zeros((int(np.shape(processed_array)[0]/5)+1,int(np.shape(processed_array)[1])))
for i in range(np.shape(nto1)[0]-1):
    nto1[i] = np.abs(processed_array[5*i])\
        +np.abs(processed_array[5*i+1])\
        +np.abs(processed_array[5*i+2])\
        +np.abs(processed_array[5*i+3])\
        +np.abs(processed_array[5*i+4])
    nto1[i] *= 1/5

image6 = Image.fromarray((nto1/64)[:,:-180])
image6.show()






#%%
#FOR DOPPLER CENTROID
prf = 1679.9 #Hz
azimuth_array = sample_array.transpose()

f_vector = sp.fft.fftfreq(number_of_lines1,1/prf)
mean_azimuth_spectrum = np.zeros(number_of_lines1,dtype=complex)
count=0
for azimuth_line in azimuth_array:  
    mean_azimuth_spectrum+=np.abs(sp.fft.fft(azimuth_array[count]))
    count+=1
mean_azimuth_spectrum /= len(azimuth_array)

plt.plot(f_vector,np.abs(mean_azimuth_spectrum),'.')
plt.plot(100,1500,'x')
plt.xlabel("f")
plt.ylabel("Magnitude")
plt.title("Average Azimuth Spectrum")

plt.show()
# %%

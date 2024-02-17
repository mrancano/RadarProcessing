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

#plt.plot(t_vec,np.abs(sample_array[1]))
#plt.show()``
#plt.plot(t_vec,np.abs(tools.compressSignal(chirp_vec,sample_array[1])))

#plt.show()


# %%
#PROBLEM 1
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
plt.xlabel("f")
plt.ylabel("Magnitude")
plt.title("Average Azimuth Spectrum")

plt.show()

#Problem 2

N_pulses = 64
doppler_freq = -300

steered_sample_array = copy.copy(compressedarray[(1-1)*N_pulses:1*N_pulses]) #COPY COMPRESSED
meanunshifted = np.zeros(N_pulses)
meanshifted = np.zeros(N_pulses)
#FINDING DOPPLER CENTROID
f_vector = sp.fft.fftfreq(N_pulses,1/prf)
for bin in range(number_of_samples1):
    az = steered_sample_array[:,bin]
    meanunshifted+=np.abs(sp.fft.fft(az))/number_of_samples1
    for k in enumerate(az):
        steered_sample_array[k[0],bin] = k[1]*np.exp(-1j*2*np.pi*doppler_freq*k[0]/prf)
    meanshifted+=np.abs(sp.fft.fft(steered_sample_array[:,bin]))/number_of_samples1


plt.plot(f_vector,meanunshifted,'.')
plt.plot(f_vector,meanshifted,'r.')





plt.xlabel("f")
plt.ylabel("Magnitude")
plt.title("Average Az Spectrum for 64 Pulse Patch")
plt.legend(["Unshifted","Shifted"])
plt.show()

#PATCHES
#GENERATE PATCHES
number_of_patches = np.floor((number_of_lines1/N_pulses))

list_of_patches = []

for i in range(int(number_of_patches)):
    if i%10 == 0: print(i)
    list_of_patches.append(tools.dopplerShiftPatch(N_pulses,doppler_freq,compressedarray,prf,i,number_of_samples1))
list_of_patches.pop(0)
#%%
#Single look
image5 = Image.fromarray(np.abs(sp.fft.fft(list_of_patches[1],axis=0))/16)
image5.show()

#DO FFT OF DOPPLER SHIFTED PATCH (AZIMUTH VECTOR)

#MULTILOOK 4 pixels in azimuth (add 4 and find average becomes new pixel)
print(len(list_of_patches))
print(list_of_patches[1])

    
# %%


pixel_spacing = 3.5

image_height = int(len(list_of_patches)*pixel_spacing+N_pulses)

multi_look_image = np.zeros((image_height,number_of_samples1))
multi_look_count = np.zeros(image_height)

list_of_patch_positions = []


for patch in enumerate(list_of_patches):
    
    start_new_patch = round(patch[0]*pixel_spacing)
    list_of_patch_positions.append(start_new_patch)
    
    multi_look_image[start_new_patch:(start_new_patch+N_pulses)]+=np.abs(np.fft.fft(patch[1],axis=0))
    multi_look_count[start_new_patch:(start_new_patch+N_pulses)]+=1

for i in range(len(multi_look_count)):
    if (multi_look_count[i]!=0):
        multi_look_image[i]/=multi_look_count[i]



image4 = Image.fromarray(multi_look_image/16)

image4.show()

print(np.shape(multi_look_image))

fourto1 = np.zeros((np.shape(multi_look_image)[0],int(np.shape(multi_look_image)[1]/4)+1))
for i in range(np.shape(multi_look_image)[0]):
    az = multi_look_image[i]
    new = 0
    for j in range(len(az)):
        new+=az[j]
        if j % 4 == 0:
            fourto1[i,int(j/4)] = new/4
            new = 0
image6 = Image.fromarray(fourto1/16)
image6.show()





# %%

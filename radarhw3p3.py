#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from PIL import Image 


#PROBLEM 1-------------------------
s1 = 4.189166*10**11 #Chirp Slope
tao = 37.12*10**-6 #interval

sampling_rate1 = 18.96*10**6
number_of_samples1 = int(9806*0.5)

def f1(t,s,tao=tao):
    if t>tao/2:
        re = 0
        im = 0
        return re,im
    re = np.cos(np.pi*s*t**2)
    im = np.sin(np.pi*s*t**2)

    return re,im

t_vector2 = np.arange(-tao/2,-tao/2+1/sampling_rate1*number_of_samples1,1/sampling_rate1)
f_vector2 = []

for el in t_vector2:
    f_vector2.append(f1(el,s1))

f_vector2 = np.array(f_vector2)

f_vector2_combined = f_vector2[:,0]+1j*f_vector2[:,1]

ft_vector2 = sp.fft.fftfreq(number_of_samples1,1/sampling_rate1)

ffreq_vector2 = sp.fft.fft(f_vector2_combined)

for el in enumerate(ffreq_vector2):
    ffreq_vector2[el[0]]=20*np.log10(np.linalg.norm(el[1]))




plt.figure("Freq")
plt.plot(ft_vector2,np.abs(ffreq_vector2),'.')
plt.ylabel("20Log(|S(f)|)")
plt.xlabel("f")
plt.title("Chirp")

plt.show()


number_of_lines = 1024
nheader_bytes = 412
nsample_bytes = 9806



f = open('ersdata','rb')





sample_array = np.zeros((number_of_lines,number_of_samples1),dtype=complex)
sample_spectra = np.zeros((number_of_lines,number_of_samples1),dtype=complex)

mean_spectra = np.zeros(number_of_samples1,dtype=complex)

for i in range(number_of_lines):
    header = f.read(nheader_bytes)
    sample_bytes = f.read(nsample_bytes)

    for el in enumerate(sample_bytes): #CHECK SUPPOSED SIZE
        if el[0]%2 != 0:
            sample_array[i,int(el[0]*0.5)] +=1j*(el[1]-15.5)
        else:
            sample_array[i,int(el[0]*0.5)] +=(el[1]-15.5)

    

    

    

    sample_spectra[i]=np.fft.fft(sample_array[i])

    mean_spectra+=np.abs(sample_spectra[i])

mean_spectra/=number_of_lines

plt.plot(t_vector2,np.abs(sample_array[1]))
plt.show()
        


compressedarray = np.zeros((number_of_lines,number_of_samples1),dtype=complex)
for i in range(number_of_lines):
    #compressedsample = np.correlate(sample_array[i],f_vector2_combined_reversed,'same')
    Rv = sp.fft.fft(f_vector2_combined)
    Sv = sp.fft.fft(sample_array[i])
    Svc = np.conj(Rv)*Sv
    
    compressedsample = sp.fft.ifft(Svc)
    
    compressedarray[i]= compressedsample


for el in enumerate(mean_spectra):
    mean_spectra[el[0]]=20*np.log10(np.linalg.norm(el[1]))


#plt.plot(t_vector2,f_vector2_combined,'.')
#plt.plot(t_vector2,compressedarray[100],'.')

plt.plot(ft_vector2,np.abs(ffreq_vector2),'r.')
plt.plot(ft_vector2,np.abs(mean_spectra),'b.')

plt.xlabel("Freq")
plt.ylabel("20Log(|S(f)|)")

plt.legend(["Reference","Mean of 1024 Spectra"])

#plt.ylim(-30,50)

plt.show()

f.close()    

image = Image.fromarray(np.abs(compressedarray))

image.show()


# %%

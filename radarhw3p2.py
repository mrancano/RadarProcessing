#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#PROBLEM 1
s1 = 10**12 #Chirp Slope
tao = 10**-5 #interval

sampling_rate1 = 10**8
number_of_samples1 = 2048

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


syntheticdatarecord = np.zeros(number_of_samples1*2,np.complex128)

syntheticdatarecord[100:(100+2048)] += f_vector2_combined
syntheticdatarecord[400:(400+2048)] += 5*f_vector2_combined
syntheticdatarecord[500:(500+2048)] += 2*f_vector2_combined

plt.plot(t_vector2,syntheticdatarecord[:number_of_samples1])
plt.plot(t_vector2,syntheticdatarecord[:number_of_samples1].imag)

plt.xlabel("t")
plt.ylabel("Signal")
plt.legend(["Real","Imaginary"])
plt.title("Synthetic Data Record")

plt.show()


Rv = sp.fft.fft(f_vector2_combined)
Sv = sp.fft.fft(syntheticdatarecord[:2048])
Svc = np.conj(Rv)*Sv
rc = sp.fft.ifft(Svc)



#plt.plot(t_vector2,compressedrecord)
plt.plot(t_vector2,np.abs(rc))
plt.xlabel("t")
plt.ylabel("Processed Signal")
plt.title("Chirp Discrimination")
#plt.plot(np.arange(-2047, 2048),compressedrecord)

plt.show()

# %%

#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#PROBLEM 1
s1 = 10**12 #Chirp Slope
tao = 10**-5 #interval

sampling_rate1 = 10**8
number_of_samples1 = 2048

omega_c = 200*np.pi
A_c = 1

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

'''
plt.figure("Problem1")
plt.plot(t_vector2,f_vector2[:,0],'r-')
plt.plot(t_vector2,f_vector2[:,1],'b-')
plt.legend(["Real","Imaginary"])
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Chirp")

plt.show()
'''

compressed_chirp1 = np.correlate(f_vector2_combined,f_vector2_combined,'same')

plt.figure("Chirp")
plt.plot(t_vector2,10*np.log10(np.abs(compressed_chirp1)))



plt.xlim(0.48e-5,0.57e-5)
plt.xlabel("t")
plt.ylabel("Response in dB for perfect reference")

plt.show()


#Problem 1b.
sref = 1.03*10**12 #Chirp Slope



f_vectorref = []

for el in t_vector2:
    f_vectorref.append(f1(el,sref))

f_vectorref = np.array(f_vectorref)



f_vectorref_combined = f_vectorref[:,0]+1j*f_vectorref[:,1]

compressed_chirp2 = np.correlate(f_vectorref_combined,f_vector2_combined,'same')

plt.figure("Chirp1")
plt.plot(t_vector2,10*np.log10(compressed_chirp2))

plt.xlabel("t")
plt.ylabel("Response in dB for s=1.03e12")

plt.xlim(0.48e-5,0.57e-5)

plt.show()
# %%
#

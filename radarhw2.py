#%%
import numpy as np
import matplotlib.pyplot as plt


#PROBLEM 1
sampling_rate1 = 5 #Hz

def f1(t):
    re = np.cos(np.pi*20*t)
    im = np.sin(np.pi*20*t)

    return re,im

t_vector1 = np.arange(0,1,1/sampling_rate1)

plt.figure("Real1")
plt.plot(t_vector1,f1(t_vector1)[0],'r-') #REAL
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Real")
plt.figure("Imaginary1")
plt.plot(t_vector1,f1(t_vector1)[1],'b-') #IMAG
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Imaginary")

plt.show()
# %%
import scipy as sp
#PROBLEM 2
s = 10**12 #Chirp Slope
tao = 10**-5 #interval

sampling_rate2 = 10**8
number_of_samples2 = 2048

def f2(t,s=s,tao=tao):
    if t>tao/2:
        re = 0
        im = 0
        return re,im
    re = np.cos(np.pi*s*t**2)
    im = np.sin(np.pi*s*t**2)

    return re,im

t_vector2 = np.arange(-tao/2,-tao/2+1/sampling_rate2*number_of_samples2,1/sampling_rate2)
f_vector2 = []

for el in t_vector2:
    f_vector2.append(f2(el))

f_vector2 = np.array(f_vector2)

plt.figure("Problem2")
plt.plot(t_vector2,f_vector2[:,0],'r-')
plt.plot(t_vector2,f_vector2[:,1],'b-')
plt.legend(["Real","Imaginary"])
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Chirp")

f_vector2_combined = f_vector2[:,0]+1j*f_vector2[:,1]

ft_vector2 = sp.fft.fftfreq(number_of_samples2,1/sampling_rate2)

ffreq_vector2 = sp.fft.fft(f_vector2_combined)


for el in enumerate(ffreq_vector2):
    ffreq_vector2[el[0]]=20*np.log10(np.linalg.norm(el[1]))




plt.figure("Freq")
plt.plot(ft_vector2,ffreq_vector2,'.')
plt.ylabel("20Log(|S(f)|)")
plt.xlabel("f")
plt.title("Chirp")
plt.show()





chirp_vec , t_vec = createchirp(s,tao,sampling_rate2,number_of_samples2)

print(chirp_vec)

plt.figure("Problem2")
plt.plot(t_vec,chirp_vec.real,'r-')
plt.plot(t_vec,chirp_vec.imag,'b-')
plt.legend(["Real","Imaginary"])
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Chirp")
plt.show()


# %%
#PROBLEM 3
s3 = 80
tao3 = 1

sampling_rate3 = 1000
number_of_samples3 = 2000
A_c = 1
omega_c = 200*np.pi


def f3(t,s=s3,tao=tao3):
    if t>tao/2:
        re = 0
        im = 0
        return re,im
    re = np.cos(np.pi*s*t**2)
    im = np.sin(np.pi*s*t**2)

    return re,im

t_vector3 = np.arange(-tao3/2,-tao3/2+1/sampling_rate3*number_of_samples3,1/sampling_rate3)
print(len(t_vector3))

f_vector3 = []

for el in t_vector3:
    f_vector3.append(f3(el))

f_vector3 = np.array(f_vector3)

f_vector3_combined = f_vector3[:,0]+f_vector3[:,1]*1j


re_vector = f_vector3[:,0]*A_c*np.cos(t_vector3*omega_c)
im_vector = f_vector3[:,1]*A_c*-np.sin(t_vector3*omega_c)

transmit_vector = re_vector+im_vector

plt.figure("IQ TRANSMIT")
plt.plot(t_vector3,transmit_vector,'-')
plt.xlabel("t")
plt.ylabel("f(t)")
plt.title("Transmitted Signal")
plt.show()

#%%
#PROBLEM 4
ft_vector3 = sp.fft.fftfreq(number_of_samples3,1/sampling_rate3)

freq_f_vector3 = sp.fft.fft(f_vector3_combined)


for el in enumerate(freq_f_vector3):
    freq_f_vector3[el[0]]=np.linalg.norm(el[1])

freq_transmit_vector = sp.fft.fft(transmit_vector)

for el in enumerate(freq_transmit_vector):
    freq_transmit_vector[el[0]]=np.linalg.norm(el[1])

plt.figure("Normal vs. Modulated")
plt.plot(ft_vector3,freq_f_vector3,'.')
plt.plot(ft_vector3,freq_transmit_vector,'.')
plt.legend(["Original","Transmitted"])
plt.xlabel("f")
plt.ylabel("|S(f)|")
plt.show()

#PROBLEM 4 part 2
received_vector = transmit_vector

I_prev = received_vector*2*np.cos(t_vector3*omega_c)
Q_prev = received_vector*-2*np.sin(t_vector3*omega_c)



fI_prev = sp.fft.fft(I_prev)
fQ_prev = sp.fft.fft(Q_prev)

#LOW PASS FILTER
LPFreq  = omega_c/(2*np.pi)
for el in enumerate(ft_vector3): 
    if np.abs(el[1])>LPFreq:
        fI_prev[el[0]] = 0
        fQ_prev[el[0]] = 0

I = sp.fft.ifft(fI_prev)
Q = sp.fft.ifft(fQ_prev)

demodulated = I+1j*Q



spectrum = sp.fft.fft(demodulated)

for el in enumerate(spectrum):
    spectrum[el[0]]=np.linalg.norm(el[1])

plt.plot(ft_vector3,freq_f_vector3)
plt.plot(ft_vector3,spectrum)

plt.xlabel("f")
plt.ylabel("|S(f)|")

plt.legend(["Transmitted","Received"])
plt.xlim((-100,100))



plt.show()


# %%

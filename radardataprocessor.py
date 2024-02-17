import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
from PIL import Image
#Class I made to modularise the code so I can focus on the question at hand -
class radartools():
    def __init__(self,name):
        self.name = name
    
    def createchirp(self,s,tao,f_sample,total_samples): #s:chirp slope tao:pulselength
        t_vec = np.arange(-tao / 2, -tao / 2 + 1 / f_sample * total_samples, 1 / f_sample)
    
        # Calculate the chirp signal directly using NumPy
        chirp_vec = np.where(t_vec <= tao / 2,np.exp(1j * np.pi * s * t_vec ** 2),0)
    
        return chirp_vec[:total_samples], t_vec[:total_samples]
    
    def compressSignal(self,reference,signal):
        Rv = sp.fft.fft(reference)
        Sv = sp.fft.fft(signal)
        Svc = np.conj(Rv)*Sv
        compressedsample = sp.fft.ifft(Svc)

        return compressedsample
    
    def createSampleArray(self,filename,number_of_lines,number_of_samples,header_length,sample_length):
        f = open(filename,'rb')

        sample_array = np.zeros((number_of_lines,number_of_samples),dtype=complex)
        sample_spectra = np.zeros((number_of_lines,number_of_samples),dtype=complex)

        for i in range(number_of_lines):
            header = f.read(header_length)
            sample_bytes = f.read(sample_length*2)

            for el in enumerate(sample_bytes): #CHECK SUPPOSED SIZE
                if el[0]%2 != 0:
                    sample_array[i,int(el[0]*0.5)] +=1j*(el[1]-15.5)
                else:
                    sample_array[i,int(el[0]*0.5)] +=(el[1]-15.5)

            sample_spectra[i]=np.fft.fft(sample_array[i])

        f.close()

        return sample_array,sample_spectra
    
    def rangeProcessArray(self,sample_array,number_of_lines,number_of_samples,reference):
        compressedarray = np.zeros((number_of_lines,number_of_samples),dtype=complex)
        for i in range(number_of_lines):
            compressedarray[i] = self.compressSignal(reference,sample_array[i])

        return compressedarray
    
    def rawdataArray(self,filename,number_of_lines,number_of_samples,header_length):
        f = open(filename,'rb')

        raw_data = np.zeros((number_of_lines,(number_of_samples*2+header_length)),dtype=np.uint8)

        for i in range(number_of_lines):
            line = f.read((number_of_samples*2+header_length))
            if i%(int(number_of_lines/10)) == 0:
                print(i)

            for el in enumerate(line):
                if el[0]>=header_length:
                    raw_data[i,el[0]] = 8*el[1]
                else:
                    raw_data[i,el[0]] = el[1]
        f.close()
        return raw_data

    def dopplerShiftPatch(self,N_pulses,doppler_freq,range_comp_array,prf,patch_number,number_of_samples):
        steered_sample_array = copy.copy(range_comp_array[(patch_number-1)*N_pulses:patch_number*N_pulses])
        #meanunshifted = np.zeros(N_pulses)
        #meanshifted = np.zeros(N_pulses)
    
        f_vector = sp.fft.fftfreq(N_pulses,1/prf)
        for bin in range(number_of_samples):
            az = steered_sample_array[:,bin]
            #meanunshifted+=np.abs(sp.fft.fft(az))/number_of_samples
            for k in enumerate(az):
                steered_sample_array[k[0],bin] = k[1]*np.exp(-1j*2*np.pi*doppler_freq*k[0]/prf) #CHANGE TO AVOID LOOP
            #meanshifted+=np.abs(sp.fft.fft(steered_sample_array[:,bin]))/number_of_samples



        

        #image3 = Image.fromarray(np.abs(steered_sample_array))
        #image3.show()

        return steered_sample_array

    

        
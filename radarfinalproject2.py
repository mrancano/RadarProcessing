#%%
import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
from PIL import Image
import radardataprocessor

tools = radardataprocessor.radartools("FinalProject")
filename = 'finalData'
file_length = 23249174
number_of_samples1 = 1391
number_of_lines1 = 8357
header_length = 0

f = open(filename,'rb')

raw_data = np.zeros((number_of_lines1,(number_of_samples1*2+header_length)),dtype=np.uint8)
#GENERATE RAW DATA IN CORRECT ORIENTATION
for i in range(number_of_samples1*2):
    line = f.read((number_of_lines1))
    

    for el in enumerate(line):
            raw_data[el[0],i] = el[1]

f.close()


image1 = Image.fromarray(raw_data)
image1.show()

raw_data.tofile("test1")
sample_array1,sample_spectra = tools.createSampleArray("test1",number_of_lines1 ,1185,412,1185)
image4 = Image.fromarray(8*np.abs(sample_array1))
image4.show()



#GET RID OF HEADER, COMBINE COMPLEX NUMBERS AND SHIFT

noheaderaw = raw_data[:,412:]
actualsamples = 1185
#number_of_samples1 = 1185 #Redefine number of samples without header
sample_array = np.zeros((number_of_lines1,actualsamples),dtype=complex)

for i in range(number_of_lines1):
    sample_bytes = noheaderaw[i]
    if i%(int(number_of_lines1/10)) == 0:
        print(i)
    for el in enumerate(sample_bytes): #CHECK SUPPOSED SIZE
        if el[0]%2 != 0:
            sample_array[i,int(el[0]*0.5)] +=1j*(el[1]-15.5)
        else:
            sample_array[i,int(el[0]*0.5)] +=(el[1]-15.5)




image2 = Image.fromarray(8*np.abs(sample_array))
image2.show()

#create chirp
s = -1.037037*10**12 #Chirp Slope
tao = 2.7*10**-5 #interval
sampling_rate1 = 32*10**6
actualsamples = 1185
chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,actualsamples)
print("Chirp Created")

compressedarray = tools.compressSignal(chirp_vec,sample_array1)

image2 = Image.fromarray(np.abs(compressedarray)/4)
image2.show()

#%%
import time
def getindeces(jaz0,ir0):
     '''
     r_near = 849780.1928671876
     dr = 4.6875
     dx = 3.340838432057467
     wavelength = 0.2360571
     l = 8.9
    '''


     dr = 6.25 #slant_range_spacing 
     dx = 1 #azimuth_pixel_spacing
     r_near = 4653
     wavelength = 0.25
     l = 2


     r0 = r_near + (ir0-1)*dr
     x0 = (jaz0-1)*dx

     w = wavelength*r0/(l)

     jazmin = round((x0-w/2)/dx+1)
     
     jazmax = round((x0+w/2)/dx+1)

     azindeces = np.arange(jazmin,jazmax+1,1)

     rindices = np.zeros((len(azindeces),1))

     rindices = (np.sqrt((r_near+(ir0-1)*dr)**2+(azindeces-jaz0)**2*dx**2)-r_near)/dr+1

     return azindeces-1,rindices-1

start = time.time()

azmax = 4978 #4978
azmin = 3377 ##3377
ranmax = 321
ranmin = 0

validaz = np.arange(azmin,azmax,1)
validrange = np.arange(ranmin,ranmax,1)
r_near = 849780.1928671876
dr = 4.6875
wavelength = 0.2360571

S = np.zeros_like(compressedarray)
for jaz0 in validaz:
     print(jaz0)
     for ir0 in validrange:
          azhist,rhist = getindeces(jaz0,ir0)
          for i in range(len(azhist)):
               jaz = azhist[i]
               ir = rhist[i]
               r = r_near+dr*ir
               
               S[jaz0,int(ir0)]+= compressedarray[jaz,int(ir)]*np.exp(1j*4*np.pi*wavelength*r)

end = time.time()

image3 = Image.fromarray(np.abs(S)/np.mean(np.abs(S))*64)
image3.show()



print(end-start)

# %%
image3 = Image.fromarray(np.abs(S[azmin:azmax,ranmin:ranmax])/np.mean(np.abs(S[azmin:azmax,ranmin:ranmax]))*64)
image3.show()

# %%

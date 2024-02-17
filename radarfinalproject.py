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

#%%
image1 = Image.fromarray(raw_data)
image1.show()
#%%
raw_data.tofile("test1")
sample_array1,sample_spectra = tools.createSampleArray("test1",number_of_lines1 ,1185,412,1185)
image4 = Image.fromarray(8*np.abs(sample_array1))
image4.show()


#%%
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



# %%
image2 = Image.fromarray(8*np.abs(sample_array))
image2.show()
# %%

#create chirp
s = -1.037037*10**12 #Chirp Slope
tao = 2.7*10**-5 #interval
sampling_rate1 = 32*10**6
actualsamples = 1185
chirp_vec , t_vec = tools.createchirp(s,tao,sampling_rate1,actualsamples)
print("Chirp Created")

compressedarray = tools.compressSignal(chirp_vec,sample_array1)
#%%
image2 = Image.fromarray(np.abs(compressedarray)/4)
image2.show()

# %%
import time
def getcurveindices(jaz0,ir0):
    dr = 4.6875 #slant_range_spacing 
    dx = 3.340838432057467 #azimuth_pixel_spacing
    r_near = 849780.1928671876
    wavelength = 0.2360571
    l = 8.9

    '''
    dr = 6.25 #slant_range_spacing 
    dx = 1 #azimuth_pixel_spacing
    r_near = 4653
    wavelength = 0.25
    l = 2
    '''
    r_0 = r_near + ir0*dr #index w/out -1 because Python zero-indexing
    x0 = jaz0*dx
    w = r_0*wavelength/l #azimuth beamwidth
    
    jazmin = round((x0-w/2)/dx)
    jazmax = round((x0+w/2)/dx)

    azindices = np.arange(jazmin,jazmax+1,1)
    rindices = (np.sqrt(r_0**2+(azindices-jaz0)**2*dx**2)-r_near)/dr

    return azindices,rindices

print(getcurveindices(99,199)[1][getcurveindices(99,199)[0]==449])

#%%
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


S = np.zeros((azmax-azmin,ranmax-ranmin),dtype=complex)
S1 = np.zeros((azmax-azmin,ranmax-ranmin),dtype=complex)

for jaz0 in validaz:
    if jaz0%10 ==0:
        print(jaz0)
    for ir0 in validrange:
        azindices,rindices = getcurveindices(jaz0,ir0)
        r =  r_near + dr*rindices
        #r1 = r_near + dr*np.floor(rindices)
        #r2 = r_near + dr*np.ceil(rindices)

        phasexp =  np.exp(1j*4*np.pi/wavelength*r)

        #phaseexp1 = np.exp(1j*4*np.pi/wavelength*r1)
        #phaseexp2 = np.exp(1j*4*np.pi/wavelength*r2)

        remainder = rindices-np.floor(rindices)
        remainderprime = 1-remainder
        
        lower = remainderprime*compressedarray[azindices,np.int8(np.floor(rindices))]
        upper = remainder*compressedarray[azindices,np.int8(np.ceil(rindices))]

        #p_lower = np.dot(lower,phaseexp1)
        #p_upper = np.dot(upper,phaseexp2)
        
        interp = lower+upper

        S[jaz0-azmin,ir0-ranmin] = np.dot(interp,phasexp)

        #S1[jaz0-azmin,ir0-ranmin] = p_lower+p_upper
    

end = time.time()

image3 = Image.fromarray(np.abs(S)/np.mean(np.abs(S))*64)
image3.show()



print(end-start)




# %%
image3 = Image.fromarray(np.abs(S)/np.mean(np.abs(S))*64)
image3.show()
# %%

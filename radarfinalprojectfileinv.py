#%%
import numpy as np
import scipy as sp
import copy
import matplotlib.pyplot as plt
from PIL import Image
import radardataprocessor

tools = radardataprocessor.radartools("FinalProject")
filename = 'finalData'
number_of_samples1 = 1185
number_of_lines1 = 8357
header_length = 412

raw_data = tools.rawdataArray(filename,number_of_lines1,number_of_samples1,header_length)

image1 = Image.fromarray(raw_data)
image1.show()
# %%
f1 = open(filename,'rb')
f2 = open("finalDataTranspose1",'wb')
f3 = open("finalDataTranspose2",'wb')

f1.read(header_length*number_of_lines1)


#assuming the file goes, azimuth line, azimuth line ....
for i in range(number_of_lines1):
    print(i)
    line = f1.read(number_of_samples1*2)
    f2.write(line)
#assuming the file goes, range line, range line
f1.seek(0)
f1.read(header_length*number_of_lines1)
temp_matrix = np.zeros((number_of_lines1,number_of_samples1*2),dtype=np.bytes_)
for i in np.arange(number_of_samples1*2,step=2):
    print(i)
    for j in range(number_of_lines1):
        number = f1.read(2)
        temp_matrix[j,i]=number[0]
        temp_matrix[j,i+1]=number[1]

for i in range(number_of_lines1):
    for j in range(number_of_samples1 * 2):
        f3.write(temp_matrix[i, j])




f1.close()
f2.close()
f3.close()

# %%

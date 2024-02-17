#%%
filename = "finalData"

with open(filename,'rb') as file:
    file.seek(0,2)
    file_length = file.tell()
    
print(file_length)


# %%
import numpy as np
import radardataprocessor
from PIL import Image
file_length = 143367000
number_of_lines = 20481
number_of_samples = 3500
header_length = 0


e = radardataprocessor.radartools('fileinv')

raw_data = e.rawdataArray(filename,number_of_lines,number_of_samples,header_length)

image1 = Image.fromarray(raw_data)

image1.show()

# %%
raw_data=raw_data.T

# %%
image2 = Image.fromarray(raw_data[:])
image2.show()
# %%
line = raw_data[:,112:]

print(np.size(line))
print(np.sum(line)/np.size(line))
# %%
line = raw_data[-4512:,112:144]

print(np.size(line))
print(np.sum(line)/np.size(line))
# %%

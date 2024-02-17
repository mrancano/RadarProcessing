#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from PIL import Image 

number_of_lines = 1024
nheader_bytes = 412
nsample_bytes = 9806

f = open('ersdata','rb')

raw_data = np.zeros((number_of_lines,(nsample_bytes+nheader_bytes)),dtype=np.uint8)

for i in range(number_of_lines):
    line = f.read((nsample_bytes+nheader_bytes))

    for el in enumerate(line):
        if el[0]>=412:
            raw_data[i,el[0]] = 8*el[1]
        else:
            raw_data[i,el[0]] = el[1]

image = Image.fromarray(raw_data)

image.show()
        
# %%

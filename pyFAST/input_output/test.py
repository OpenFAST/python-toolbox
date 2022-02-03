import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import pyFAST
from pyFAST.input_output.fast_output_file import FASTOutputFile, writeBinary


outb='C:/Work/DigiTwin-Stiesdal/code/F5T1RNA/Main_Spar_ED.outb'
outb='./Main_Spar_ED.outb'
outb2='./Main_Spar_ED2.outb'
outb3='./Main_Spar_ED3.outb'
outb22='./Main_Spar_ED22.outb'
outb33='./Main_Spar_ED33.outb'

# f = FASTOutputFile(outb)
# print(f.info)
# print(f.data.shape)

#     info = {'name': os.path.splitext(os.path.basename(filename))[0],
#             'description': DescStr,
#             'attribute_names': ChanName,
#             'attribute_units': ChanUnit}
# 


time = np.linspace(0,10,100)
C1   = np.cos(time)
C2   = time*0

colNames = ['Time','cos','cst']
colUnits = ['(s)','(deg)','(-)']

datalist = [time, C1, C2]

data     = np.column_stack([time, C1, C2])

print(np.min(data, axis=0))

# 
# nT, numOutWithTime = np.shape(f.data)
# 
# datalist = [[None]*nT for i in range(numOutWithTime)]
# 
# for j in range(numOutWithTime): 
#     datalist[j] = f.data[:,j]
# 
# 

# writeBinary2(outb2, datalist, f.info['attribute_names'], f.info['attribute_units'], fileID=2, descStr='')
# writeBinary2(outb2, datalist, colNames, colUnits, fileID=2, descStr='')
writeBinary(outb3, data, colNames, colUnits, fileID=2, descStr='')
# writeBinary(outb2, f.data, f.info['attribute_names'], f.info['attribute_units'], fileID=2, descStr='')

 
f3= FASTOutputFile(outb3)
f3.write(outb33)

# print(f2.toDataFrame())


if __name__ == '__main__':
    pass

#!/usr/bin/env python
import os
import pandas as pd


dfCD   = pd.read_csv('CD.csv',header=None)
dfCL   = pd.read_csv('CL.csv',header=None)
dfArea = pd.read_csv('Area.csv',header=None)

CDData   = dfCD.to_numpy()
CLData   = dfCL.to_numpy()
AreaData = dfArea.to_numpy()

num_rows, num_cols = CDData.shape

cdValues = CDData[:,20]
clValues = CLData[:,20]
areaValues = AreaData[:,20]


minCD = 10000
for i in range(num_rows):
    if(cdValues[i] < minCD):
        if(clValues[i] >= 0.32):
            if(areaValues[i] >= 0.081):
                bestIndex = i
                minCD = cdValues[i] 

print("Best index = ",bestIndex)
print("CD = ",cdValues[bestIndex])
print("CL = ",clValues[bestIndex])
print("Area = ",areaValues[bestIndex])


for i in range(num_rows):
    if(i != bestIndex):
        for k in range(20):
            CDData[i,k+21] = 0.0
            
df = pd.DataFrame(CDData)
filename = "CD.csv"
df.to_csv(filename, sep=',', index=False,header=None, encoding='utf-8')            

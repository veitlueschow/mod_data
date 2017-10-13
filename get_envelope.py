import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

u_x = [cline1[0,0],]
u_y = [cline1[0,1],]

for y in range(1,cline1.shape[0]-1):
    if (sign(cline1[y,1]-cline1[y-1,1])==1) and (sign(cline1[y,1]-cline1[y+1,1])==1):
        u_x.append(cline1[y,0])
        u_y.append(cline1[y,1])

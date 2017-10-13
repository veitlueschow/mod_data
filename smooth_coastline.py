from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d

def smooth(x,window_len=11,window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    #if x.ndim != 1:
        #raise ValueError, "smooth only accepts 1 dimension arrays."

    #if x.size < window_len:
        #raise ValueError, "Input vector needs to be bigger than window size."


    #if window_len<3:
        #return x


    #if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        #raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

filename='box_wmo.nc'
fh_wmo 		= Dataset(filename,mode='r')
filename='box_alat.nc'
fh_alat 	= Dataset(filename,mode='r')

wmo_ 	= fh_wmo.variables["wmo"]
lat_ 	= fh_wmo.variables["lat"]
lon_ 	= fh_wmo.variables["lon"]
depth_ 	= fh_wmo.variables["depth_15"]


var	= wmo_[0,56,:,:].copy()
lat	= lat_[:,:].copy()
lon	= lon_[:,:].copy()
depth = depth_[:].copy()

window_len = 101

#def smooth_coastline(var,lat,lon,window_len): 
coastline_	 = np.zeros((var.shape))
coastline	 = np.zeros((var.shape))
for y in range(var.shape[0]):
	for x in range(var.shape[1]):
		if var.mask[y,x] == True:
			coastline_[y,x] = 1
# delete stuff 1: southern part

for y in range(430,var.shape[0]):
	for x in range(200,var.shape[1]-10):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-6] == 1 and coastline_[y,x-5] == 1 and coastline_[y,x-4] == 1 and coastline_[y,x-3] == 1 and coastline_[y,x-2] == 1 and coastline_[y,x-1] == 1 and coastline_[y,x] == 1 and coastline_[y,x+1] == 0:
			coastline[y,x] = 1

for y in range(405,430):
	for x in range(380,80,-1):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-1] == 0 and coastline_[y,x] == 1:
			coastline[y,x] = 1
			
for y in range(280,405):
	for x in range(380,80,-1):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y-1,x] == 0 and coastline_[y,x] == 1:
			coastline[y,x] = 1
			
for y in range(150,280):
	for x in range(380,80,-1):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-1] == 0 and coastline_[y,x] == 1:
			coastline[y,x] = 1

for y in range(20,150):
	for x in range(80,480):
		if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y-2,x] == 1 and coastline_[y-1,x] == 1 and coastline_[y,x] == 1 and coastline_[y+1,x] == 0:
			coastline[y,x] = 1
for y in range(0,20):
	if np.count_nonzero(coastline[y,:]) == 0:
		tmp = 510-y*10
		coastline[y,tmp] = 1
	
for y in range(1,var.shape[0]):
	if np.count_nonzero(coastline[y,:]) == 0:
		coastline[y,np.flatnonzero(coastline[y-1,:])] = 1.
			
for y in range(1,var.shape[0]):
	tmp1 = np.flatnonzero(coastline[y,:])
	tmp2 = np.flatnonzero(coastline[y-1,:])	
	if abs(np.nonzero(coastline[y,:]) - np.flatnonzero(coastline[y-1,:])) > 70:
		coastline[y,tmp1] = 0
		coastline[y,tmp2] = 1

# Now devide into three parts

# First part from 100 - 400 (north america) ----------------------------

coastline1 = np.zeros((300,var.shape[1]))
j = 0
for y in range(100,400):
	coastline1[j,:] = coastline[y,:]
	j = j+1

cline1 = np.zeros((np.count_nonzero(coastline1),2))
j = 0
for y in range(100,400):
	#if np.count_nonzero(coastline1[j,:]) != 0:
	tmp = np.flatnonzero(coastline1[j,:])
	cline1[j,0] = lat[y,300]
	cline1[j,1] = lon[y,tmp]
	j = j +1

# try to find envelope

cine1_y = cline1[:,1].copy()
cine1_x = cline1[:,0].copy()

q_u = np.zeros(cine1_y.shape)
u_y = [cine1_y[0],]
u_x = [cine1_x[0],]

for y in range(1,len(cine1_y)-1):
    if (np.sign(cine1_y[y]-cine1_y[y-1])==1) and (np.sign(cine1_y[y]-cine1_y[y+1])==1):
        u_y.append(cine1_y[y])
        u_x.append(cine1_x[y])

u_x.append(cine1_x[len(cine1_y)-1])
u_y.append(cine1_y[len(cine1_y)-1])

u_p = interp1d(u_x,u_y, kind = 'cubic',bounds_error = False, fill_value=0.0)

for y in range(len(cine1_y)):
    q_u[y] = u_p(cine1_x[y])

#window_len = 71
tmp = smooth(q_u,window_len,"hanning") # only one possible smoothing configuration!
#tmp2 = interp1d(cline[:,1],cline[:,0])
cline_s1 = cline1.copy()
cline_s1[:,1]= tmp[window_len/2-1:cline1.shape[0]+window_len/2-1]

j=0
for y in range(100,400):
	tmp = np.flatnonzero(coastline1[j,:])
	cline_s1[j,0] = lat[y,300]
	j=j+1
	
# Second part from 400 - 700 (central america)

coastline2= np.zeros((300,var.shape[1]))
j = 0
for y in range(400,700):
	coastline2[j,:] = coastline[y,:]
	j = j+1

cline2 = np.zeros((np.count_nonzero(coastline2),2))
j = 0
for y in range(400,700):
	#if np.count_nonzero(coastline1[j,:]) != 0:
	tmp = np.flatnonzero(coastline2[j,:])
	cline2[j,0] = lat[y,300]
	cline2[j,1] = lon[y,tmp]
	j = j +1
	
cine2_y = cline2[:,1].copy()
#cine2_x = cline2[:,0].copy()

#window_len = 71
tmp = smooth(cine2_y,window_len,"hanning") # only one possible smoothing configuration!
#tmp2 = interp1d(cline[:,1],cline[:,0])
cline_s2 = cline2.copy()
cline_s2[:,1]= tmp[window_len/2-1:cline2.shape[0]+window_len/2-1]

j=0
for y in range(400,700):
	tmp = np.flatnonzero(coastline2[j,:])
	cline_s2[j,0] = lat[y,300]
	j=j+1
	
# Third part from 700 - 1100 (south america) ---------------------------

coastline3 = np.zeros((400,var.shape[1]))
j = 0
for y in range(700,1100):
	coastline3[j,:] = coastline[y,:]
	j = j+1

cline3 = np.zeros((np.count_nonzero(coastline3),2))
j = 0
for y in range(700,1100):
	#if np.count_nonzero(coastline1[j,:]) != 0:
	tmp = np.flatnonzero(coastline3[j,:])
	cline3[j,0] = lat[y,300]
	cline3[j,1] = lon[y,tmp]
	j = j +1

# try to find envelope

cine3_y = cline3[:,1].copy()

#window_len = 71
tmp = smooth(cine3_y,window_len,"hanning") # only one possible smoothing configuration!
#tmp2 = interp1d(cline[:,1],cline[:,0])
cline_s3 = cline3.copy()
cline_s3[:,1]= tmp[window_len/2-1:cline3.shape[0]+window_len/2-1]

j=0
for y in range(700,1100):
	tmp = np.flatnonzero(coastline3[j,:])
	cline_s3[j,0] = lat[y,300]
	j=j+1

# Put everything together. Check that no offset exists

tmp1 = cline_s1[cline_s1.shape[0]-1,1] - cline_s2[0,1]
if tmp1 <= 0:
	cline_s1[:,1] = cline_s1[:,1] - tmp1
elif tmp1 > 0:
	cline_s2[:,1] = cline_s2[:,1] + tmp1

tmp2 = cline_s2[cline_s2.shape[0]-1,1] - cline_s3[0,1]
#if tmp2 <= 0:
	#cline_s2[:,1] = cline_s2[:,1] - tmp2
#elif tmp2 > 0:
cline_s3[:,1] = cline_s3[:,1] + tmp2

cline = np.concatenate((cline1,cline2,cline3),axis=0)
cline_s = np.concatenate((cline_s1,cline_s2,cline_s3),axis=0)

#return [cline,cline_s]

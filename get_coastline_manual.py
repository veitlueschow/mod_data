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

def smooth_coastline(filename):
	#filename = "box3_1941.5"
	fh = Dataset(filename,mode='r')
	
	uko_ 	= fh.variables["uko"]
	vke_ 	= fh.variables["vke"]
	lat_2_ 	= fh.variables["lat_2"]
	lon_2_ 	= fh.variables["lon_2"]
	
	uko = uko_[0,0,:,:].copy()
	vke = vke_[0,0,:,:].copy()
	lat_2 = lat_2_[:,:].copy()
	lon_2 = lon_2_[:,:].copy()
	#a = uko.mask - vke.mask[:,0:570]
	
	coastline_	 = np.zeros((uko.shape))
	coastline	 = np.zeros((uko.shape))
	for y in range(uko.shape[0]):
		for x in range(uko.shape[1]):
			if uko.mask[y,x] == True:
				coastline_[y,x] = 1
	# delete stuff 1: southern part
	
	for y in range(430,uko.shape[0]):
		for x in range(200,uko.shape[1]-10):
			if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-6] == 1 and coastline_[y,x-5] == 1 and coastline_[y,x-4] == 1 and coastline_[y,x-3] == 1 and coastline_[y,x-2] == 1 and coastline_[y,x-1] == 1 and coastline_[y,x] == 1 and coastline_[y,x+1] == 0:
				coastline[y,x] = 1
	
	for y in range(150,430):
		for x in range(380,80,-1):
			#if np.count_nonzero(coastline[y,:]) == 0 and coastline_[y,x-3] == 0 and coastline_[y,x-2] == 0 and coastline_[y,x-1] == 0 and coastline_[y,x] == 1:
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
		
	
	for y in range(1,uko.shape[0]):
		if np.count_nonzero(coastline[y,:]) == 0:
			coastline[y,np.flatnonzero(coastline[y-1,:])] = 1.
				
	for y in range(1,uko.shape[0]):
		tmp1 = np.flatnonzero(coastline[y,:])
		tmp2 = np.flatnonzero(coastline[y-1,:])	
		if abs(np.nonzero(coastline[y,:]) - np.flatnonzero(coastline[y-1,:])) > 100:
			coastline[y,tmp1] = 0
			coastline[y,tmp2] = 1
			
	cline = np.zeros((np.count_nonzero(coastline),2))
	j = 0
	for y in range(coastline.shape[0]):
		if np.count_nonzero(coastline[y,:]) != 0:
			tmp = np.flatnonzero(coastline[y,:])
			cline[j,0] = lat_2[y,tmp]
			cline[j,1] = lon_2[y,tmp]
			j = j +1
	
	cline[315:390,1] = cline[315:390,1] + 2
	
	cline = cline[0:j,:].copy()
	window_len = 71
	tmp = smooth(cline[:,1],window_len,"hanning") # only one possible smoothing configuration!
	#tmp2 = interp1d(cline[:,1],cline[:,0])
	cline_s = cline.copy()
	cline_s[:,1]= tmp[window_len/2-1:cline.shape[0]+window_len/2-1]
	cline_tmp = cline.copy()
	for y in range(cline_s.shape[0]):
		if (cline_s[y,1]-cline[y,1]) > 0:
			cline_tmp[y,1] =cline_s[y,1]
	
	tmp = smooth(cline_tmp[:,1],window_len,"hanning")
	cline_s2 = cline.copy()
	cline_s2[:,1]= tmp[window_len/2-1:cline.shape[0]+window_len/2-1]
	
	for y in range(cline_s.shape[0]):
		if (cline_s2[y,1]-cline_tmp[y,1]) > 0:
			cline_tmp[y,1] =cline_s2[y,1]
	tmp = smooth(cline_tmp[:,1],window_len,"hanning")
	cline_s2 = cline.copy()
	cline_s2[:,1]= tmp[window_len/2-1:cline.shape[0]+window_len/2-1]

	for y in range(cline_s.shape[0]):
		if (cline_s2[y,1]-cline_tmp[y,1]) > 0:
			cline_tmp[y,1] =cline_s2[y,1]
	tmp = smooth(cline_tmp[:,1],window_len,"hanning")
	cline_s2 = cline.copy()
	cline_s2[:,1]= tmp[window_len/2-1:cline.shape[0]+window_len/2-1]
	
	window_len = 11
	tmp = smooth(cline[:,1],window_len,"hanning") # only one possible smoothing configuration!
	#tmp2 = interp1d(cline[:,1],cline[:,0])
	cline_s3 = cline.copy()
	#cline_s3[:,1]= tmp[window_len/2-1:cline.shape[0]+window_len/2-1]	
	
	#return [cline_s2, cline_tmp, coastline, coastline_]
	return [cline_s2, cline_s3,coastline]
	
def arc_length(cline_s):
	arc = np.zeros((cline_s.shape[0],1))
	for y in range(1,cline_s.shape[0]):
		arc[y] = arc[y-1]+np.sqrt((cline_s[y,1]-cline_s[y-1,1])**2+(cline_s[y,0]-cline_s[y-1,0])**2)
	return arc

def compute_normals(cline_s):
	
	tangent = np.ones((cline_s.shape[0],2))
	normal = np.ones((cline_s.shape[0],2))
	for y in range(1,cline_s.shape[0]-1):
		tangent[y,0] = cline_s[y+1,0]-cline_s[y,0]
		tangent[y,1] = cline_s[y+1,1]-cline_s[y,1]
		normal[y,0]  = -cline_s[y+1,1]+cline_s[y,1]
		normal[y,1]  = cline_s[y+1,0]-cline_s[y,0]
	
	tangent[0,0] = tangent[1,0]
	tangent[0,1] = tangent[1,1]
	tangent[cline_s.shape[0]-1,0] = tangent[cline_s.shape[0]-2,0]
	tangent[cline_s.shape[0]-1,1] = tangent[cline_s.shape[0]-2,1]
	
	normal[0,0] = normal[1,0]
	normal[0,1] = normal[1,1]
	normal[cline_s.shape[0]-1,0] = normal[cline_s.shape[0]-2,0]
	normal[cline_s.shape[0]-1,1] = normal[cline_s.shape[0]-2,1]
	
	# Normalize Vectors
	
	for y in range(cline_s.shape[0]):
		normal[y,0] = normal[y,0] / math.sqrt(normal[y,1]**2 + normal[y,0]**2)
		normal[y,1] = normal[y,1] / math.sqrt(normal[y,1]**2 + normal[y,0]**2)
		tangent[y,0] = tangent[y,0] / math.sqrt(tangent[y,1]**2 + tangent[y,0]**2)
		tangent[y,1] = tangent[y,1] / math.sqrt(tangent[y,1]**2 + tangent[y,0]**2)
	
	return [tangent, normal]

def check_neighbours(x,y,cline_s,dist,lon_2,lat_2):
	tmp1 = dist
	for j in range(y-100,y+100):
		if j>=0 and j <= (cline_s.shape[0]-1):
			tmp = np.sqrt((lon_2[j,x]-cline_s[j,1])**2+(lat_2[y,x]-lat_2[j,x])**2)
			#print y-j
			if tmp <= tmp1:
				tmp1 = tmp
				ymin = j
	return [tmp1, ymin]

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math

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

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y


#def smooth_coastline(filename):
filename = "box_-30_-15N_1941.5"
fh = Dataset(filename,mode='r')

uko_ 	= fh.variables["uko"]

uko = uko_[0,0,:,:].copy() # Choose timestep

# Find coastline. We got 200 longitudes values and 186 lats

coastline	 = np.zeros((uko.shape))

# First check horizontal direction
for y in range(uko.shape[0]):
	for x in range(1,uko.shape[1]-1):
		if uko.mask[y,x] == True and uko.mask[y,x+1] == False:
			coastline[y,x] = 1.
			break

cline = np.zeros((uko.shape[0]))

for y in range(uko.shape[0]):
	tmp = np.flatnonzero(coastline[y,:])
	cline[y] = tmp[0]
	
tmp = smooth(cline,25,"hanning") # only one possible smoothing configuration!
cline_s = cline.copy()
cline_s[:]= tmp[12:cline.shape[0]+12] # improve!
#return cline_s
	
def compute_normals(cline_s):
	
	tangent = np.ones((cline_s.shape[0],2))
	normal = np.ones((cline_s.shape[0],2))
	for x in range(1,cline_s.shape[0]-1):
		tangent[x,0] = cline_s[x+1]-cline_s[x-1]
		tangent[x,1] = 2.
		normal[x,0]  = 2.
		normal[x,1]  = -cline_s[x+1]+cline_s[x-1]
	
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

def check_neighbours(x,y,cline_s,dist):
	tmp1 = dist
	for j in range(y-100,y+100):
		if j>=0 and j <= (cline_s.shape[0]-1):
			tmp = np.sqrt((x-cline_s[j])**2+(y-j)**2)
			#print y-j
			if tmp <= tmp1:
				tmp1 = tmp
				ymin = j
	return [tmp1, ymin]
	
	
	
# Find nearest coast point

#def find_nearest(uko,cline_s,x)
	
	#for x in range(uko.shape[1]):
		#for y in range(uko.shape[0]):
			#if uko.mask[y,x] == False:
				#dist = abs(x -cline_s[y])
				#[minimaldist, ymin] = check_neighbours(x,y,cline_s,dist)
	#return ymin

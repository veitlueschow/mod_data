from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
import matplotlib as mpl
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable



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

def myround(x, prec=1, base=0.1):
  return round(base * round(float(x)/base),prec)
  

def netread_data(filename,varname):
	fh= Dataset(filename,mode='r')
	data_ 	= fh.variables[varname]
	data	= data_[0,:,:,:].copy()
	return data

def netread_grid(filename,latn,lonn,depthn):
	fh= Dataset(filename,mode='r')
	lat_ 	= fh.variables[latn]
	lat	= lat_[:,:].copy()
	lon_ 	= fh.variables[lonn]
	lon	= lon_[:,:].copy()
	depth_ 	= fh.variables[depthn]
	depth	= depth_[:].copy()
	return lat, lon, depth

def case_region(case):
	if case == "11N":
	    lx=200 # for average 11N - 18N
	    rx=350 # for average
	    llon=208 # for plot
	    rlon=240 # for plot
	if case == "10S":
	    lx=380 # for average 10S - 18S
	    rx=600 # for average
	    llon=471 # for plot
	    rlon=502 # for plot
	
	if case == "28S-40S":
	    lx= 20 # for average 10S - 18S
	    rx= 180 # for average
	    llon=97 # for plot
	    rlon=127 # for plot
	
	if case == "1S-18N":
	    lx= 10 # for average 10S - 18S
	    rx= 250 # for average
	    llon=40 # for plot
	    rlon=90 # for plot
	
	if case == "5S-10S":
	    lx= 5 # for average 5S - 10S
	    rx= 130 # for average
	    llon=67 # for plot
	    rlon=97 # for plsot
	
	if case == "1S":
	    lx= 10 # for average 10S - 18S
	    rx= 190 # for average
	    llon=30 # for plot
	    rlon=60 # for plot
	    
	if case == "1S-10S":
	    lx= 70 # for average 1S-10S
	    rx= 250 # for average
	    llon=105 # for plot
	    rlon=135 # for plot
	
	if case == "26N-29N":
	    lx= 10 # for average 1S-10S
	    rx= 130 # for average
	    llon=63 # for plot
	    rlon=90 # for plot
	    
	if case == "20N-24N":
	    lx= 3 # for average 20N-24N
	    rx= 110 # for average
	    llon=7 # for plot
	    rlon=36 # for plot

	if case == "10N-14N":
	    lx= 5 # for average 20N-24N
	    rx= 49 # for average
	    llon=7 # for plot
	    rlon=30 # for plot
	    
	return lx,rx,llon,rlon 
	
def case_region_meancirc(case):
	if case == "11N":
	    lx=200 # for average 11N - 18N
	    rx=350 # for average
	    llon=202 # for plot
	    rlon=250 # for plot
	if case == "10S":
	    lx=380 # for average 10S - 18S
	    rx=600 # for average
	    llon=471 # for plot
	    rlon=520 # for plot
	
	if case == "28S-40S":
	    lx= 20 # for average 10S - 18S
	    rx= 180 # for average
	    llon=97 # for plot
	    rlon=127 # for plot
	
	if case == "1S-18N":
	    lx= 10 # for average 10S - 18S
	    rx= 250 # for average
	    llon=40 # for plot
	    rlon=90 # for plot
	
	if case == "5S-10S":
	    lx= 5 # for average 5S - 10S
	    rx= 130 # for average
	    llon=67 # for plot
	    rlon=97 # for plsot
	
	if case == "1S":
	    lx= 10 # for average 10S - 18S
	    rx= 190 # for average
	    llon=30 # for plot
	    rlon=60 # for plot
	    
	if case == "1S-10S":
	    lx= 70 # for average 1S-10S
	    rx= 250 # for average
	    llon=105 # for plot
	    rlon=135 # for plot
	
	if case == "26N-29N":
	    lx= 10 # for average 1S-10S
	    rx= 130 # for average
	    llon=63 # for plot
	    rlon=90 # for plot
	    
	if case == "20N-24N":
	    lx= 3 # for average 20N-24N
	    rx= 110 # for average
	    llon=7 # for plot
	    rlon=36 # for plot

	if case == "10N-14N":
	    lx= 5 # for average 20N-24N
	    rx= 49 # for average
	    llon=7 # for plot
	    rlon=30 # for plot
	    
	return lx,rx,llon,rlon 

def case_path(case,time):
	if case=="10S":
	    path = "/work/mh0256/m300522/data_storm/eddies/"+time+"/10S-18S/aa"
	
	if case=="11N":
	    path = "/work/mh0256/m300522/data_storm/eddies/"+time+"/11N-18N/aa"
	
	if case=="28S-40S":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/28S-40S/aa"
	
	if case=="1S-18N":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/1S-18N/aa"
	    
	if case=="5S-10S":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/5S-10S/aa"
	    
	if case=="1S":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/1S-8N/aa"
	
	if case=="1S-10S":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/1S-10S/aa"
	
	if case=="26N-29N":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/26N-29N/aa"

	if case=="20N-24N":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/20N-24N/aa"

	if case=="10N-14N":
	    path = "/work/mh0256/m300522/data_storm/eddies/2000s/10N-14N/aa"
	return path

def depths_xz_sections(ztop,zbot,depth):
	tmp1 = min(depth[:], key=lambda x:abs(x-ztop))
	tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
	ktop = tmp[0]
	tmp1 = min(depth[:], key=lambda x:abs(x-zbot))
	tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
	kbot = tmp[0]
	tmp1 = min(depth[:], key=lambda x:abs(x-2000.))
	tmp = np.where(np.around(depth[:],decimals=1)==tmp1)[0]
	k2k = tmp[0]
	return ktop,kbot,k2k
	
def clines_from_uko_vke(kbot,ktop,k2k,uko,uko2,vke2):
	clines2 = np.zeros((kbot-ktop,uko.shape[1]))

	maxpos = np.zeros((vke2.shape[1]))
	tmp = vke2**2+uko2**2
	for j in range(vke2.shape[1]):
	        tmp1 = tmp[k2k,j,:].max()
	        tmp2 = np.where(tmp[k2k,j,:]==tmp1)[0]
	        maxpos[j] = tmp2[0]
	
	for k in range(0,kbot-ktop):
	    clines2[k,:] = maxpos[:]
	return clines2,maxpos

def plot_rho_upar(ax,x,z,data_upar,data_rho,title,pos):
	#ax.set_title(title)
	ax.set_xlabel("Distance from coast [km]",fontsize=15)	
	ax.patch.set_facecolor('grey')	
	w=(-0.2,-0.16,-0.1,-0.07,-0.04,-0.01,0.02,0.1)
	fvke = ax.contour(x, -z, data_upar, w,colors='k', linewidths=2)
	plt.clabel(fvke, fontsize=9, inline=1)
	w = (1027.35,1027.5,1027.65,1027.75,1027.8,1027.825,1027.85,1027.875)
	frho = ax.contour(x,-z,data_rho, w, colors='grey', linewidths=2)
	plt.clabel(frho, fontsize=10, inline=1,colors="grey",styles="oblique")
	#ax.set_yticklabels((-2500,-2000,-1500,-1000),fontsize=18)
	#if pos == "left":
		#plt.ylabel("Depth [m]", fontsize=15)
	#if pos == "center":
		#plt.yticks((),fontsize=18)
	#if pos == "right":
		#plt.yticks((),fontsize=18)
	return

def plot_rho_upar_meancirc(ax,x,z,data_upar,data_rho,title):
	#ax.set_title(title)
	ax.set_xlabel("Distance from coast [km]",fontsize=15)	
	ax.patch.set_facecolor('grey')	
	w=(-0.08,0.02)
	fvke = ax.contour(x, -z, data_upar, w,colors='k', linewidths=3)
	plt.clabel(fvke, fontsize=9, inline=1)
	w = (1027.35,1027.5,1027.65,1027.75,1027.8,1027.825,1027.85,1027.87,1027.88)
	frho = ax.contour(x,-z,data_rho, w, colors='grey', linewidths=2)
	plt.clabel(frho, fontsize=10, inline=1,colors="grey",styles="oblique")
	return

def plot_colorbar_all(fig,ax,fdat):
	cb_coord = [0.16,0.03,0.7,0.02]
	cbar_ax = fig.add_axes(cb_coord)
	cb = fig.colorbar(fdat,ax=ax,orientation="horizontal",format='%.0e',cax=cbar_ax)
	tick_locator = ticker.MaxNLocator(nbins=4)
	#cb =fig.colorbar(fdat,cax=ax,orientation='vertical')
	cb.locator = tick_locator
	cb.update_ticks()
	cb.ax.tick_params(labelsize=15)
	return

def plot_colorbar_single(fig,ax,fdat):
	divider = make_axes_locatable(ax)
	cax = divider.new_vertical(size="3%", pad=0.5, pack_start=True)
	fig.add_axes(cax)
	cb = fig.colorbar(fdat, cax=cax, orientation="horizontal",format='%.0e')
	tick_locator = ticker.MaxNLocator(nbins=5)
	#cb =fig.colorbar(fdat,cax=ax,orientation='vertical')
	cb.locator = tick_locator
	cb.update_ticks()
	return
	
def mask_according_to(ref,data):
	data_ = data.copy()
	output =  np.ma.masked_where(ref.mask==True, data_, copy=True)
	
	output_ = output.copy()
	output  = np.ma.masked_where(output == 0, output_,copy=True)
	
	return output

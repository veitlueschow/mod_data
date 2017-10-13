import get_coastline_part as coasti
from netCDF4 import Dataset
import numpy as np
import operator
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
class Args: pass

def prepare_data_alt(uko,vke,upar,uper,minimal,cline_s,minimaly):
	along 	= np.zeros((uko.shape[0],400)) # take only 100 across points
	across	= np.zeros((uko.shape[0],400))
	data1	= np.zeros((uko.shape[0],400))
	data2	= np.zeros((uko.shape[0],400))
	data3	= np.zeros((uko.shape[0],400))
	data4	= np.zeros((uko.shape[0],400))
	arc = coasti.arc_length(cline_s)
	
	for y in range(10,uko.shape[0]-10):
		i = 0
		for x in range(5,uko.shape[1]-5):
			if minimal[y,x] >= -1 and minimal[y,x] < 6 and i < 400 and uko.mask[y,x] == False:
				a = arc[minimaly[y,x]]
				along[y,i]=a[0]
				across[y,i]	= minimal[y,x]
				data1[y,i] 	= upar[y,x]
				data2[y,i] 	= uper[y,x]
				data3[y,i] 	= uko[y,x]
				data4[y,i] 	= vke[y,x]
				i = i+1

	tmp = data1[10,:]
	index, value = min(enumerate(tmp), key=operator.itemgetter(1))
	ref = across[10,index]

	
	acr = across.copy()
	for y in range(across.shape[0]):
		tmp = data1[y,:]
		index, value = min(enumerate(tmp), key=operator.itemgetter(1))
		diff = across[y,index] - ref
		acr[y,:] = across[y,:] -diff
				
	data1[data1 == 0] = np.nan
	data2[data2 == 0] = np.nan
	data3[data3 == 0] = np.nan
	data4[data4 == 0] = np.nan
	#along[along == 0] = np.nan
	#across[across == 0] = np.nan
	
	return [across,acr,along,data1,data2,data3,data4]

def prepare_data(uko,upar,minimal,cline_s,minimaly,rho):
	arc = coasti.arc_length(cline_s)
	
	# New coordiantes
	along 	= np.zeros((minimal.shape[0],1000))
	across	= np.zeros((minimal.shape[0],1000))	
	
	# Variables that carry the new coordiantes
	raw 	= Args()
	full	= Args()
	raw.upar	= np.zeros((minimal.shape[0],1000))
	raw.rho		= np.zeros((minimal.shape[0],1000))
	full.upar	= np.zeros((minimal.shape[0],1000))
	full.rho	= np.zeros((minimal.shape[0],1000))
		
	for y in range(len(arc)):
		along[y,:] = arc[y]
	
	for x in range(1000):
		across[:,x] = -1+x*0.01
	
	mini = np.around(minimal,2)
	
	for y in range(len(arc)):
		for x in range(minimal.shape[1]):
			if uko.mask[y,x] == False and mini[y,x] < 10 and mini[y,x] > -2 and mini[y,x] != 0:
				a = across[y,:]-mini[y,x]
				tmp = np.where(abs(a)<0.005)
				tmp = tmp[0]
				if len(tmp) != 0:
					tmp = tmp[0]
				#print tmp
					tmp2 = np.int16(minimaly[y,x])
					raw.upar[tmp2,tmp] = upar[y,x]
					raw.rho[tmp2,tmp]  = rho[y,x]
				#else:
					#print("no hit")
	data_raw_upar = raw.upar.copy()
	data_raw_rho = raw.rho.copy()
	for y in range(len(arc)):
		c1 = np.zeros((raw.upar.shape[1]))
		c1[:] = raw.upar[y,:]
		c2 = np.zeros((raw.rho.shape[1]))
		c2[:] = raw.rho[y,:]
		np.reshape(c1, len(c1))
		bad_indexes = (c1==0)
		c1[bad_indexes] = np.nan
		c2[bad_indexes] = np.nan
		good_indexes = np.logical_not(bad_indexes)
		good_data1 = c1[good_indexes].copy()
		good_data2 = c2[good_indexes].copy()
		if good_data1.size:
			interpolated1 = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data1)
			c1[bad_indexes] = interpolated1
			interpolated2 = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data2)
			c2[bad_indexes] = interpolated2
			full.upar[y,:] = c1.copy()
			full.rho[y,:] = c2.copy()
		
		tmp = np.nonzero(data_raw_upar[y,:])
		tmp = tmp[0]
		if tmp.size:
			tmp1 = min(tmp)
			full.upar[y,0:tmp1-1] = 0.
			tmp2 = max(tmp)
			full.upar[y,tmp2+1:full.upar.shape[1]] = 0.
		
		tmp = np.nonzero(data_raw_rho[y,:])
		tmp = tmp[0]
		if tmp.size:
			tmp1 = min(tmp)
			full.rho[y,0:tmp1-1] = np.nan
			tmp2 = max(tmp)
			full.rho[y,tmp2+1:full.rho.shape[1]] = np.nan
	
	tmp = data_raw_upar[10,:]
	index, value = min(enumerate(tmp), key=operator.itemgetter(1))
	ref = across[10,index]

	
	acr = across.copy()
	for y in range(across.shape[0]):
		tmp = data_raw_upar[y,:]
		index, value = min(enumerate(tmp), key=operator.itemgetter(1))
		diff = across[y,index] - ref
		acr[y,:] = across[y,:] -diff
	#acr=data_raw.rho - full.rho
	#acr=0
	return [along,across,acr,raw,full]
	



def along_across(lon,lat,data):
	data[data == 0] = np.nan
	fig = plt.figure(1)
	plt.set_cmap('bwr')
	v = np.linspace(-.115, 0.115, 150, endpoint=True)
	plt.contourf(lon,lat,data,v,vmin = -0.115, vmax = 0.115,extend='both')
	#plt.clim(-0.015,0.015)
	cbar = plt.colorbar(ticks=[-0.1,0,0.1])
	#cbar.set_clim(-0.1,0.1)
	return fig

def rho_along_across(lon,lat,data):
	data[data == 0] = np.nan
	fig = plt.figure(1)
	plt.set_cmap('Blues')
	v = np.linspace(1027.73, 1027.85, 150, endpoint=True)
	plt.contourf(lon,lat,data,v,vmin = 1027.73, vmax = 1027.85,extend='both')
	#plt.clim(-0.015,0.015)
	#cbar = plt.colorbar(ticks=[1027.84,1027.85],format="%.2f")
	cbar = plt.colorbar()
	#cbar = plt.colorbar(ticks=[1027.7,1027.75,1027.8,1027.85])
	#cbar.set_clim(-0.1,0.1)
	return fig

##def plots(lon_2,lat_2,upar,uper):

#figure3 = plt.figure(3)
##figure3.suptitle('Coast parallel', fontsize=20)
#decompo 	= plt.contour(lon_2,lat_2,upar,50)
#plt.clim(-0.15,0.15)
#CB		= plt.colorbar(decompo) 
#plt.show()
##figure3.savefig(filename+"_upar.png")

#figure4 = plt.figure(4)
#figure4.suptitle('Coast perpendicular', fontsize=20)
#decompo 	= plt.contour(lon_2,lat_2,uper,50)
#CB		= plt.colorbar(decompo) 
#plt.show()
##return[figure3,figure4]
##figure4.savefig(filename+"_uper.png")

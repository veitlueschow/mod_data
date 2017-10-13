import get_coastline_part as coasti
from netCDF4 import Dataset
import numpy as np
import operator
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
class Args: pass

def prepare_data_alt(inn,gulf,cline_s):
	along 	= np.zeros((inn.uko.shape[0],400)) # take only 100 across points
	across	= np.zeros((inn.uko.shape[0],400))
	data1	= np.zeros((inn.uko.shape[0],400))
	data2	= np.zeros((inn.uko.shape[0],400))
	data3	= np.zeros((inn.uko.shape[0],400))
	data4	= np.zeros((inn.uko.shape[0],400))
	arc = coasti.arc_length(cline_s)
	
	for y in range(10,inn.uko.shape[0]-10):
		i = 0
		for x in range(5,inn.uko.shape[1]-5):
			if gulf.minimal[y,x] >= -1 and gulf.minimal[y,x] < 6 and i < 400 and inn.uko.mask[y,x] == False:
				a = arc[gulf.minimaly[y,x]]
				along[y,i]=a[0]
				across[y,i]	= gulf.minimal[y,x]
				data1[y,i] 	= gulf.upar[y,x]
				data2[y,i] 	= gulf.uper[y,x]
				data3[y,i] 	= inn.uko[y,x]
				data4[y,i] 	= inn.vke[y,x]
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
	
	return [across,acr,along,data1,data2,data3,data4]

def prepare_data(inn,gulf,cline_s):
	arc = coasti.arc_length(cline_s)
	
	# New coordiantes
	along 	= np.zeros((gulf.minimal.shape[0],1000))
	across	= np.zeros((gulf.minimal.shape[0],1000))	
	
	# Variables that carry the new coordiantes
	raw 	= Args()
	full	= Args()
	raw.upar	= np.zeros((gulf.minimal.shape[0],1000))
	raw.rho		= np.zeros((gulf.minimal.shape[0],1000))
	raw.po		= np.zeros((gulf.minimal.shape[0],1000))
	full.upar	= np.zeros((gulf.minimal.shape[0],1000))
	full.rho	= np.zeros((gulf.minimal.shape[0],1000))
	full.po		= np.zeros((gulf.minimal.shape[0],1000))
		
	for y in range(len(arc)):
		along[y,:] = arc[y]
	
	for x in range(1000):
		across[:,x] = -1+x*0.01
	
	mini = np.around(gulf.minimal,2)
	#j=0
	for y in range(len(arc)):
		for x in range(gulf.minimal.shape[1]):
			if inn.uko.mask[y,x] == False and mini[y,x] < 10 and mini[y,x] > -2 and mini[y,x] != 0:
				a = across[y,:]-mini[y,x]
				tmp = np.where(abs(a)<0.005)
				tmp = tmp[0]
				if len(tmp) != 0:
					tmp = tmp[0]
				#print tmp
					tmp2 = np.int16(gulf.minimaly[y,x])
					if raw.upar[tmp2,tmp] == 0:
						raw.upar[tmp2,tmp] = gulf.upar[y,x]
						raw.rho[tmp2,tmp]  = inn.rho[y,x]
						raw.po[tmp2,tmp]  = inn.po[y,x]
				#else:
					#print("no hit")
					
	#print("doppelt", j)
	data_raw_upar = raw.upar.copy()
	data_raw_rho = raw.rho.copy()
	for y in range(len(arc)):
		c1 = np.zeros((raw.upar.shape[1]))
		c1[:] = raw.upar[y,:]
		c2 = np.zeros((raw.rho.shape[1]))
		c2[:] = raw.rho[y,:]
		c3 = np.zeros((raw.po.shape[1]))
		c3[:] = raw.po[y,:]
		np.reshape(c1, len(c1))
		bad_indexes = (c1==0)
		c1[bad_indexes] = np.nan
		c2[bad_indexes] = np.nan
		c3[bad_indexes] = np.nan
		good_indexes = np.logical_not(bad_indexes)
		good_data1 = c1[good_indexes].copy()
		good_data2 = c2[good_indexes].copy()
		good_data3 = c3[good_indexes].copy()
		if good_data1.size:
			interpolated1 = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data1)
			c1[bad_indexes] = interpolated1
			interpolated2 = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data2)
			c2[bad_indexes] = interpolated2
			interpolated3 = np.interp(bad_indexes.nonzero()[0], good_indexes.nonzero()[0], good_data3)
			c3[bad_indexes] = interpolated3
			full.upar[y,:] = c1.copy()
			full.rho[y,:] = c2.copy()
			full.po[y,:] = c3.copy()
		
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
			full.po[y,0:tmp1-1] = np.nan
			tmp2 = max(tmp)
			full.rho[y,tmp2+1:full.rho.shape[1]] = np.nan
			full.po[y,tmp2+1:full.po.shape[1]] = np.nan
	
	tmp = data_raw_upar[10,:]
	index, value = min(enumerate(tmp), key=operator.itemgetter(1))
	ref = across[10,index]

	
	acr = across.copy()
	ind = 0
	for y in range(across.shape[0]):
		tmp = data_raw_upar[y,:]
		index, value = min(enumerate(tmp), key=operator.itemgetter(1))
		diff = across[y,index] - ref
		acr[y,:] = across[y,:] -diff
		ind = ind + index
	ind = ind /(across.shape[0])
	raw.ind = np.int16(ind)

	return [along,across,acr,raw,full]
	



def upar_along_across(lon_,lat_,data_,i1,i2):
	lon = lon_[i1:i2,:]
	lat = lat_[i1:i2,:]
	data = data_[i1:i2,:].copy()
	data[data == 0] = np.nan
	fig = plt.figure(figsize=(18,6))
	plt.title("Parallel Flow")
	plt.set_cmap('bwr')
	v = np.linspace(-.115, 0.115, 150, endpoint=True)
	lon = lon*110.6
	lat = lat*110.6
	lat = lat-80
	lon = lon-min(lon[:,1])
	plt.contourf(lon,lat,data,v,vmin = -0.115, vmax = 0.115,extend='both')
	#plt.clim(-0.015,0.015)
	cbar = plt.colorbar(ticks=[-0.1,0,0.1])
	cbar.ax.tick_params(labelsize=17)
	#cbar.set_clim(-0.1,0.1)
	plt.ylim([-10,900])
	plt.xlabel('Along coast distance [km]',fontsize=20)
	plt.ylabel('Distance from coast [km]',fontsize=20)
	plt.xticks(fontsize = 17)
	plt.yticks(fontsize = 17)
	plt.savefig("upar.png",bbox_inches='tight')
	return fig

def rho_along_across(lon_,lat_,data_,ind,i1,i2):
	lon = lon_[i1:i2,:]
	lat = lat_[i1:i2,:]
	data = data_[i1:i2,:].copy()
	data[data == 0] = np.nan
	fig = plt.figure(figsize=(18,6))
	plt.set_cmap('Blues')
	plt.title("Potential Density",fontsize=20)
	# Find the range of the colorbar at vertain across-level
	tmp = data[:,ind] 
	mini = np.nanmin(tmp)
	maxi = np.nanmax(tmp)
	tmp = maxi-mini
	tmp = tmp*0.1
	mini = mini-tmp
	maxi = maxi+tmp
	v = np.linspace((mini),(maxi), 150, endpoint=True)
	lon = lon*110.6
	lat = lat*110.6
	lat = lat-80
	lon = lon-min(lon[:,1])
	plt.contourf(lon,lat,data,v,vmin = (mini), vmax = (maxi),extend='both')
	cbar = plt.colorbar(ticks=[1027.8,1027.81,1027.82,1027.83,1027.84],format="%.2f")
	cbar.ax.tick_params(labelsize=17)
	#cbar = plt.colorbar(format="%.2f")
	plt.ylim([-10,900])
	plt.xlabel('Along coast distance [km]',fontsize=20)
	plt.ylabel('Distance from coast [km]',fontsize=20)
	plt.xticks(fontsize = 17)
	plt.yticks(fontsize = 17)
	plt.savefig("rho.png",bbox_inches='tight')
	return fig

def po_along_across(lon_,lat_,data_,ind,i1,i2):
	lon = lon_[i1:i2,:]
	lat = lat_[i1:i2,:]
	data = data_[i1:i2,:].copy()
	data[data == 0] = np.nan
	fig = plt.figure(figsize=(18,6))
	plt.set_cmap('Blues')
	plt.title("Sea water pressure")
	# Find the range of the colorbar at certain across-level
	tmp = data[:,ind] 
	mini = np.nanmin(tmp)
	maxi = np.nanmax(tmp)
	tmp = maxi-mini
	tmp = tmp
	v = np.linspace((mini-tmp),(maxi+tmp), 250, endpoint=True)
	plt.contourf(lon,lat,data,v,vmin = (mini-tmp), vmax = (maxi+tmp),extend='both')
	#cbar = plt.colorbar(ticks=[1027.84,1027.85],format="%.2f")
	cbar = plt.colorbar(format="%.2f")
	plt.savefig("po.png",bbox_inches='tight')
	return fig

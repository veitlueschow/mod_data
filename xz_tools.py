from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import numpy.ma as ma
class Args: pass

def meri_average(vke,maxpos,var,ktop,kbot,lx,rx,llon,rlon):
	sum1	= np.zeros((vke.shape[0],vke.shape[2]))
	sum2 	= np.zeros((vke.shape[0],vke.shape[2]))
	ref_i 	= np.int64(maxpos[0])
	for i in range(lx,rx):
		for k in range(ktop,kbot):
			count = 0
			for j in range(vke.shape[1]):
				diff_i = np.int64(maxpos[j] - ref_i)
				if var.mask[k,j,i+diff_i] == False:
					count = count+1
					sum1[k,i] = sum1[k,i] + var[k,j,i+diff_i]

				sum2[k,i] = sum1[k,i] / count
	mean = ma.masked_array(data=sum2, mask = var.mask[:,0,:])

	data_plt = mean[ktop:kbot,llon:rlon]
	data_plt = ma.masked_invalid(data_plt,copy=False)

	return mean, data_plt

def along_average(upar,var,ktop,kbot,k2k,llon,rlon,maxpos,lat,lat_from,lat_to):
	tmp4 = np.where(np.around(lat[:,0],decimals=1)==lat_from)[0]
	tmp5 = tmp4[0]

	tmp4 = np.where(np.around(lat[:,0],decimals=1)==lat_to)[0]
	tmp6 = tmp4[0]

	sum1	= np.zeros((upar.shape[0],upar.shape[2]))
	sum2 	= np.zeros((upar.shape[0],upar.shape[2]))
	ref_i 	= np.int64(maxpos[0])
	for k in range(ktop,kbot):
		maxpos = np.zeros((upar.shape[1]))
		#tmp = vke2**2+uko2**2
		for j in range(upar.shape[1]):
			tmp1 = upar[k,j,:].min()
			tmp2 = np.where(upar[k,j,:]==tmp1)[0]
			maxpos[j] = tmp2[0]
		for i in range(upar.shape[2]):
			count = 0
			#for j in range(3,upar.shape[1]-2):
			for j in range(tmp5,tmp6):
			#for j in range(15,50):
				diff_i = np.int64(maxpos[j] - ref_i)
				if i+diff_i < upar.shape[2]:
					if var.mask[k,j,i+diff_i] == False:
						count = count + 1
						sum1[k,i] = sum1[k,i] + var[k,j,i+diff_i]
			if count > 5:
				sum2[k,i] = sum1[k,i] / count
	#mean = ma.masked_array(data=sum2, mask = var.mask[:,7,:])
	mean = sum2.copy()
	data_plt_ = mean[ktop:kbot,llon:rlon]
	data_plt = ma.masked_invalid(data_plt_,copy=False)
	return mean,data_plt

def along_average_dd(uko2,vke2,var,ktop,kbot,k2k,llon,rlon,maxpos,lat,lat_from,lat_to):
	tmp4 = np.where(np.around(lat[:,0],decimals=1)==lat_from)[0]
	tmp5 = tmp4[0]

	tmp4 = np.where(np.around(lat[:,0],decimals=1)==lat_to)[0]
	tmp6 = tmp4[0]

	sum1	= np.zeros((uko2.shape[0],uko2.shape[2]))
	sum2 	= np.zeros((uko2.shape[0],uko2.shape[2]))
	ref_i 	= np.int64(maxpos[0])
	for k in range(ktop,kbot):
		maxpos = np.zeros((uko2.shape[1]))
		tmp = vke2**2+uko2**2
		for j in range(tmp.shape[1]):
			tmp1 = tmp[k,j,:].max()
			tmp2 = np.where(tmp[k,j,:]==tmp1)[0]
			maxpos[j] = tmp2[0]
		for i in range(uko2.shape[2]):
			count = 0
			for j in range(tmp5,tmp6):
			#for j in range(15,50):
				diff_i = np.int64(maxpos[j] - ref_i)
				if i+diff_i < uko2.shape[2]:
					if var.mask[k,j,i+diff_i] == False:
						count = count + 1
						sum1[k,i] = sum1[k,i] + var[k,j,i+diff_i]
			if count > 5.:
				sum2[k,i] = sum1[k,i] / count
	#mean = ma.masked_array(data=sum2, mask = var.mask[:,7,:])
	mean = sum2.copy()
	data_plt_ = mean[ktop:kbot,llon:rlon]
	data_plt = ma.masked_invalid(data_plt_,copy=False)

	data_plt__ = np.ma.masked_where(data_plt==0, data_plt, copy=True)
	return mean,data_plt__


def meri_fit(vel,rho,lat,depth,shift):
	# mask regions that don't belong to the atlantic
	for j in range(0,450):
	    for i in range(0,120):
	        vel[:,j,i] = 0.

	for j in range(450,vel.shape[1]):
	    for i in range(0,250):
	        vel[:,j,i] = 0.
	vel1 = np.ma.masked_where(vel == 0, vel,copy=True)

	# Find maximum and write to 1D array
	#shift = 0
	#vel_ = vel1[k,:,:].copy()
	#rho_ = rho[k,:,:].copy()

	tmp = np.zeros((vel.shape[0],vel.shape[1]))
	for k in range(35,69):
	    for j in range(vel.shape[1]):
	        tmp2 = np.max(vel1[k,j,:])
	        tmp[k,j] = np.where(vel1[k,j,:]==tmp2)[0]
	    tmp3 = np.int64(tmp)

	rho_m = np.zeros((rho.shape[0],rho.shape[1]))
	for k in range(35,69):
	    for j in range(rho.shape[1]):
	        if tmp3[k,j]+shift < rho.shape[2]:
	            rho_m[k,j] = rho[k,j,tmp3[k,j]+shift]
	        else:
	            rho_m[k,j] = np.nan


	# Try now to vertically integrate from level 38 to 65. This is between 800m and 3200m, DWBC depth.
	# How to get layer depth?
	sum_rho = np.zeros(rho_m.shape[1])
	for z in range(38,68):
	    dz = depth[z+1]-depth[z]
	    sum_rho[:] = sum_rho[:] + dz*rho_m[z,:]

	# Now get correct lat coords, between -30 and 30 only!
	y = []
	rho_y = []
	for j in range(lat.shape[0]):
	    if np.abs(lat[j,tmp3[50,j]]) <= 25.:
	        y.append(lat[j,tmp3[50,j]])
	        rho_y.append(sum_rho[j])

	y = np.array(y)
	rho_y =np.array(rho_y)

	rho_fit = np.zeros(y.shape[0])
	idx = np.isfinite(rho_y)
	rho_fitparams = np.polyfit(y[idx],rho_y[idx],1)
	m = rho_fitparams[0]
	n = rho_fitparams[1]
	rho_fit[:] = n + m*np.float64(y[:])

	return rho_fit, rho_y, y, m

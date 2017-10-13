from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
import tools
import numpy.ma as ma

# This collection of funtions processes the vertical coordinate


# Get layer thickness

#filename='/work/mh0256/m300522/meta_storm/ddpo.nc'
def get_w_var_at_p_point(ifile,varname,ofile):
	filename1='/work/mh0256/m300522/meta_storm/ddpo.nc'
	fh_ddpo= Dataset(filename1,mode='r')
	
	ddpo_ = fh_ddpo.variables["ddpo"]
	depth_p_ = fh_ddpo.variables["depth_2"]
	lat_ 	= fh_ddpo.variables["lat"]
	lon_ 	= fh_ddpo.variables["lon"]
	
	ddpo = ddpo_[0,:,:,:].copy()
	depth_p = depth_p_[:].copy()
	lon	= lon_[:,:].copy()
	lat	= lat_[:,:].copy()
	
	# Get variable at w-point, it usually has one level more!
	
	#filename="/work/mh0256/m300522/data_storm/tape/1960s/wo_tm.nc"
	fh_wo = Dataset(ifile,mode="r")
	
	wo_ = fh_wo.variables[varname]
	depth_wo_ = fh_wo.variables["depth_4"]
	
	wo = wo_[0,:,:,:].copy()
	depth_wo = depth_wo_[:].copy()
	
	
	# Now compute new values for the variable at the w point at every depth_p point
	wo_p = ma.masked_array(data=np.zeros((ddpo.shape)), mask = ddpo.mask)
	
	for k in range(1,wo.shape[0]):
		ik = k-1
		wo_p[ik,:,:] = 0.5*(wo[ik]+wo[k])
	
	wo_p = ma.masked_array(data=wo_p, mask = ddpo.mask)
	
	write.write_3D(ofile,wo_p,varname,lat,lon,depth_p)
	
#def vertical_derrivative(filename,varname,ofile):

filename = "w+rho+.nc"
varname = "wrho_eddy"

filename1='/work/mh0256/m300522/meta_storm/ddpo.nc'
fh_ddpo= Dataset(filename1,mode='r')
ddpo_ = fh_ddpo.variables["ddpo"]
ddpo = ddpo_[0,:,:,:].copy()

fh_wrho = Dataset(filename,mode="r")
wrho_ = fh_wrho.variables[varname]
wrho = wrho_[0,:,:,:].copy()
dz_wrho = np.zeros((wrho.shape))
ddz = np.zeros((wrho.shape[0]))

#~ for j in range(wrho.shape[1]):
	#~ for i in range(wrho.shape[2]):
ddz[0] = ddpo[0,1488,2282] / 4.

for k in range(wrho.shape[0]-1):
	print(k)
	ik = k+1
	#if wrho.mask[k,j,i] == False and wrho.mask[ik,j,i] == False:
	dz = (ddpo[k,1488,2282] + ddpo[ik,1488,2282]) / 2
	dz_wrho[ik,:,:] = (wrho[k,:,:] - wrho[ik,:,:]) / dz
	ddz[ik] = ddz[k] + dz/2.
				
dz_wrho[0,:,:] = -wrho[0,:,:] * 2. / ddpo[0,1,1]

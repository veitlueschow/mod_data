from netCDF4 import Dataset
import numpy as np
class Args: pass

# Now write netCDF file for v
def write_depth_meridional(ofile,data,data_name,lat,depth):
	f1 = Dataset(ofile,'w',clobber=True, format='NETCDF4') #'w' stands for write
	ys 	= f1.createDimension('y', data.shape[1])
	zs 	= f1.createDimension('z', data.shape[0])
	meris 	= f1.createVariable('lat', np.float64,('z','y'))
	depths   = f1.createVariable('depth',np.float64,('z','y'))
	data1		= f1.createVariable(data_name,np.float64,('z','y'))
	
	meris[:,:]	= lat.copy()
	depths[:,]	= depth.copy()
	data1[:,:]	= data.copy()
	#lons[:] = np.arange(48, 51.1, 0.1)
	
	f1.Conventions = 'CF-1.0'
	
	f1.close()
	
def write_3D(ofile,data,data_name,lat,lon,depth):
	f1 = Dataset(ofile,'w',clobber=True, format='NETCDF4') #'w' stands for write
	f1.createDimension('x', data.shape[2])
	f1.createDimension('y', data.shape[1])
	f1.createDimension('depth', data.shape[0])
	lato 	= f1.createVariable('lat', np.float64,('y','x'))
	lono	= f1.createVariable('lon',np.float64,('y','x'))
	deptho   = f1.createVariable('depth',np.float64,('depth',))
	data1		= f1.createVariable(data_name,np.float64,('depth','y','x'), fill_value=9e33)
	data1.coordinates = "lat lon" 
	
	lato[:]	= lat.copy()
	lono[:]	= lon.copy()
	deptho[:]	= depth.copy()
	data1[:]	= data.copy()
	#lons[:] = np.arange(48, 51.1, 0.1)
	
	f1.Conventions = 'CF-1.0'
	
	f1.close()

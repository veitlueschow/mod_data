from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
import amoc_plots as aplot
class Args: pass
inn = Args()
pas = Args()

def myround(x, prec=2, base=1.):
  return round(base * round(float(x)/base),prec)

def full_slice(vke, ddue, dlxv, depth,rbek, ij):
# ij is index of the latitude at which the profile should be taken
	# Create coordinates	
	full_slice = np.zeros((80))
	
	for k in range(full_slice.shape[0]):
		print(k)
		for i in range(vke.shape[2]):
			if vke.mask[k,ij,i] == False:
				if rbek[ij,i] == 4. or rbek[ij,i] == 5. or rbek[ij,i] == 3. or rbek[ij,i] == 2. or rbek[ij,i] == 1.:
					full_slice[k] = full_slice[k] + vke[k,ij,i]*dlxv[ij,i]*ddue[k,ij,i] # Output in Sverdrup [m^3 s^-1]
	
	#write.write_depth_meridional("transport_meridional.nc",tr_mer,'transport_mer',y,z)
	
	return full_slice

def dwbc_slice(vke, ddue, dlxv, depth,rbek, ij):
# ij is index of the latitude at which the profile should be taken
	# Create coordinates	
	dwbc_slice = np.zeros((80))
	
	for k in range(dwbc_slice.shape[0]):
		print(k)
		#for i in range(680,725): # for ij = 1039
		for i in range(455,490): # for ij = 723
			if vke.mask[k,ij,i] == False:
				if rbek[ij,i] == 4. or rbek[ij,i] == 5. or rbek[ij,i] == 3. or rbek[ij,i] == 2. or rbek[ij,i] == 1.:
					dwbc_slice[k] = dwbc_slice[k] + vke[k,ij,i]*dlxv[ij,i]*ddue[k,ij,i] # Output in Sverdrup [m^3 s^-1]
	
	#write.write_depth_meridional("transport_meridional.nc",tr_mer,'transport_mer',y,z)
	
	return dwbc_slice
	

def plot(dwbc_slice,full_slice,depth,lat):
	x1 = dwbc_slice[3:71] # 70 is index of approx 4 km depth
	x2 = full_slice[3:71]
	y = -depth[3:71]
	
	figure = plt.figure()
	plt.plot(x1,y,label='DWBC',lw=2)
	plt.plot(x2,y,label='Whole Basin',lw=2)
	plt.axis([-2e6,2e6,-4000,0])
	
	plt.title('Meridional Transport at latitude ' + str(lat))
	plt.legend(loc='upper left',frameon=False)
	plt.xlim(-2e6,2e6)
	plt.xlabel('Meridional Transport [Sv]')
	plt.ylabel('Depth')
	plt.xticks(np.linspace(-2e6,2e6,5,endpoint=True),fontsize=14)
	plt.yticks(np.linspace(-4000,0,5,endpoint=True),fontsize=14)
	plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
	plt.axvline(x=0.,linewidth=1,linestyle="dashed",color="k")
	

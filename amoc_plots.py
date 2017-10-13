from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import get_coastline_amoc as coasti
from imp import reload
import write_netCDF as write
class Args: pass
read = Args()
pas = Args()

def myplot_diverging(y,z,data,lim):

	f3 = plt.figure()
	plt.set_cmap('coolwarm')
	v = np.linspace(-lim, lim, 100, endpoint=True)
	plt.contourf(y,-z,data[:,:],v,extend="both")
	plt.colorbar()
	
	return f3
	
def myplot_sequential(y,z,data):

	low 	= np.min(data)
	high 	= np.max(data)
	f3 = plt.figure()
	plt.set_cmap('rainbow')
	v = np.linspace(low, high, 20, endpoint=True)
	plt.contour(y,-z,data[:,:],v,extend="both")
	plt.colorbar()
	
	return f3

def myplot_multi(x,z,rho,x3,z3,vke,x2,z2,urho):
	f1 = plt.figure()
	
	fvke = plt.contour(x3, -z3, vke, 8,colors='k', linewidths=3)  # negative contours will be dashed by default
	plt.clabel(fvke, fontsize=9, inline=1)
	frho = plt.contour(x,-z,rho, 20, colors='white', linewidths=2)
	plt.clabel(frho, fontsize=9, inline=1)
	high 	= np.max(abs(urho))	
	high = high - 0.1*high
	v = np.linspace(-high, high, 50, endpoint=True)
	furho = plt.contourf(x2,-z2,urho,v,extend="both")
	furho.set_cmap('coolwarm')
	plt.colorbar(format='%.0e')
	plt.show()
	
	return f1, high
	
def myplot_multi2(x,z,rho,x3,z3,vke,x2,z2,urho,high,title):
	f1 = plt.figure()
	plt.title(title)
	fvke = plt.contour(x3, -z3, vke, 8,colors='k', linewidths=3)  # negative contours will be dashed by default
	plt.clabel(fvke, fontsize=9, inline=1)
	frho = plt.contour(x,-z,rho, 20, colors='white', linewidths=2)
	plt.clabel(frho, fontsize=9, inline=1)
	#high 	= np.max(abs(urho))	
	#high = high - 0.1*high
	v = np.linspace(-high, high, 50, endpoint=True)
	furho = plt.contourf(x2,-z2,urho,v,extend="both")
	furho.set_cmap('coolwarm')
	plt.colorbar(format='%.0e')
	plt.show()
	
	return f1
	
def myplot_quiver(x3,z3,vke,x,z,u,w,title):
	f1 = plt.figure()
	plt.title(title)
	fvke = plt.contourf(x3, -z3, vke, 40) 
	fvke.set_cmap('coolwarm')
	plt.colorbar(format='%.0e')
	fquiver = plt.quiver(x,-z,u,w,pivot='mid')
	plt.show()
	
	return f1

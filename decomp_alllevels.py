import sys
sys.path.insert(0, '/home/mpim/m300522/mod_data')
from netCDF4 import Dataset
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy.ma as ma
from imp import reload
import tools
	
def get_normals_dd(lat,lon,ktop,kbot,uko2,vke2):
	
	maxpos = np.zeros((kbot-ktop,uko2.shape[1]))
	tmp = vke2**2+uko2**2
	for k in range(0,18):
		ik = k+ktop
		for j in range(uko2.shape[1]):
			tmp1 = tmp[ik,j,:].max()
			tmp2 = np.where(tmp[ik,j,:]==tmp1)[0]
			maxpos[k,j] = 1*tmp2[0]
	
	for k in range(18,kbot-ktop):
		maxpos[k,:] = 1*maxpos[17,:]
	
	clines = 1*maxpos
	# Smooth coastline
	window_len=5
	clines2 = np.zeros((clines.shape[0],clines.shape[1]))
	for k in range(clines.shape[0]):
		tmp = tools.smooth(clines[k,:],window_len,"flat")
		clines2[k,:] = 1*tmp[np.int64(window_len/2-1):np.int64(clines.shape[1]+window_len/2-1)]

	#clines2 = 1*clines

	tangent = np.ones((kbot-ktop,clines.shape[1],2))
	normal = np.ones((kbot-ktop,clines.shape[1],2))
	for j in range(0,clines.shape[1]-1):
		# the tangent component points northward
		# the normal component points towards the coast!!!
		tangent[:,j,0] = 9900 # ty
		tangent[:,j,1] = 10855*(-clines2[:,j+1]+clines2[:,j]) # tx
		normal[:,j,0]  = -1*tangent[:,j,1] # ny
		normal[:,j,1]  = 1*tangent[:,j,0] # nx

	#tangent[:,0,0] = 1*tangent[:,1,0]
	#tangent[:,0,1] = 1*tangent[:,1,1]
	tangent[:,clines.shape[1]-1,0] = 1*tangent[:,clines.shape[1]-2,0]
	tangent[:,clines.shape[1]-1,1] = 1*tangent[:,clines.shape[1]-2,1]
	
	#normal[:,0,0] = 1*normal[:,1,0]
	#normal[:,0,1] = 1*normal[:,1,1]
	normal[:,clines.shape[1]-1,0] = 1*normal[:,clines.shape[1]-2,0]
	normal[:,clines.shape[1]-1,1] = 1*normal[:,clines.shape[1]-2,1]

	#tangent_ = tangent.copy()

	# Normalize Vectors
	normal_ = np.zeros((normal.shape))
	tangent_ = np.zeros((tangent.shape))
	
	for j in range(clines.shape[1]):
		normal_[:,j,0] = 1*normal[:,j,0] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		normal_[:,j,1] = 1*normal[:,j,1] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		tangent_[:,j,0] = 1*tangent[:,j,0] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		tangent_[:,j,1] = 1*tangent[:,j,1] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		
	return tangent_, normal_

def get_normals(clines,lat,lon,ktop,kbot):
	
	# Smooth coastline
	window_len=21
	clines2 = np.zeros((clines.shape[0],clines.shape[1]))
	for k in range(clines.shape[0]):
		tmp = tools.smooth(clines[k,:],window_len,"flat")
		clines2[k,:] = 1*tmp[np.int64(window_len/2-1):np.int64(clines.shape[1]+window_len/2-1)]


	tangent = np.ones((clines.shape[0],clines.shape[1],2))
	normal = np.ones((clines.shape[0],clines.shape[1],2))
	for j in range(0,clines.shape[1]-1):
		# the tangent component points northward
		# the normal component points towards the coast!!!
		tangent[:,j,0] = 9900 # ty
		tangent[:,j,1] = 10855*(-clines2[:,j+1]+clines2[:,j]) # tx
		normal[:,j,0]  = -1*tangent[:,j,1] # ny
		normal[:,j,1]  = 1*tangent[:,j,0] # nx

	#tangent[:,0,0] = 1*tangent[:,1,0]
	#tangent[:,0,1] = 1*tangent[:,1,1]
	tangent[:,clines.shape[1]-1,0] = 1*tangent[:,clines.shape[1]-2,0]
	tangent[:,clines.shape[1]-1,1] = 1*tangent[:,clines.shape[1]-2,1]
	
	#normal[:,0,0] = 1*normal[:,1,0]
	#normal[:,0,1] = 1*normal[:,1,1]
	normal[:,clines.shape[1]-1,0] = 1*normal[:,clines.shape[1]-2,0]
	normal[:,clines.shape[1]-1,1] = 1*normal[:,clines.shape[1]-2,1]

	#tangent_ = tangent.copy()

	# Normalize Vectors
	normal_ = np.zeros((normal.shape))
	tangent_ = np.zeros((tangent.shape))
	
	for j in range(clines.shape[1]):
		normal_[:,j,0] = 1*normal[:,j,0] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		normal_[:,j,1] = 1*normal[:,j,1] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		tangent_[:,j,0] = 1*tangent[:,j,0] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		tangent_[:,j,1] = 1*tangent[:,j,1] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		
	return tangent_, normal_



	
def decomposition(uko,vke,tangent,normal,lon,lat,clines,ktop,kbot):

	uper = np.zeros((uko.shape))
	upar = np.zeros((uko.shape))
	#global minimal
	#global minimaly
	minimal = np.zeros((uko.shape))
	minimalj = np.zeros((uko.shape))
	for k in range(ktop,kbot):
		ik = k - ktop
		for j in range(1,uko.shape[1]-1):
			for i in range(uko.shape[2]):
				if uko.mask[k,j,i] == False and vke.mask[k,j,i] == False: # positive distances
					#nx = normal[ik,j-1,1]/3 + normal[ik,j,1]/3 + normal[ik,j+1,1]/3
					#ny = normal[ik,j-1,0]/3 + normal[ik,j,0]/3 + normal[ik,j+1,0]/3
					#tx = tangent[ik,j-1,1]/3 + tangent[ik,j,1]/3 + tangent[ik,j+1,1]/3
					#ty = tangent[ik,j-1,0]/3 + tangent[ik,j,0]/3 + tangent[ik,j+1,0]/3
					#nx = normal[ik,j,1]
					#ny = normal[ik,j,0]
					#tx = tangent[ik,j,1]
					#ty = tangent[ik,j,0]				
					nx = normal[ik,j-1,1]/6 + normal[ik,j,1]*2/3 + normal[ik,j+1,1]/6
					ny = normal[ik,j-1,0]/6 + normal[ik,j,0]*2/3 + normal[ik,j+1,0]/6
					tx = tangent[ik,j-1,1]/6 + tangent[ik,j,1]*2/3 + tangent[ik,j+1,1]/6
					ty = tangent[ik,j-1,0]/6 + tangent[ik,j,0]*2/3 + tangent[ik,j+1,0]/6
					
					uper[k,j,i] = nx*uko[k,j,i] + ny*vke[k,j,i]
					upar[k,j,i] = ty*vke[k,j,i] + tx*uko[k,j,i]
	
		upar[:,0,:] = 1*upar[:,1,:]
		upar[:,upar.shape[1]-1,:] = 1*upar[:,upar.shape[1]-2,:]
		
		uper[:,0,:] = 1*uper[:,1,:]
		uper[:,uper.shape[1]-1,:] = 1*uper[:,uper.shape[1]-2,:]
	
	#upar = -upar
	uper1 = 1*uper  # uper now points away from the coast for positive values
	upar_ = ma.masked_where(upar==0,upar,copy=False)
	uper_ = ma.masked_where(uper1==0,uper1,copy=False)
	return upar_, uper_
	
def decomposition_test(uko,vke,tangent,normal,lon,lat,clines,ktop,kbot):
	
	#clines_lat = np.zeros((60,clines.shape[1]))
	#clines_lon = np.zeros((60,clines.shape[1]))

	#for k in range(60):
		#for j in range(clines.shape[1]):
			#clines_lat[k,j] = lat[j,np.int64(clines[k,j])]
			#clines_lon[k,j] = lon[j,np.int64(clines[k,j])]

	uper = np.zeros((uko.shape))
	upar = np.zeros((uko.shape))
	#global minimal
	#global minimaly
	minimal = np.zeros((uko.shape))
	minimalj = np.zeros((uko.shape))
	for k in range(ktop,kbot):
		ik = k - ktop
		for j in range(1,uko.shape[1]-1):
			for i in range(uko.shape[2]):
				if uko.mask[k,j,i] == False and vke.mask[k,j,i] == False: # positive distances
					nx = normal[ik,j-1,1]/3 + normal[ik,j,1]/3 + normal[ik,j+1,1]/3
					ny = normal[ik,j-1,0]/3 + normal[ik,j,0]/3 + normal[ik,j+1,0]/3
					tx = tangent[ik,j-1,1]/3 + tangent[ik,j,1]/3 + tangent[ik,j+1,1]/3
					ty = tangent[ik,j-1,0]/3 + tangent[ik,j,0]/3 + tangent[ik,j+1,0]/3
					
					
					uper[k,j,i] = nx*uko[k,j,i] + ny*vke[k,j,i]
					upar[k,j,i] = ty*vke[k,j,i] + tx*uko[k,j,i]
	
		upar[:,0,:] = 1*upar[:,1,:]
		upar[:,upar.shape[1]-1,:] = 1*upar[:,upar.shape[1]-2,:]
		
		uper[:,0,:] = 1*uper[:,1,:]
		uper[:,uper.shape[1]-1,:] = 1*uper[:,uper.shape[1]-2,:]
	
	#upar = -upar
	uper1 = 1*uper  # uper now points away from the coast for positive values
	upar_ = ma.masked_where(upar==0,upar,copy=False)
	uper_ = ma.masked_where(uper1==0,uper1,copy=False)
	return upar_, uper_

def decomposition_little_sophisticated(uko,vke,tangent,normal,lon,lat,clines,ktop,kbot):
	
	clines_lat = np.zeros((kbot-ktop,clines.shape[1]))
	clines_lon = np.zeros((kbot-ktop,clines.shape[1]))

	for k in range(0,kbot-ktop):
		for j in range(clines.shape[1]):
			clines_lat[k,j] = lat[j,np.int64(clines[k,j])]
			clines_lon[k,j] = lon[j,np.int64(clines[k,j])]

	uper = np.zeros((uko.shape))
	upar = np.zeros((uko.shape))
	#global minimal
	#global minimaly
	minimal = np.zeros((uko.shape))
	minimalj = np.zeros((uko.shape))
	for k in range(ktop,kbot):
		ik = k -ktop
		for j in range(uko.shape[1]):
			for i in range(uko.shape[2]):
				if uko.mask[k,j,i] == False and vke.mask[k,j,i] == False: # positive distances
					dist = abs(lon[j,i] -clines_lon[ik,j])
					[minimaldist, jmin] = check_neighbours(i,j,clines_lon[ik,:],dist,lon,lat)
					#minimal[k,j,i] = minimaldist
					minimalj[k,j,i] = jmin
					#print minimaldist
					if (jmin-1) >0 and (jmin+1) < uko.shape[1]:
						nx = normal[ik,jmin-1,1]/6 + normal[ik,jmin,1]*2/3 + normal[ik,jmin+1,1]/6
						ny = normal[ik,jmin-1,0]/6 + normal[ik,jmin,0]*2/3 + normal[ik,jmin+1,0]/6
						tx = tangent[ik,jmin-1,1]/6 + tangent[ik,jmin,1]*2/3 + tangent[ik,jmin+1,1]/6
						ty = tangent[ik,jmin-1,0]/6 + tangent[ik,jmin,0]*2/3 + tangent[ik,jmin+1,0]/6
					else:
						nx = normal[ik,jmin,1]
						ny = normal[ik,jmin,0]
						tx = tangent[ik,jmin,1]
						ty = tangent[ik,jmin,0]
					uper[k,j,i] = +nx*uko[k,j,i] + ny*vke[k,j,i]
					upar[k,j,i] = ty*vke[k,j,i] + tx*uko[k,j,i]


	
	#upar = -upar
	uper1 = 1*uper # uper now points away from the coast for positive values
	upar_ = ma.masked_where(upar==0,upar,copy=False)
	uper_ = ma.masked_where(uper1==0,uper1,copy=False)
	return upar_, uper_, minimalj

def decomposition_sophisticated(uko,vke,tangent,normal,lon,lat,clines):
	
	clines_lat = np.zeros((60,clines.shape[1]))
	clines_lon = np.zeros((60,clines.shape[1]))

	for k in range(60):
		for j in range(clines.shape[1]):
			clines_lat[k,j] = lat[j,np.int64(clines[k,j])]
			clines_lon[k,j] = lon[j,np.int64(clines[k,j])]

	uper = np.zeros((uko.shape))
	upar = np.zeros((uko.shape))
	#global minimal
	#global minimaly
	minimal = np.zeros((uko.shape))
	minimalj = np.zeros((uko.shape))
	for k in range(20,80):
		ik = k -20
		for j in range(uko.shape[1]):
			for i in range(uko.shape[2]):
				if uko.mask[k,j,i] == False and vke.mask[k,j,i] == False and lon[j,i] >= (clines_lon[ik,j]-1.): # positive distances
					dist = abs(lon[j,i] -clines_lon[ik,j])
					[minimaldist, jmin] = check_neighbours(i,j,clines_lon[ik,:],dist,lon,lat)
					minimal[k,j,i] = minimaldist
					minimalj[k,j,i] = jmin
					#print minimaldist
					if (jmin-1) >0 and (jmin+1) < uko.shape[1]:
						nx = normal[ik,jmin-1,1]/6 + normal[ik,jmin,1]*2/3 + normal[ik,jmin+1,1]/6
						ny = normal[ik,jmin-1,0]/6 + normal[ik,jmin,0]*2/3 + normal[ik,jmin+1,0]/6
						tx = tangent[ik,jmin-1,1]/6 + tangent[ik,jmin,1]*2/3 + tangent[ik,jmin+1,1]/6
						ty = tangent[ik,jmin-1,0]/6 + tangent[ik,jmin,0]*2/3 + tangent[ik,jmin+1,0]/6
					else:
						nx = normal[ik,jmin,1]
						ny = normal[ik,jmin,0]
						tx = tangent[ik,jmin,1]
						ty = tangent[ik,jmin,0]
					uper[k,j,i] = nx*uko[k,j,i] + ny*vke[k,j,i]
					upar[k,j,i] = ty*vke[k,j,i] + tx*uko[k,j,i]
					
					
				elif uko.mask[k,j,i] == False and vke.mask[k,j,i] == False and lon[j,i] < (clines_lon[ik,j]) and lon[j,i] >= (clines_lon[ik,j]-1): # negative distances
					dist = abs(lon[j,i] -clines_lon[ik,j])
					[minimaldist, jmin] = check_neighbours(i,j,clines_lon[ik,:],dist,lon,lat)
					minimal[k,j,i] = -minimaldist
					#print "negative", x,y, -minimaldist
					minimalj[k,j,i] = jmin
					#print minimaldist
					if (jmin-1) >0 and (jmin+1) < uko.shape[1]:
						nx = normal[ik,jmin-1,1]/4 + normal[ik,jmin,1]*2 + normal[ik,jmin+1,1]/4
						ny = normal[ik,jmin-1,0]/4 + normal[ik,jmin,0]/2 + normal[ik,jmin+1,0]/4
						tx = tangent[ik,jmin-1,1]/4 + tangent[ik,jmin,1]/2 + tangent[ik,jmin+1,1]/4
						ty = tangent[ik,jmin-1,0]/4 + tangent[ik,jmin,0]/2 + tangent[ik,jmin+1,0]/4
					else:
						nx = normal[ik,jmin,1]
						ny = normal[ik,jmin,0]
						tx = tangent[ik,jmin,1]
						ty = tangent[ik,jmin,0]
					uper[k,j,i] = nx*uko[k,j,i] + ny*vke[k,j,i]
					upar[k,j,i] = ty*vke[k,j,i] + tx*uko[k,j,i]


	
	#upar = -upar
	uper1 = 1*uper # uper now points away from the coast for positive values
	upar_ = ma.masked_where(upar==0,upar,copy=False)
	uper_ = ma.masked_where(uper1==0,uper1,copy=False)
	return upar_, uper_, minimalj

def get_normals_sophisticated(clines,lat,lon):
	clines_lat_ = np.zeros((60,clines.shape[1]))
	clines_lon_ = np.zeros((60,clines.shape[1]))
	clines_lat = np.zeros((60,clines.shape[1]))
	clines_lon = np.zeros((60,clines.shape[1]))

	for k in range(clines.shape[0]):
		for j in range(clines.shape[1]):
			clines_lat_[k,j] = lat[j,np.int64(clines[k,j])]
			clines_lon_[k,j] = lon[j,np.int64(clines[k,j])]
	
	# Smooth coastline
	window_len=31
	for k in range(clines.shape[0]):
		tmp = tools.smooth(clines_lat_[k,:],window_len,"flat")
		clines_lat[k,:] = tmp[np.int64(window_len/2-1):np.int64(clines_lat_.shape[1]+window_len/2-1)]
		
		tmp = tools.smooth(clines_lon_[k,:],window_len,"flat")
		clines_lon[k,:] = tmp[np.int64(window_len/2-1):np.int64(clines_lon_.shape[1]+window_len/2-1)]

	tangent = np.ones((60,clines.shape[1],2))
	normal = np.ones((60,clines.shape[1],2))
	for j in range(1,clines.shape[1]-1):
		tangent[:,j,0] = -clines_lat[:,j+1]+clines_lat[:,j] # ty
		tangent[:,j,1] = -clines_lon[:,j+1]+clines_lon[:,j] # tx
		normal[:,j,0]  = -tangent[:,j,1] # ny
		normal[:,j,1]  = tangent[:,j,0] # nx

	tangent[:,0,0] = tangent[:,1,0]
	tangent[:,0,1] = tangent[:,1,1]
	tangent[:,clines.shape[1]-1,0] = tangent[:,clines.shape[1]-2,0]
	tangent[:,clines.shape[1]-1,1] = tangent[:,clines.shape[1]-2,1]
	
	normal[:,0,0] = normal[:,1,0]
	normal[:,0,1] = normal[:,1,1]
	normal[:,clines.shape[1]-1,0] = normal[:,clines.shape[1]-2,0]
	normal[:,clines.shape[1]-1,1] = normal[:,clines.shape[1]-2,1]

	# Normalize Vectors
	normal_ = np.zeros((normal.shape))
	tangent_ = np.zeros((tangent.shape))
	
	for j in range(clines.shape[1]):
		normal_[:,j,0] = 1*normal[:,j,0] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		normal_[:,j,1] = 1*normal[:,j,1] / np.sqrt(normal[:,j,1]*normal[:,j,1] + normal[:,j,0]*normal[:,j,0])
		tangent_[:,j,0] = 1*tangent[:,j,0] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		tangent_[:,j,1] = 1*tangent[:,j,1] / np.sqrt(tangent[:,j,1]*tangent[:,j,1] + tangent[:,j,0]*tangent[:,j,0])
		
	return tangent_, normal_

def check_neighbours(x,y,clines_lon,dist,lon,lat):
	tmp1 = dist
	for j in range(y-10,y+10):
		if j>=0 and j <= (clines_lon.shape[0]-1):
			tmp = np.sqrt((lon[j,x]-clines_lon[j])**2+(lat[y,x]-lat[j,x])**2)
			#print y-j
			if tmp <= tmp1:
				tmp1 = tmp
				ymin = j
	return [tmp1, ymin]

def read_normal_get_tangent(nx,ny,clines2,k2k):
	normal = np.zeros((60,clines2.shape[1],2))
	tangent = np.zeros((60,clines2.shape[1],2))
	
	for j in range(clines2.shape[1]):
		normal[:,j,0] = ny[k2k,j,clines2[0,j]+5]
		normal[:,j,1] = -nx[k2k,j,clines2[0,j]+5]

	tangent[:,:,1] = -normal[:,:,0]
	tangent[:,:,0] = normal[:,:,1]
	
	return normal, tangent

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math

# define search function
def find(x,y,coastline):
		above = np.flatnonzero(coastline[y-1,:])
		below = np.flatnonzero(coastline[y+1,:])
		if np.count_nonzero(coastline[y+1,:]) == 0:
			below = np.array([0])
		return [above, below]
		
def get_normal(x,y,max_a,max_b):
	normal = -(max_a - max_b) / 2.
	return normal;
	


def normals(filename):
	
	fh = Dataset(filename,mode='r')
	
	uko_ 	= fh.variables["uko"]
	lat_2_ 	= fh.variables["lat_2"]
	lon_2_ 	= fh.variables["lon_2"]
	
	uko = uko_[0,0,:,:].copy()
	
	lat_2 = lat_2_[:,:].copy()
	lon_2 = lon_2_[:,:].copy()
	
	# Find coastline. We got 200 longitudes values and 186 lats
	
	coastline	 = np.zeros((uko.shape))
	
	# First check horizontal direction
	for y in range(uko.shape[0]):
		for x in range(1,uko.shape[1]-1):
			if uko.mask[y,x] == True and uko.mask[y,x+1] == False:
				coastline[y,x] = 1.
				break
	
	# Add stuff in vertical direction
	for x in range(uko.shape[1]):
		if np.count_nonzero(coastline[:,x]) == 0:
			for y in range(1,uko.shape[0]-1):
				if uko.mask[y,x] == True and uko.mask[y-1,x] == True and uko.mask[y-1,x] == False and uko.mask[y+1,x]:
					coastline[y,x] = 1.
					break
					
	for y in range(1,uko.shape[0]-1):
		tmp = np.flatnonzero(coastline[y,:])
		if np.count_nonzero(coastline[y,:]) > 1 and (tmp[-1]-tmp[-1-1]) != 1:
			coastline[y,tmp[-1]] = 0
#			print "Island deleted y", y
	
	for x in range(1,uko.shape[1]-1):
		tmp = np.flatnonzero(coastline[:,x])
		if np.count_nonzero(coastline[:,x]) > 1 and (tmp[-1]-tmp[-1-1]) != 1:
			coastline[y,tmp[-1]] = 0
			#print "Island deleted x", x
			
	# Find orientation of coast line for each point
	# The normal vector has a x component in position "0" and the y comp at "1"
	normal = np.ones((uko.shape[0],2))
	normal[0,:] = 1.
			
	for y in range(1,uko.shape[0]-1):
		tmp = np.flatnonzero(coastline[y,:])
		if np.count_nonzero(coastline[y,:]) == 1:
			for x in range(1,uko.shape[1]-1):
				if coastline[y,x] == 1.:
					[above, below] = find(x,y,coastline)
					normal[y,1] = get_normal(x,y,max(above),max(below))
	#				print y,above,x,below
	#				print lon[y,x],lat[y,x], normal
		elif np.count_nonzero(coastline[y,:]) == 2:
			[above, below] = find(x,y,coastline)
			tmp = np.flatnonzero(coastline[y,:])
	#		print "two points", y
			if max(below) <= min(tmp): # works only in this region!
				normal[y,1] = -get_normal(x,y,max(below),max(tmp))
		elif np.count_nonzero(coastline[y,:]) > 2: # horizontal line
			[above, below] = find(x,y,coastline)
			tmp = np.flatnonzero(coastline[y,:])
	#		print "more points", y
			normal[y,0] = 0.
			normal[y,1] = -1.
	#	print lat[y,1], normal[y], y
	normal[180:186,1] = -0.5 # End of the window
	normal[0,1] = 0
	
	# normalize the vector 
	for y in range(coastline.shape[0]):
		normal[y,0] = normal[y,0] / math.sqrt(normal[y,1]**2 + normal[y,0]**2)
		normal[y,1] = normal[y,1] / math.sqrt(normal[y,1]**2 + normal[y,0]**2)
	
	
	# Now construct tangential vector which is orthogonal to normal
	
	tangent = normal.copy()
	tangent[:,0] = -1*normal[:,1]
	tangent[:,1] = normal[:,0]
	
	for y in range(coastline.shape[0]):
		tangent[y,0] = tangent[y,0] / math.sqrt(tangent[y,1]**2 + tangent[y,0]**2)
		tangent[y,1] = tangent[y,1] / math.sqrt(tangent[y,1]**2 + tangent[y,0]**2)	
	
	return [normal, tangent, coastline];

#!/usr/bin/env python

from codes.mesh import makemesh
### from codes.review3d import *
import scipy
import scipy.spatial

machine_eps = eps = np.finfo(float).eps

"""
pseudocode:

read in the structure
get the average surface
get the mesh for the average surface with the average box
get the maximum fluctuation
set the pivot to the middle height in the box given by the maximum fluctuation plus ten percent
take a frame with its box
pretend it is in the average box and literalize
for each point on the average surface, find the corresponding nearest point on the literalized frame
given this mapping, take the projection of that point onto the tangent plane, still in the average box vector
	using the average box is a necessary approximation otherwise we would have to redo the mesh each time
		although that's another way to go

To measure fluctuations we loop over all mesh points and find the point nearest them on the 
instantaneous surface. Then we project the instantaneous point into that mesh point.
The key thing is that everything is done from the perspective of the mesh points on the average surface.
"""

###---TOOLS

def literalize(heights,vec):
	"""Hasty function to turn heights into points."""
	nx,ny = heights.shape
	return np.concatenate((np.concatenate(np.transpose(np.meshgrid(*[np.linspace(0,vec[ii],n) 
		for ii,n in enumerate([nx,ny])]))).T,[heights.reshape(-1).T])).T

def height_recenter(pts,pivot,maxflux): 
	"""Recenter all heights to the mean of the average height, shifted up one half box vector."""
	pts[...,2] += maxflux - pivot
	return pts

def boxstuff(pts,vec): 
	"""Ensure that the points are inside the box."""
	return pts-(pts>vec)*vec+(pts<np.array([0.,0.,0.]))*vec

def aind(x): 
	"""Convert back to advanced indexing."""
	return tuple(x.T)

def vecnorm(x): return x/np.linalg.norm(x)
def planeproject(x,n): return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)

def get_normal_fluctuation(hover,target,normal,vec):
	"""
	Given a target point, an instantaneous point, and the vertex normal, compute the distance to the tangent plane.
	"""
	vector = hover - target
	vector = vector - vec*(vector>(vec/2.)) + vec*(vector<(-1*vec/2.))
	projected = planeproject(vector,normal)
	#---get the sign of the projection
	plane_point = vector+projected
	sign = 1.0-2.0*(np.arccos(np.dot(vecnorm(normal),vecnorm(vector)))>np.pi/2.)
	return sign*np.linalg.norm(plane_point)

def inflate_lateral(source,inflate_factor):
	"""Add extra rows and columns to heights."""
	return source[np.meshgrid(*[np.arange(-inflate_factor,i+inflate_factor+1)%i for i in source.shape])]	

###---MAIN

if 'data' not in globals(): (data,calc),sns = plotload(plotname),work.sns()

do_inflate = False
sn = sns[0]
dat = data[sn]['data']
vecs = dat['vecs']
surfs = dat['mesh'].mean(axis=0)
#---average surface
surf_average_base = surfs.mean(axis=0)
#---inflate the average surface in XY
inflate_factor = 10
if do_inflate: surf_average = inflate_lateral(surf_average_base,inflate_factor)
else: surf_average = surf_average_base
#---height of the average surface
pivot = surf_average.mean()
#---standardized box vectors for all calculations (see notes above)
mvec_base = vecs.mean(axis=0)
#---get height fluctuations to set the half box height
maxflux = surfs.ptp()*1.1/2.
#---new standard box vectors have the correct height and inflated XY dimensions
inflate_factors = np.array(surf_average.shape).astype(float)/np.array(surf_average_base.shape)
if do_inflate: mvec = np.array([mvec_base[0]*inflate_factors[0],mvec_base[1]*inflate_factors[1],maxflux*2.])
else: mvec = np.array([mvec_base[0],mvec_base[1],maxflux*2.])
#---compute a reference surface in absolute points
#---we use vertical center so that all heights are shifted center of the new box given by twice maxflux
surf = boxstuff(height_recenter(literalize(surf_average,mvec),pivot=pivot,maxflux=maxflux),mvec)
#---make the reference mesh (slow step)
if 'mesh' not in globals(): 
	status('making mesh (somewhat slow)',tag='compute')1
	mesh = makemesh(surf,mvec)
	status('mesh is ready',tag='compute')

#---loop over frames
fr = 0

#---! getting: calcs/codes/mesh.py:24: RuntimeWarning: invalid value encountered in divide ... in vecnorm

#---inflate the instantaneous surface
this_surf_inflated = surfs[fr]#inflate_lateral(surfs[fr],inflate_factor)
#---find the points on the instantaneous surface which are nearest the points on the regular grid on the average
#---convert instantaneous points to XYZ with the reference box vectors mvec
instant_all = boxstuff(height_recenter(literalize(this_surf_inflated,mvec),pivot=pivot,maxflux=maxflux),mvec)
#---after literalizing the inflated points, we take only the points which are relevant to the base structure
#---! is the order correct?
if do_inflate:
	source = surf_average_base
	inds = np.concatenate(np.transpose(np.meshgrid(*[np.arange(-inflate_factor,i+inflate_factor+1) for i in source.shape])))
	base = np.where(np.all((np.all(inds>0,axis=1),np.all(np.array(source.shape)>=inds,axis=1)),axis=0))[0]
	instant = instant_all[base]
else: instant = instant_all

#---note that we make a tree from the instantaneous points then probe over the average surface
#---! more efficient to do this in reverse, however it might not cover all of the average/reference points?
#---prepare a KDTree. we use a fudge factor of 1000 epsilon to avoid angry errors about being outside the box
tree = scipy.spatial.ckdtree.cKDTree(instant,boxsize=np.concatenate((mvec,mvec))+1000.*eps)
#---find the nearest reference points for each instantaneous point
close,nns = tree.query(surf,k=1)
#---given a mapping between instantaneous point and target position (on XY), project the instantaneous point
#---...onto the tangent plane given by the reference point. note that this is obviously a minor approximation in which
#---...we assume that the point is hovering "above" the reference point close enough that the projection onto that 
#---...tangent plane is correct. a more literal form of this might actually try to find the exact distance to the 
#---...triangle adjacent to the nearest reference vertex, but this would require adding extra neighbor information
#---...and I think it takes the surface a bit too literally.
#---! note that we could use the real points instead of regular grid points for the instantaneous point?
deviations = np.array([
	get_normal_fluctuation(
		normal=mesh['vertnorms'][index],
		target=surf[index],
		hover=instant[nns][index],
		vec=mvec) 
	for ii,index in enumerate(nns)])
#---corners fail for whatever reason. could not get the do_inflate method working
deviations[np.isnan(deviations)] = 0.0

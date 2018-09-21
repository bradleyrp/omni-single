#!/usr/bin/env python

import MDAnalysis
import numpy as np
from numpy import linalg
import time
from joblib import Parallel,delayed
#from joblib.pool import has_shareable_memory
from base.tools import status,framelooper
#from codes.mesh import *
#from codes.head_tail import *

def head_angle(grofile,trajfile,**kwargs):

	"""
	Routine which computes the tilt and rotation angle of a PtdIns ring. Returns these angles
	along with the orientation of the tilt and a kernel density estimate comparison of tilt and rotation.
	"""

	#---unpack
	sn = kwargs['sn']
	work = kwargs['workspace']
	perform_kernel = kwargs.get('perform_kernel',False)
	
	#---prepare universe	
	uni = MDAnalysis.Universe(grofile,trajfile)
	nframes = len(uni.trajectory)
	#---MDAnalysis uses Angstroms not nm
	lenscale = 10.

	#---selections
	atomlist = ['C11','C14','C13','C15']
	all_resnames = work.vars['selectors']['resnames_PIP2']
	residue_list = [i for i in list(set(uni.residues.resnames)) if i in all_resnames]
	sels = [uni.select_atoms('resname '+work.meta[sn]['ptdins_resname']+' and name '+n)
		for n in atomlist]
	nmols = len(sels[0])

	#---preload all headgroup coordinates
	coords = np.zeros((len(atomlist),nframes,nmols,3))
	for fr in range(nframes):
		status('[LOAD] frame',i=fr,looplen=nframes)
		uni.trajectory[fr]
		for aa,a in enumerate(atomlist): 
			coords[aa,fr] = sels[aa].positions

	imono = kwargs['upstream']['lipid_abstractor']['monolayer_indices']
	#---assume resids are in the same order as incoming monolayer_indices, offset by 1
	#resids = sels[0].resids-1
	#---get the residue ids from the lipid_abstractor calculation
	resids = kwargs['upstream']['lipid_abstractor']['resids']
	imono_this = np.array([np.where(resids==i)[0][0] for i in sels[0].resids])
	leaflets = imono[imono_this]
	#---assume imono==0 is the "top" leaflet and set the normalvector accordingly
	#---! is this worth checking? 
	normalvec = np.zeros((len(imono_this),3))
	normalvec[np.where(leaflets==0),2] = 1
	normalvec[np.where(leaflets==1),2] = -1

	#---geometry functions
	cat = lambda a,b : np.concatenate((a,b))
	def vecnorm(vec): return vec/np.linalg.norm(vec)
	def planeproject(x,n): return x-np.dot(x,n)/np.linalg.norm(n)*vecnorm(n)
	def vecangle(v1,v2):
		"""
		Compute the anglebetween two vectors
		"""
		v1n,v2n = vecnorm(v1),vecnorm(v2)
		dotp = np.dot(v1n,v2n)
		angle = np.arccos(dotp)*(180./np.pi)	
		if np.isnan(angle): return (0.0 if (v1n==v2n).all() else np.pi*(180/np.pi))
		return angle

	#---compute ring vectors
	vec14 = coords[1]-coords[0]
	vec35 = coords[3]-coords[2]

	#---get orientations by projecting to the XY plane and checking for Z sign of 3,5-1,4 cross
	status('[COMPUTE] orientations')
	proj14 = np.concatenate((vec14[...,:2].T,[np.zeros((nmols,nframes))])).T
	proj35 = np.concatenate((vec35[...,:2].T,[np.zeros((nmols,nframes))])).T
	#---the orientations flag below tells you whether the head group is "tilted back"
	#---...which is to say that the cross product between the 1,4-3,5 vectors projected onto the 
	#---...XY plane has a negative Z coordinate, which can only happen if the ring is tilted in the 
	#---...opposite or "back" direction (compared to the usual forward one in which the left-handed
	#---...cross between 1,4-3,5 faces "up" indicating that the tilt is forward)
	orients = ((np.cross(proj14,proj35)[...,2]*np.tile(normalvec[:,2],(nframes,1)))<0.)
	#---we repeat this procedure to find the left-right orientation of the rotation
	proj_35_on_14 = np.array([[cat(planeproject(vec35[fr,n],vec14[fr,n])[:2],[0])
		for n in range(nmols)] 
		for fr in range(nframes)])
	proj_z_on_14 = np.array([[cat(planeproject(normalvec[n],vec14[fr,n])[:2],[0])
		for n in range(nmols)] 
		for fr in range(nframes)])			
	orients_rot = ((np.cross(proj_35_on_14,proj_z_on_14)[...,2]*np.tile(normalvec[:,2],(nframes,1)))<0.)
	
	#---compute angles
	#---note that we added normalvec multiplication below to flip the sign on the opposite bilayer
	#---...most likely because this 
	status('[COMPUTE] angles')
	angles14 = np.array([[vecangle(vec14[fr,n],normalvec[n]) 
		for n in range(nmols)] 
		for fr in range(nframes)])\
		*normalvec[:,2]*(1-2*orients)
	angles35 = np.array([[vecangle(vec35[fr,n],normalvec[n]) 
		for n in range(nmols)] 
		for fr in range(nframes)])\
		*normalvec[:,2]*(1-2*orients_rot)

	#---kernel density estimates
	if perform_kernel:
		try: bins = headerdat['calculations']['hta']['kde_bins']
		except: bins = 20.
		status('[COMPUTE] kernel density estimate bins = '+str(bins))
		data = vstack((array(reshape(angles14,-1)),array(reshape(angles35,-1)))).T
		kde = stats.kde.gaussian_kde(data.T)
		X, Y = mgrid[0:180:complex(0,bins), 0:180:complex(0,bins)]
		grid_coords = append(X.reshape(-1,1),Y.reshape(-1,1),axis=1)
		z = kde(grid_coords.T)
		z = z.reshape(bins,bins)

	#---package
	result = {}
	result['theta'] = np.array(angles14)
	result['phi'] = np.array(angles35)
	if perform_kernel:
		result['kernel_x'] = X
		result['kernel_y'] = Y
		result['kernel_z'] = z
	attrs = {}
	attrs['frames'] = nframes
	if perform_kernel: attrs['kernel_bins'] = bins
	return result,attrs

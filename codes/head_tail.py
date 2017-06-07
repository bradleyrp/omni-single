#!/usr/bin/python

from numpy.linalg import *
from scipy import stats
import numpy as np

def unit_vector(vector):
	""" Returns the unit vector of the vector.  """
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	""" Returns the angle in DEGREES (DRS) between vectors 'v1' and 'v2'::
	"""

	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	dotp = np.dot(v1_u, v2_u)
	cab = np.cross(v1_u,v2_u)
	angle = np.arccos(dotp)*(180./pi)	
	if np.isnan(angle):
		if (v1_u == v2_u).all():
			return 0.0
		else:
			return np.pi*(180/pi)
	return angle

def angle_between_ptdins(v1,v2):
	""" Returns the angle in DEGREES (DRS) between vectors 'v1' and 'v2'::
	"""

	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	dotp = np.dot(v1_u, v2_u)
	cab = np.cross(v1_u,v2_u)
	angle = np.arccos(dotp)*(180./pi)	
	if np.isnan(angle):
		if (v1_u == v2_u).all():
			return 0.0
		else:
			return np.pi*(180/pi)
	return angle

def generate_axes_pad(nrows=1,ncols=1,v_pad=0,h_pad=0,figsize=None):
	fig = plt.figure(figsize=figsize)
	gs = gridspec.GridSpec(nrows,ncols,hspace=v_pad,wspace=h_pad)
	axes = [[fig.add_subplot(gs[r,c]) for c in range(ncols)] for r in range(nrows)]
	return fig,gs,axes
	
def compute_hta(simname,grofile=None,trajfile=None,
	metadat=None,focus=None,panel=None,headerdat=None,**kwargs):

	"""
	This function returns the head group tilt and rotation angle relative 
	to the z axis of the simulation (not relative to a local membrane 
	normal or projection of the z axis onto the current residue). It 
	returns non-flattened arrays so that post-processing has access to the residue IDs.
	"""

	mset = SimSetMembrane()
	mset.load_trajectory((grofile,trajfile))
	status('Running hta.py using the global z axis for all dot products \
		and returning lipid-wise resolution.')
	z = asarray([0,0,1])
	c11c14_pi2p, c13c15_pi2p,p_pi2p = [], [], []
	theta, phi = [], []
	residue_list = []
	
	for fr in range(mset.nframes)[:]:
		status('Computing theta, phi',i=fr,looplen=mset.nframes)
		mset.gotoframe(fr) 
		c11c14_pi2p.append(mset.universe.selectAtoms('resname '+ \
			metadat[simname]['ptdins_resname'] + \
			' and (name C11 or name C14)').coordinates()) # Tilt, theta
		c13c15_pi2p.append(mset.universe.selectAtoms( 'resname '+ \
			metadat[simname]['ptdins_resname'] + \
			' and (name C13 or name C15)').coordinates()) # Rotation, phi
		p_pi2p.append(mset.universe.selectAtoms( 'resname '+ \
			metadat[simname]['ptdins_resname'] + \
			' and (name P)').coordinates()) # Rotation, phi
		residue_list.append(mset.universe.selectAtoms( 'resname '+ \
			metadat[simname]['ptdins_resname'] + \
			' and (name C13 or name C15)').resids())
		theta_by_frame = []
		phi_by_frame = []
		# Skip 2 because I'm selecting C11 and C14 above.
		for j in range(len(c11c14_pi2p[0]))[::2]:
			s = list(c11c14_pi2p[fr][j:j+2])
			c11c14 = transpose([[a-b] for a,b in zip(s[0],s[1])])[0]
			c11c14_hat = c11c14 / norm(c11c14)
			s = list(c13c15_pi2p[fr][j:j+2])
			c13c15 = transpose([[a-b] for a,b in zip(s[0],s[1])])[0]
			theta_by_frame.append(angle_between(c11c14,z))
			phi_by_frame.append(angle_between(c13c15,z))
		theta.append(theta_by_frame)	
		phi.append(phi_by_frame)
	
	# Keep flattened versions for compatibility
	all_theta = array([item for sublist in theta for item in sublist])
	all_phi = array([item for sublist in phi for item in sublist])
	print('Flattened theta, phi: '+str(len(all_theta))+' ' + str(len(all_phi)))
	bins = 64
	status('Computing kernel density estimate')
	data = vstack((array(all_theta),array(all_phi))).T
	kde = stats.kde.gaussian_kde(data.T)
	X, Y = np.mgrid[0:180:complex(0,bins), 0:180:complex(0,bins)]
	grid_coords = np.append(X.reshape(-1,1),Y.reshape(-1,1),axis=1)
	z = kde(grid_coords.T)
	z = z.reshape(bins,bins)
	result = {}
	result['theta'] = all_theta
	result['phi'] = all_phi
	result['theta_nonflat'] = array(neighborpack(theta))
	result['phi_nonflat'] = array(neighborpack(phi))

	result['kernel_x'] = X
	result['kernel_y'] = Y
	result['kernel_z'] = z
	result['resids'] = residue_list
	
	if 'times' in kwargs: result['times'] = array(kwargs['times'])

	attrs = {}
	attrs['frames'] = mset.nframes
	attrs['kernel_bins'] = bins
	return result,attrs
	
def head_tail(grofile,trajfile,headerdat=None,**kwargs):

	"""
	Routine which computes the tilt and rotation angle of a PtdIns ring. Returns these angles
	along with the orientation of the tilt and a kernel density estimate comparison of tilt and rotation.
	"""
	
	#---prepare
	sn = simname
	mset = SimSetMembrane()
	mset.load_trajectory((grofile,trajfile))
	mset.identify_monolayers(metadat[sn]['director'])
	residues = [i for i in list(set(mset.universe.residues.resnames())) if i in all_resnames]
	mset.identify_residues(residues)

	resid_by_monolayer = [sort(j) for j in 
		[mset.monolayer_by_resid[m][residues.index(metadat[sn]['ptdins_resname'])] 
		for m in range(2)]]
	#---zero means upper monolayer
	imono = array([0 if i in resid_by_monolayer[0] else 1 for i in concatenate(resid_by_monolayer)])
	
	#---selections
	atomlist = ['C11','C14','C13','C15']
	residue_list = [i for i in list(set(mset.universe.residues.resnames())) if i in all_resnames]
	sels = [mset.universe.selectAtoms(
		'resname '+metadat[simname]['ptdins_resname']+' and name '+n)
		for n in atomlist]
	nmols = len(sels[0])
	nframes = mset.nframes

	#---preload all coordinates
	coords = zeros((len(atomlist),nframes,nmols,3))
	for fr in range(mset.nframes):
		status('[LOAD] frame',i=fr,looplen=mset.nframes)
		mset.gotoframe(fr)
		for aa,a in enumerate(atomlist): 
			coords[aa,fr] = sels[aa].coordinates()

	#---geometry functions
	cat = lambda a,b : concatenate((a,b))
	def vecnorm(vec): return vec/linalg.norm(vec)
	def planeproject(x,n): return x-dot(x,n)/linalg.norm(n)*vecnorm(n)
	def vecangle(v1,v2):
		"""
		Compute the anglebetween two vectors
		"""
		v1n,v2n = vecnorm(v1),vecnorm(v2)
		dotp = dot(v1n,v2n)
		angle = arccos(dotp)*(180./pi)	
		if isnan(angle): return (0.0 if (v1n==v2n).all() else pi*(180/pi))
		return angle

	#---compute ring vectors
	vec14 = coords[1]-coords[0]
	vec35 = coords[3]-coords[2]
	
	#---determine normal vector depending on leaflet
	normalvec = array([array([0,0,1]) if j in resid_by_monolayer[0] else array([0,0,-1]) 
		for j in sels[0].resids()-1])

	#---get orientations by projecting to the XY plane and checking for Z sign of 3,5-1,4 cross
	status('[COMPUTE] orientations')
	proj14 = concatenate((vec14[...,:2].T,[zeros((nmols,nframes))])).T
	proj35 = concatenate((vec35[...,:2].T,[zeros((nmols,nframes))])).T
	#---the orientations flag below tells you whether the head group is "tilted back"
	#---...which is to say that the cross product between the 1,4-3,5 vectors projected onto the 
	#---...XY plane has a negative Z coordinate, which can only happen if the ring is tilted in the 
	#---...opposite or "back" direction (compared to the usual forward one in which the left-handed
	#---...cross between 1,4-3,5 faces "up" indicating that the tilt is forward)
	orients = ((cross(proj14,proj35)[...,2]*tile(normalvec[:,2],(nframes,1)))<0.)
	#---we repeat this procedure to find the left-right orientation of the rotation
	proj_35_on_14 = array([[cat(planeproject(vec35[fr,n],vec14[fr,n])[:2],[0])
		for n in range(nmols)] 
		for fr in range(nframes)])
	proj_z_on_14 = array([[cat(planeproject(normalvec[n],vec14[fr,n])[:2],[0])
		for n in range(nmols)] 
		for fr in range(nframes)])			
	orients_rot = ((cross(proj_35_on_14,proj_z_on_14)[...,2]*tile(normalvec[:,2],(nframes,1)))<0.)
	
	#---compute angles
	#---note that we added normalvec multiplication below to flip the sign on the opposite bilayer
	#---...most likely because this 
	status('[COMPUTE] angles')
	angles14 = array([[vecangle(vec14[fr,n],normalvec[n]) 
		for n in range(nmols)] 
		for fr in range(nframes)])\
		*normalvec[:,2]*(1-2*orients)
	angles35 = array([[vecangle(vec35[fr,n],normalvec[n]) 
		for n in range(nmols)] 
		for fr in range(nframes)])\
		*normalvec[:,2]*(1-2*orients_rot)

	#---kernel density estimates
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
	result['theta'] = array(angles14)
	result['phi'] = array(angles35)
	result['kernel_x'] = X
	result['kernel_y'] = Y
	result['kernel_z'] = z
	result['resids0'] = array(resid_by_monolayer[0])
	result['resids1'] = array(resid_by_monolayer[1])
	if 'times' in kwargs: result['times'] = array(kwargs['times'])
	attrs = {}
	attrs['frames'] = mset.nframes
	attrs['kernel_bins'] = bins
	return result,attrs

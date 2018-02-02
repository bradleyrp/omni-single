#!/usr/bin/env python

"""
MESH ROUTINES

Useful for making triangulated surfaces.
"""

import numpy as np
import scipy
import scipy.spatial
import scipy.interpolate
from numpy import linalg
machine_eps = eps = np.finfo(float).eps

#---import from omni.base
from base.tools import status
from base.timer import time_limit,TimeoutException
from base.compute_loop import basic_compute_loop

_not_reported = ['triarea','vecnorm','facenorm','torusnorm','reclock','beyonder','makemesh','rotation_matrix']
_shared_extensions = ['vecnorm','rotation_matrix','makemesh']

#---geometry functions
triarea = lambda a : linalg.norm(np.cross(a[1]-a[0],a[2]-a[0]))/2.
vecnorm = lambda vec: np.array(vec)/linalg.norm(vec)
facenorm = lambda a: np.cross(a[1]-a[0],a[2]-a[0])

def centroid(coords,masses,divider):
	"""
	Compute the centroid of a collection of molecules given their atomistic XYZ coordinates, 
	their masses, and a list of lists containing the indices for the atoms in each molecule.
	This function refers to masses and divider from the global namespace. 
	Note: using the global namespace improves speed from 2.3s to 1.2s for a test system.
	"""
	return np.array([np.dot(masses[r].T,coords[r])/np.sum(masses[r]) for r in divider])

def rotation_matrix(axis,theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians using Euler-Rodrigues formula.
	"""
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	if np.all(axis==0): return np.identity(3) 
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b, c, d = -axis*np.sin(theta/2)
	aa, bb, cc, dd = a*a, b*b, c*c, d*d
	bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
	return np.array([[aa+bb-cc-dd,2*(bc+ad),2*(bd-ac)],[2*(bc-ad),aa+cc-bb-dd,2*(cd+ab)],
		[2*(bd+ac),2*(cd-ab),aa+dd-bb-cc]])

def reclock(points,ords=[0,1,2]):
	"""
	Reorders R2 points in clockwise order.
	"""
	rels = points[ords]-np.mean(points[ords],axis=0)
	return np.argsort(np.arctan2(rels[:,0],rels[:,1]))[::-1]

def beyonder(points,vecs,dims=[0,1],growsize=0.2,new_only=False,growsize_nm=None,return_ids=False):
	"""
	Given points and box vectors, this function will expand a bilayer patch across periodic boundaries. It is 
	meant to create inputs for meshing algorithms which later store data under the correct topology. The dims
	argument names the wrapping dimensions, which is typically XY.
	"""
	over = np.array([(i,j,k) for i in [-1,0,1] for j in [-1,0,1] for k in [-1,0,1]])
	over = np.array([list(j) for j in list(set([tuple(i) for i in over[:,:len(dims)]])) 
		if j != tuple(np.zeros(len(dims)))])
	vec2 = np.array([(vecs[i] if i in dims else 0.) for i in range(3)])	
	over2 = np.array([[(o[i] if i in dims else 0.) for i in range(3)] for o in over])
	alls = np.concatenate([points]+[points+i*vec2 for i in over2])
	#---save indices for extra points
	inds = np.concatenate([np.arange(len(points)) for i in range(len(over2)+1)])
	#---bugfix applied 2015.07.20 to fix incorrect ghost index numbering
	if 0: inds = np.concatenate([np.arange(len(points))]+[np.ones(len(over))*j 
		for j in np.arange(len(points))])
	#---override the proportional growsize by a literal size
	if growsize_nm != None: growsize = max([growsize_nm/vecs[i] for i in dims])
	valids = np.where(np.all([np.all((alls[:,d]>-1*vec2[d]*growsize,alls[:,d]<vec2[d]*(1.+growsize)),axis=0) 
		for d in dims],axis=0))[0]
	if not return_ids: return alls[valids]
	else: return alls[valids],inds[valids].astype(int)

def torusnorm(pts1,pts2,vecs):
	"""
	Compute distances between points on a torus.
	"""
	cd = np.array([scipy.spatial.distance.cdist(pts1[:,d:d+1],pts2[:,d:d+1]) for d in range(2)])
	cd[0] -= (cd[0]>vecs[0]/2.)*vecs[0]
	cd[1] -= (cd[1]>vecs[1]/2.)*vecs[1]
	cd2 = linalg.norm(cd,axis=0)
	return cd2

def makemesh(pts,vec,growsize=0.2,curvilinear_neighbors=10,
	curvilinear=True,debug=False,growsize_nm=None,excise=True,areas_only=False):
	"""
	Function which computes curvature and simplex areas on a standard mesh.
	"""
	if debug: 
		import time
		st = time.time()
		def checkpoint(note):
			print(note)
			st = time.time()

	nmol = len(pts)
	pts = pts
	vec = vec
	if debug: 
		import time
		st = time.time()
		print("[STATUS] start makemesh %0.2f"%(time.time()-st))
	ptsb,ptsb_inds = beyonder(pts,vec,growsize=growsize,growsize_nm=growsize_nm,return_ids=True)
	if debug: print("[STATUS] project curvilinear="+str(curvilinear)+" %0.2f"%(time.time()-st))
	#---if curvilinear then use the isomap otherwise project onto the xy plane
	if curvilinear: 
		import sklearn
		from sklearn import manifold
		proj = manifold.Isomap(curvilinear_neighbors,2).fit_transform(ptsb)
	else: proj = ptsb[...,:2]
	if debug: checkpoint("[STATUS] delaunay %0.2f"%(time.time()-st))
	if debug: checkpoint("[STATUS] shape="+str(np.shape(ptsb)))
	dl = scipy.spatial.Delaunay(proj)
	if debug: checkpoint("[STATUS] reclock %0.2f"%(time.time()-st))
	simplices = np.array([a[reclock(ptsb[a])] for a in dl.simplices])
	#---rework simplices and ptsb to exclude superfluous points
	if debug: checkpoint("[STATUS] trim %0.2f"%(time.time()-st))
	#---relevants is a unique list of simplices with exactly one member that is equal to a core vertex point
	relevants = np.unique(np.concatenate([simplices[
		np.where(np.sum(simplices==i,axis=1)==1)[0]] for i in range(nmol)]))
	points = ptsb[relevants]
	ghost_indices = ptsb_inds[relevants]
	ptsb = points

	if debug: checkpoint("[STATUS] simplices %0.2f"%(time.time()-st))
	simplices = np.array([[np.where(relevants==r)[0][0] for r in s] 
		for s in simplices if np.all([r in relevants for r in s])])
	#---end rework
	if debug: checkpoint("[STATUS] areas %0.2f"%(time.time()-st))
	areas = np.array([triarea(ptsb[a]) for a in simplices])
	if areas_only: return {'simplices':simplices,'areas':areas,'nmol':nmol,'vec':vec,'points':points}
	if debug: checkpoint("[STATUS] facenorms %0.2f"%(time.time()-st))
	facenorms = np.array([vecnorm(facenorm(ptsb[a])) for a in simplices])	
	if debug: checkpoint("[STATUS] vertex-to-simplex %0.2f"%(time.time()-st))
	v2s = [np.where(np.any(simplices==i,axis=1))[0] for i in range(nmol)]
	if debug: checkpoint("[STATUS] vertex normals %0.2f"%(time.time()-st))
	vertnorms = np.array([vecnorm(np.sum(facenorms[ind]*\
		np.transpose([areas[ind]/np.sum(areas[ind])]),axis=0)) for ind in v2s])
	principals = np.zeros((nmol,2))
	nl = []
	if debug: checkpoint("[STATUS] curvatures %0.2f"%(time.time()-st))
	for v in range(nmol):
		neighbors = np.unique(simplices[np.where(np.any(simplices==v,axis=1))[0]])
		neighbors = neighbors[neighbors!=v]
		nl.append(neighbors)
		edges = ptsb[neighbors]-ptsb[v]
		weights = [areas[sl]/2./np.sum(areas[v2s[v]]) for sl in v2s[v]]
		tijs = [vecnorm(np.dot(np.identity(3)-np.outer(vertnorms[v],
			vertnorms[v].T),ab)) for ab in edges]
		kijs = [np.dot(vertnorms[v].T,ab)/linalg.norm(ab)**2 for ab in edges]
		ct = np.sum([weights[ind]*kijs[ind]*np.outer(tijs[ind],tijs[ind]) 
			for ind,i in enumerate(v2s[v])],axis=0)
		wsign = 1-2*(linalg.norm(np.array([1,0,0])+\
			vertnorms[v])<linalg.norm(np.array([1,0,0])-vertnorms[v]))
		wvi = vecnorm(np.array([1,0,0])+wsign*vertnorms[v])
		hm = np.identity(3)-2*np.outer(wvi,wvi.T)
		hhm = np.dot(np.dot(hm.T,ct),hm)
		principals[v] = -1*hhm[1,1],-1*hhm[2,2]
	if debug: checkpoint("[STATUS] PBC neighborlist %0.2f"%(time.time()-st))
	#---neighborlist under PBCs
	checksubssort,nlsubs = np.where(torusnorm(points[nmol:],points[:nmol],vec)==0)
	#if not all(checksubssort==np.arange(len(points)-nmol)): raise Exception('torusnorm lookup fail')
	try: nlpbc = [[(i if i<nmol else nlsubs[i-nmol]) for i in n] for n in nl]
	except: nlpbc = []
	gauss = (3*principals[:,0]-principals[:,1])*(3*principals[:,1]-\
		principals[:,0])
	mean = 1./2*((3*principals[:,0]-principals[:,1])+\
		(3*principals[:,1]-principals[:,0]))
	if debug: checkpoint("[STATUS] complete %0.2f"%(time.time()-st))

	if debug:
		import matplotlib as mpl;import matplotlib.pylab as plt
		plt.scatter(points[:,0],points[:,1])
		plt.show()
		import pdb;pdb.set_trace()

	return {'nmol':nmol,'vec':vec,'simplices':simplices,'points':points,
		'areas':areas,'facenorms':facenorms,'vertnorms':vertnorms,'principals':principals,
		'ghost_ids':ghost_indices,'gauss':gauss,'mean':mean}

def identify_lipid_leaflets_legacy(pts,vec,monolayer_cutoff,
	monolayer_cutoff_retry=True,max_count_asymmetry=0.05,pbc_rewrap=True,
	topologize_tolerance=None,topologize_time_limit=30):
	"""
	Identify leaflets in a bilayer by consensus.
	Note that the time limit on the topologize call was increased from 10 to 30 for large systems.
	This is the legacy version of this algorithm. Previously it was recursive, lowering the cutoff by small
	increments and then calling itself again if the bilayer did not appear to be split correctly. The current
	version is called by the LeafletFinder class and throws exceptions to trigger a lower cutoff. We have
	tried to preserve the legacy version for other users, but the cluster version is more reliable.
	"""
	#---previous default was somewhat high, but typically came in from specs, and we reduced it incrementally
	if monolayer_cutoff==None: monolayer_cutoff = 2.0
	#---time limit on the tolerance checker
	try:
		with time_limit(topologize_time_limit): 
			wrapper = topologize(pts,vec,
				**({'tol':topologize_tolerance} if topologize_tolerance else {}))
	except TimeoutException, msg: 
		status('topologize failed to join the bilayer. '
			'if it is broken over PBCs e.g. a saddle, this is a serious error which may go undetected. '
			'make sure you always inspect the topology later.',tag='error')
		wrapper = np.zeros((len(pts),3))
	findframe = pts + wrapper*np.array(vec)
	status('this step is somewhat slow. it uses scipy.spatial.pdist.',tag='warning')
	pd = [scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(findframe[:,d:d+1])) 
		for d in range(3)]
	if pbc_rewrap:
		pd3pbc = np.sqrt(np.sum(np.array([pd[d]-(pd[d]>vec[d]/2.)*vec[d]+(pd[d]<-1*vec[d]/2.)*vec[d] 
			for d in range(3)])**2,axis=0))
	else: pd3pbc = pd
	nbors = np.transpose(np.where(pd3pbc<monolayer_cutoff))
	nlipids = len(pts)
	imono = np.zeros(nlipids)
	nlist = []
	for i in range(nlipids):
		status('cataloging lipids',i=i,looplen=nlipids,tag='compute')
		nlist.append(nbors[np.where(nbors[:,0]==i)[0],1])
	iref = 0
	mono = np.zeros(nlipids)
	searched = np.zeros(nlipids)
	imono[iref],searched[iref] = 1,1
	imono[nlist[iref]] = 1
	while np.any(np.all((imono==1,searched==0),axis=0)):
		for iref in np.where(np.all((imono==1,searched==0),axis=0))[0]: 
			imono[nlist[iref]] = 1
			searched[iref] = 1
	#---check that the leaflets were properly distinguished by looking at the number in each monolayer
	if np.mean(imono)==0.5: 
		status('[STATUS] perfect split is %0.5f'%np.mean(imono))
		return imono
	elif (monolayer_cutoff_retry and (np.all(np.array(imono)==0) or np.all(np.array(imono)==1) or 
		np.abs(np.mean(imono)-0.5)>=max_count_asymmetry)):
		status('[STATUS] split is %0.5f'%np.mean(imono))
		status('[STATUS] one side has %d'%np.sum(imono))
		status('[WARNING] leaflets were not distinguished')
		status('[COMPUTE] leaflets = '+str(np.sum(imono))+'/'+str(len(imono)))
		status('[WARNING] previous monolayer_cutoff = '+str(monolayer_cutoff))
		raise Exception(
			'[ERROR] failed to identify leaflets so we are returning an exception to the LeafletFinder')
	else: status('[STATUS] some lipids might be flipped %d %.5f'%(np.sum(imono),np.mean(imono)))
	return imono

def topologize(pos,vecs,tol=None):
	"""
	Join a bilayer which is broken over periodic boundary conditions by translating each point by units
	of length equal to the box vectors so that there is a maximum amount of connectivity between adjacent
	points. This method decides how to move the points by a consensus heuristic.
	This function is necessary only if the bilayer breaks over a third spatial dimension.
	Note that we changed tol from 0.05 to 0.07 on 2017.08.02. We also added a timer and pass-through for 
	belligerent systems. See RPB notes for more details.
	"""
	if tol==None: tol = 0.07
	step = 0
	kp = np.array(pos)
	natoms = len(pos)
	move_votes = np.zeros((1,natoms,3))
	while np.sum(np.abs(move_votes[0])>len(pos)/2)>len(pos)*tol or step == 0:
		move_votes = np.concatenate((move_votes,np.zeros((1,natoms,3))))
		pos = np.array(kp)
		#---! can the pdist be optimized somehow?
		pd = [scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(pos[:,d:d+1])) for d in range(3)]
		for ind,bead in enumerate(pos):
			#---compute distances to the probe
			for d in range(3):
				adds_high = np.zeros(natoms)
				adds_low = np.zeros(natoms)
				adds_high[np.where(np.all((pd[d][ind]>vecs[d]/2.,pos[ind,d]<pos[:,d]),axis=0))] = 1
				adds_low[np.where(np.all((pd[d][ind]>vecs[d]/2.,pos[ind,d]>pos[:,d]),axis=0))] = -1
				#---for each spatial dimension, we tally votes to move the point by a box vector
				move_votes[step][:,d] += adds_high+adds_low
		kp = np.array(pos)
		step += 1
		move_votes = np.concatenate((move_votes,np.zeros((1,natoms,3))))
	moved = np.transpose([np.sum([1*(move_votes[it][:,d]<(-1*natoms/2.))-1*(move_votes[it][:,d]>(natoms/2.)) 
		for it in range(step+1)],axis=0) for d in range(3)])
	return moved

class LeafletFinder:

	"""
	Manage different algorithms for detecting the separate leaflets.
	"""

	def __init__(self,atoms_separator,vecs,**kwargs):
		#---we require the separator and vectors
		self.atoms_separator = atoms_separator
		#---nframes is the number of separate frames to use to attempt the separator
		self.vecs,self.nframes = vecs,len(vecs)
		#---new flags
		self.cluster = kwargs.pop('cluster',False)
		self.scan_mode = kwargs.pop('scan_mode',False)
		#---legacy flags
		self.monolayer_cutoff = kwargs.pop('monolayer_cutoff',None)
		self.monolayer_cutoff_retry = kwargs.pop('monolayer_cutoff_retry',True)
		self.topologize_tolerance = kwargs.pop('topologize_tolerance',None)
		self.cutoff_shrink_increment = kwargs.pop('cutoff_shrink_increment',None)
		self.cutoff_min = kwargs.pop('cutoff_min',None)
		self.random_tries = kwargs.pop('random_tries',None)
		self.cluster_neighbors = kwargs.pop('cluster_neighbors',None)
		if self.cluster_neighbors==None: self.cluster_neighbors = 4
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		#---check for scikit-learn
		if self.cluster:
			try: import sklearn
			except: 
				status('cannot import scikit-learn so we will use legacy leaflet finder',tag='warning')
				self.cluster = False
		#---the persistent function tries to distinguish leaflets according to the mode
		self.persistent()

	def persistent(self):
		"""
		Try to find the leaflets by using multiple frames and multiple cutoffs.
		"""
		if self.monolayer_cutoff==None: self.monolayer_cutoff = 2.0
		#---determine the mode and retry settings
		if self.cutoff_shrink_increment==None: self.cutoff_shrink_increment = 0.01
		#---previously we reduced the cutoff to zero before trying a different frame
		if self.cutoff_min==None: self.cutoff_min = 0.8
		#---legacy mode
		if not self.cluster:
			#---try multiple times
			if self.monolayer_cutoff_retry:
				#---legacy retry mode starts high and reduces the cutoff at each step
				#---! we could implement a method that tries cutoffs above/below the start point
				cutoffs = np.arange(self.cutoff_min,
					self.monolayer_cutoff+self.cutoff_shrink_increment,self.cutoff_shrink_increment)[::-1]
			#---only try one cutoff
			else: cutoffs = [self.monolayer_cutoff]
		#---cluster mode uses a default cutoff
		else: cutoffs = [None]
		monolayer_indices = None
		#---recall that the caller provides frames for testing
		for fr in range(self.nframes):
			#---loop over cutoffs if we have multiple cutoffs
			for cutoff in cutoffs:
				if not self.cluster:
					try: 
						if not self.cluster:
							#---call the legacy leaflet finder (outside of this class)
							monolayer_indices = identify_lipid_leaflets_legacy(
								self.atoms_separator[fr],self.vecs[fr],monolayer_cutoff=cutoff)
					except: status('failed to distinguish leaflets with cluster=%s and cutoff=%s'%(
						self.cluster,cutoff),tag='error')
				else: monolayer_indices = self.identify_leaflets_cluster(
					pts=self.atoms_separator[fr],vec=self.vecs[fr])
				#---break when successful
				if type(monolayer_indices)!=bool: 
					self.monolayer_indices = monolayer_indices
					return

	def identify_leaflets_cluster(self,pts,vec,topologize_time_limit=30,max_count_asymmetry=0.05):
		"""
		Use scikit-learn clustering methods to separate leaflets.
		Note that this method can cluster a tortuous manifold and may work for complex morphologies.	
		"""
		import scipy
		import sklearn
		import sklearn.neighbors
		import sklearn.cluster
		nlipids = len(pts)
		#---time limit on the topologize function which joins broken bilayers e.g. a saddle that crosses PBCs
		try:
			with time_limit(topologize_time_limit): 
				wrapper = topologize(pts,vec,
					**({'tol':self.topologize_tolerance} if self.topologize_tolerance else {}))
		except TimeoutException, msg: 
			status('topologize failed to join the bilayer. '
				'if it is broken over PBCs e.g. a saddle, this is a serious error which may go undetected. '
				'make sure you always inspect the topology later.',tag='error')
			wrapper = np.zeros((len(pts),3))
		findframe = pts + wrapper*np.array(vec)
		#---ensure that all points are in the box
		findframe += vec*(findframe<0) - vec*(findframe>vec)
		#---previous calculation of connectivity was done manually
		if False:
			#---conservative cutoff gets lots of nearby points
			cutoff = 10.0
			cutoff_short = 2.0
			#---make a K-D tree from the points
			tree = scipy.spatial.ckdtree.cKDTree(findframe,boxsize=np.concatenate((vec,vec))+0.*eps)
			#---find the nearest reference points for each instantaneous point
			close,nns = tree.query(findframe,distance_upper_bound=cutoff,k=20)
			#---construct the neighbor list
			subjects = np.where(np.all((close<cutoff,close>0),axis=0))
			#---get the pairs of neighbors
			subjects,neighbors = subjects[0],nns[subjects]
			pds = np.ones((nlipids,nlipids))*0.0
			pds[tuple((np.arange(nlipids),np.arange(nlipids)))] = 0.0
			nears = np.where(np.all((close>0,close<=cutoff_short),axis=0))
			pds[tuple((nears[0],nns[nears]))] = 1.0#close[nears]
			pds[tuple((nns[nears],nears[0]))] = 1.0#close[nears]
		connectivity = sklearn.neighbors.kneighbors_graph(findframe,
			n_neighbors=self.cluster_neighbors,include_self=False)
		ward = sklearn.cluster.AgglomerativeClustering(n_clusters=2,
			connectivity=connectivity,linkage='complete').fit(findframe)
		imono = ward.labels_
		if np.mean(imono)==0.5: 
			status('[STATUS] perfect split is %0.5f'%np.mean(imono))
		elif (np.all(np.array(imono)==0) or np.all(np.array(imono)==1) or 
			np.abs(np.mean(imono)-0.5)>=max_count_asymmetry):
			status('[STATUS] split is %0.5f'%np.mean(imono))
			status('[STATUS] one side has %d'%np.sum(imono))
			status('[WARNING] leaflets were not distinguished')
			raise Exception('[ERROR] failed to identify leaflets. '
				'DEVELOPMENT NOTE!? use legacy or a different cutoff?')
		else: status('[STATUS] some lipids might be flipped %d %.5f'%(np.sum(imono),np.mean(imono)))
		return np.array(imono)

def makemesh_regular(data,vecs,grid):
	"""
	Generate a regular grid from a monolayer.
	"""
	data = beyonder(data,vecs,growsize=0.1)
	xypts = np.array([[i,j] for i in np.linspace(0,vecs[0],grid[0].astype(int)) 
		for j in np.linspace(0,vecs[1],grid[1].astype(int))])
	interp = scipy.interpolate.LinearNDInterpolator(data[:,0:2],data[:,2],fill_value=0.0)
	bilinear_pts = np.array([[i[0],i[1],interp(i[0],i[1])] for i in xypts])
	result = scipy.interpolate.griddata(bilinear_pts[:,0:2],bilinear_pts[:,2],bilinear_pts[:,0:2],
		method='cubic')
	#---observed that griddata returns points where we cycle through the points in the following
	#---...order:x0,y0),(x0,y1),...(x0,yn),(x1,y0),... and so on, suggesting that the following 
	#---...reshape command (which reshape function claims to use the "C" programming language convention
	#---...for reshaping objects by default, which convention has the last index changing "fastest")
	xyz_pts = np.array([[bilinear_pts[i,0],bilinear_pts[i,1],result[i]] for i in range(len(result))])
	return np.reshape(xyz_pts[:,2],grid.astype(int))

###---AVERAGE  NORMAL PROJECTION FUNCTIONS

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

def average_normal_projections(fr,mvec,pivot,maxflux,do_inflate=False):
	"""
	Projection subroutine for measure_normal_deviation_from_wavy_surface moved here for computation in parallel.
	"""
	global surf,surfs,mesh
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
		inds = np.concatenate(np.transpose(np.meshgrid(*[np.arange(-inflate_factor,i+inflate_factor+1) 
			for i in source.shape])))
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
	#---...onto the tangent plane given by the reference point. note that this is obviously a minor approximation in 
	#---...which we assume that the point is hovering "above" the reference point close enough that the projection onto
	#---...that tangent plane is correct. a more literal form of this might actually try to find the exact distance to 
	#---...the triangle adjacent to the nearest reference vertex, but this would require adding extra neighbor
	#---...information and I think it takes the surface a bit too literally.
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
	return deviations

def measure_normal_deviation_from_wavy_surface(heights,vecs,curvilinear=False):
	"""
	Given heights on a regular grid, compute the average surface and then compute the 
	"""
	global surf,surfs,mesh
	do_inflate = False
	inflate_factor = 10
	surfs = heights
	#---average surface
	surf_average_base = surfs.mean(axis=0)
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
	#---use globals for parallel
	if do_inflate: mvec = np.array([mvec_base[0]*inflate_factors[0],mvec_base[1]*inflate_factors[1],maxflux*2.])
	else: mvec = np.array([mvec_base[0],mvec_base[1],maxflux*2.])
	#---compute a reference surface in absolute points
	#---we use vertical center so that all heights are shifted center of the new box given by twice maxflux
	surf = boxstuff(height_recenter(literalize(surf_average,mvec),pivot=pivot,maxflux=maxflux),mvec)
	#---make the reference mesh (slow step)
	status('making mesh (curvilinear=%s)'%curvilinear,tag='compute')
	mesh = makemesh(surf,mvec,curvilinear=curvilinear)
	status('mesh is ready',tag='compute')
	looper = [dict(fr=fr,pivot=pivot,mvec=mvec,maxflux=maxflux) for fr in range(len(surfs))]
	incoming = basic_compute_loop(average_normal_projections,looper=looper,run_parallel=True)
	#---we must reshape and concatenate the points
	return np.reshape(incoming,(-1,)+surf_average.shape)

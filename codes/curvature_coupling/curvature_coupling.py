#!/usr/bin/env python

"""
CURVATURE-UNDULATION COUPLING (version 3b)
"""

import importlib,re
import numpy as np
from base.tools import status
from base.compute_loop import basic_compute_loop
import codes.curvature_coupling.tools as cctools
import scipy
import scipy.optimize
from codes.undulate import perfect_collapser,blurry_binner

machine_eps = eps = np.finfo(float).eps

###---AT-LARGE FUNCTIONS

def make_fields(**kwargs):
	"""Wrap the curvature field constructor for joblib."""
	return cctools.construct_curvature_field(**kwargs)

dotplace = lambda n : re.compile(r'(\d)0+$').sub(r'\1',"%3.5f"%float(n)).rjust(8)

def vecnorm(vec): return vec/np.linalg.norm(vec)

def principal_axis(pts):
	"""
	Return one of the three principal axes of a three-dimensional collection of points.
	! Put this in a central location?
	"""
	axis = vecnorm(pts[0]-pts[-1])
	eigs = np.linalg.eig(np.dot(pts.T,pts))
	principal_axis_index = np.argsort(eigs[0])[-1]
	axis = vecnorm(eigs[1][:,principal_axis_index])
	return axis

def rotation_matrix(axis,theta):
	"""
	Return the rotation matrix associated with counterclockwise rotation about
	the given axis by theta radians using Euler-Rodrigues formula.
	! Put this in a central location?
	"""
	axis = np.asarray(axis)
	theta = np.asarray(theta)
	if all(axis==0): return np.identity(3) 
	axis = axis/np.sqrt(np.dot(axis,axis))
	a = np.cos(theta/2)
	b,c,d = -axis*np.sin(theta/2)
	aa,bb,cc,dd = a*a,b*b,c*c,d*d
	bc,ad,ac,ab,bd,cd = b*c,a*d,a*c,a*b,b*d,c*d
	return np.array([[aa+bb-cc-dd,2*(bc+ad),2*(bd-ac)],[2*(bc-ad),aa+cc-bb-dd,2*(cd+ab)],
		[2*(bd+ac),2*(cd-ab),aa+dd-bb-cc]])

###---OBJECTIVE FUNCTION and associated paraphanalia

def gopher(spec,module_name,variable_name,work=None):
	"""Load an external module. Useful for changing the workflow without changing the code."""
	#---note that this function was folded into the omnicalc codebase in omni/base/tools.py
	mod = importlib.import_module(spec[module_name])
	target = mod.__dict__.get(spec[variable_name],None)
	#---the database design might need work so we always export it
	mod.work = work
	if not target: raise Exception('add %s and %s to the specs'%(module_name,variable_name))
	return target

def formulate_wavevectors(dims,vecs=None,lx=None,ly=None,lenscale=1.0):
	"""
	Generate wavevectors from box dimensions and the number of grid points.
	"""
	if type(vecs)!=type(None):
		Lx,Ly = np.mean(vecs,axis=0)[:2]
	elif type(lx)==type(None) or type(ly)==type(None):
		raise Exception('send box dimensions (lx and ly) or framewise box vectors (vecs)')
	else: Lx,Ly = lx,ly
	if len(dims)!=2: raise Exception('dims must be the number of grid points in each of two directions')
	#---formulate the wavevectors
	lenscale = 1.0
	m,n = mn = dims
	q2d = lenscale*np.array([[np.sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	q_raw = np.reshape(q2d,-1)[1:]
	area = (Lx*Ly/lenscale**2)
	return dict(wavevectors=q_raw,area=area)

def curvature_sum_function(cfs,curvatures,**kwargs):
	"""
	Curvature fields are maxed not summed.
	!? What's with the weird transpose below. Document this!!!
	!? Should this be turned into a decorator?
	"""
	method = kwargs.get('method')
	if method=='mean':
		combo = np.array([np.transpose(cfs,(1,0,2,3))[cc]*c 
			for cc,c in enumerate(curvatures)]).mean(axis=0)
	else: raise Exception('invalid method %s'%method)
	return combo

def prepare_residual_DEPRECATED(residual_form='log',weighting_scheme=None):
	"""
	DEPRECATED!
	Decorate a residual function.
	"""
	#---residual mode for the return format from the objective
	if residual_form == 'log':
		if weighting_scheme=='blurry_explicit':
			def residual(values):
				return np.sum(weights*np.log10(values.clip(min=machine_eps))**2)/float(len(values))
		else: 
			def residual(values):
				return np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
	#---deprecated residual method
	elif residual_form == 'linear':
		raise Exception('linear is not the right residual option') 
		def residual(values): 
			return np.sum((values-1.0)**2)/float(len(values))
	else: raise Exception('unclear residual form %s'%residual_form)
	return residual

def prepare_oscillator_function(reverse=False):
	"""Reversible for .... !!! mysterious reasons."""
	def oscillator_function(vibe,qs):
		"""Standard harmonic oscillator function."""
		#---vibration must be nonnegative
		vibe = np.abs(vibe)
		raw = ((vibe*qs+machine_eps)/((np.exp(vibe*qs)-1)+machine_eps))
		if not reverse: return raw
		else: return 10.**(-1.0*np.log10(raw))
	return oscillator_function

def prepare_residual():
	"""Decorate the residual function."""
	#def residual(hel,hosc): return np.mean((np.log10(hel)+np.log10(hosc))**2)
	#def residual(hel,hosc): return np.mean((np.log10(hel))**2)
	#def residual(hel,hosc): return np.mean((np.log10(hel/hosc))**2) ############### WHY WORKS?
	#see messed up vxxx def residual(hel,hosc): return np.mean((np.log10(hel*hosc-1.0))**2)
	#def residual(hel,hosc): return np.mean((np.log10(hel*(10**(-1*np.log10(hosc)))))**2)
	def residual(hel,hosc): return np.mean((np.log10(hel)-np.log10(hosc))**2)
	return residual

# shit looks the same to me really: ... x = 1./10;(10.*x);10/(10.**(-1.*np.log10(x)))

def prepare_objective(
	hqs,curvature_fields,wavevectors,area,
	curvature_sum_function,fft_function,band,residual_function,
	**kwargs):
	"""
	Package the objective function for optimization.
	"""
	#---the ndrops_uniform switch tells us how to multiplex a single field to multiple proteins
	ndrops_uniform = kwargs.pop('ndrops_uniform',0)
	#---default curvature sum method is the mean
	curvature_sum_method = kwargs.pop('curvature_sum_method','mean')
	#---binning methods are blurry (i.e. nonzero bins), perfect (i.e. one bin per unique magnitude)
	#---...or explicit, which does no binning whatsoever
	binner_method = kwargs.pop('binner_method','blurry')
	perfect_collapser = kwargs.pop('perfect_collapser',None)
	blurry_binner = kwargs.pop('blurry_binner',None)
	signterm = kwargs.pop('inner_sign',1.0)
	positive_vibe = kwargs.pop('positive_vibe',True)
	imaginary_mode = kwargs.pop('imaginary_mode','complex')
	fix_curvature = kwargs.pop('fix_curvature',None)
	oscillator_function_local = kwargs.pop('oscillator_function',None)
	if oscillator_function_local==None:
		oscillator_function_local = prepare_oscillator_function(reverse=True)
	residual = kwargs.pop('residual_function',prepare_residual())
	weighting_scheme = kwargs.pop('weighting_scheme',None)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---consistency checks
	if binner_method=='blurry' and type(blurry_binner)==type(None): 
		raise Exception('send blurry_binner')
	elif binner_method=='perfect' and type(perfect_collapser)==type(None): 
		raise Exception('send perfect_collapser')
	#---rename for the function
	cfs,q_raw = curvature_fields,wavevectors
	#---compute weights for specific weighting schemes	
	if weighting_scheme=='blurry_explicit':
		q_raw_binned_temp,_,_,q_mapping = blurry_binner(q_raw,q_raw,return_mapping=True)
		weights = 1./np.array([len(i) for i in q_mapping])[band]
	else: weights = None

	def multipliers(x,y): 
		"""Multiplying complex matrices in the list of terms that contribute to the energy."""
		return x*np.conjugate(y)

	#---fitting kappa, gamma, and the vibration with a fixed curvature
	#---note that this use-case is meant for matching to the legacy code
	if type(fix_curvature)!=type(None):
		if ndrops_uniform==0:
			raise Exception('you can only fix curvature for the uniform case right now')
		curvatures = [fix_curvature for i in range(ndrops_uniform)]
		#---generate a composite field from 
		composite = curvature_sum_function(cfs,curvatures,method=curvature_sum_method)
		#---Fourier transform		
		cqs = fft_function(composite)
		#---construct the terms in eqn 23
		termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
		#---reshape the terms into one-dimensional lists, dropping the zeroth mode
		termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
		#---! explain logic behind real/imaginary here
		if imaginary_mode=='real': termlist = [np.real(k) for k in termlist]
		elif imaginary_mode=='complex': termlist = [np.abs(k) for k in termlist]
		else: raise Exception('invalid imaginary_mode %s'%imaginary_mode)
		def objective(args,mode='residual',return_spectrum=False):
			"""
			Fit parameters are defined in sequence for the optimizer.
			They are: kappa,gamma,vibe,*curvatures-per-dimple.
			"""
			if len(args)>3: raise Exception('only send three arguments')
			#---the first three arguments are always bending rigitidy, surface tension, vibration correction
			kappa,gamma,vibe = args[:3]
			#---CONSTRAINTS are enforced here
			if positive_vibe: vibe = np.abs(vibe)
			#---constructing the elastic Hamiltonian based on the wavevectors
			hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
				+signterm*termlist[2]*q_raw**2+termlist[3])
				+gamma*area*(termlist[0]*q_raw**2))
			hosc = oscillator_function_local(vibe,q_raw)
			#---apply the vibration correction
			ratio = hel*hosc
			#---note that the band is prepared in advance above
			if binner_method=='explicit': ratio = hel
			elif binner_method=='perfect':
				q_binned,ratio,q_binned_inds = perfect_collapser(q_raw,hel)
			elif binner_method=='blurry':
				q_binned,ratio,q_binned_inds = blurry_binner(q_raw,hel)
			else: raise Exception('invalid binner_method %s'%binner_method)
			#----compute residuals with relevant wavevectors (in the band) and return
			if mode=='residual': 
				if type(weights)!=type(None): value = residual(weights*hel[band],weights*hosc[band])
				else: value = residual(hel[band],hosc[band])
			#elif mode=='ratio': value = ratio
			elif mode=='elastic': value = hel
			else: raise Exception('invalid residual mode %s'%mode)
			return value
		return objective

	def objective(args,mode='residual',return_spectrum=False):
		"""
		Fit parameters are defined in sequence for the optimizer.
		They are: kappa,gamma,vibe,*curvatures-per-dimple.
		"""
		#---the first three arguments are always bending rigitidy, surface tension, vibration correction
		kappa,gamma,vibe = args[:3]
		#---CONSTRAINTS are enforced here
		if positive_vibe: vibe = np.abs(vibe)
		#---uniform curvatures are multiplexed here
		if ndrops_uniform!=0: curvatures = [args[3] for i in range(ndrops_uniform)]
		#---one curvature per field
		else: curvatures = args[3:]
		#---generate a composite field from the individual fields
		composite = curvature_sum_function(cfs,curvatures,method=curvature_sum_method)
		#---Fourier transform		
		cqs = fft_function(composite)
		#---construct the terms in eqn 23
		termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
		#---reshape the terms into one-dimensional lists, dropping the zeroth mode
		termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
		#---! explain logic behind real/imaginary here
		if imaginary_mode=='real': termlist = [np.real(k) for k in termlist]
		elif imaginary_mode=='complex': termlist = [np.abs(k) for k in termlist]
		else: raise Exception('invalid imaginary_mode %s'%imaginary_mode)
		#---constructing the elastic Hamiltonian based on the wavevectors
		hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
			+signterm*termlist[2]*q_raw**2+termlist[3])
			+gamma*area*(termlist[0]*q_raw**2))
		#---apply the vibration correction
		hosc = oscillator_function_local(vibe,q_raw)
		#ratio = hel*hosc
		#---note that the band is prepared in advance above
		if binner_method=='explicit': ratio = hel
		elif binner_method=='perfect':
			q_binned,ratio,q_binned_inds = perfect_collapser(q_raw,hel)
		elif binner_method=='blurry':
			q_binned,ratio,q_binned_inds = blurry_binner(q_raw,hel)
		else: raise Exception('invalid binner_method %s'%binner_method)
		#----compute residuals with relevant wavevectors (in the band) and return
		if mode=='residual':
			if type(weights)!=type(None): value = residual(weights*ratio[band],weights*hosc[band])
			else: value = residual(hel[band],hosc[band])
		#elif mode=='ratio': value = ratio
		elif mode=='elastic': value = ratio
		else: raise Exception('invalid residual mode %s'%mode)
		return value

	#---return the decorated function
	return objective

###---CONTROLLER

class InvestigateCurvature:

	"""
	Manage the entire curvature coupling calculation.
	"""

	def __init__(self,**kwargs):
		"""Assemble the necessary modules and execute."""
		#---sometimes we want to instantiate this class without running the calculation
		self.do_calculation = kwargs.pop('do_calculation',True)
		#---unpack large incoming data
		self.work = kwargs.pop('work',None)
		#---! currently set for the optimizer to run one simulation but retaining the loops for later
		self.sns = [kwargs.pop('sn')]
		if not self.work: raise Exception('we require an instance of the omnicalc workspace')
		self.data = kwargs.pop('undulations',None)
		if not self.data: raise Exception('send upstream undulations data')
		self.design = kwargs.pop('design',{})
		self.fitting_parameters = kwargs.pop('fitting',{})
		self.style = self.design.get('style',None)
		self.single_mode = kwargs.pop('single_mode',True)
		if not self.style: raise Exception('invalid style %s'%self.style)
		self.data_prot_incoming = kwargs.pop('protein_abstractor',None)
		if not self.data_prot_incoming: raise Exception('send upstream protein_abstractor data')
		#---allow the user to send through previous data from the memory to avoid repeating previous steps
		remember = kwargs.pop('remember',{})
		if kwargs: raise Exception('unprocessed arguments %s'%kwargs)
		#---reformat some data if only one simulation. legacy versions of this code analyzed multiple 
		#---...simulations at once and we retain some of this design for later
		if self.single_mode:
			if not len(self.sns)==1: raise Exception('you can only have one simulation name in single_mode')
			self.data_prot_incoming = {self.sns[0]:{'data':self.data_prot_incoming}}
			self.data = {self.sns[0]:{'data':self.data}}
		#---the following loaders allow you to manipulate incoming data. defaults are provided however we 
		#---...strongly recommend that you retain the specification in the yaml file in case you change it,
		#---...so that omnicalc can distinguish different trials. obviously omnicalc cannot track changes
		#---...to the loader codes, so you should soft-code those parameters via the yaml file
		#---get the "protein" positions
		self.loader_spec_protein = self.design.get('loader_protein',
			{'module':'codes.curvature_coupling_loader','function':'curvature_coupling_loader_protein'})
		self.loader_func_protein = self.gopher(
			self.loader_spec_protein,module_name='module',variable_name='function')
		self.data_prot = self.loader_func_protein(data=dict(protein_abstractor=self.data_prot_incoming))
		#---get the "membrane" positions and send them through the FFT in the loader function
		self.loader_spec = self.design.get('loader_membrane',
			{'module':'codes.curvature_coupling_loader','function':'curvature_coupling_loader_membrane'})
		self.loader_func = self.gopher(self.loader_spec,module_name='module',variable_name='function')
		#---midplane method is "flat" by default
		midplane_method = self.design.get('midplane_method','flat')
		self.memory = self.loader_func(data=dict(undulations=self.data),midplane_method=midplane_method)
		#---load remembered data into memory. the curvature function may notice this and avoid repetition
		self.memory.update(**remember)
		#---the style should be a function in this class
		if self.do_calculation: self.finding = getattr(self,self.style)()

	def gopher(self,spec,module_name,variable_name):
		"""Load an external module. Useful for changing the workflow without changing the code."""
		return gopher(spec,module_name,variable_name,work=self.work)

	###---OPTIMIZED SOLUTIONS

	def drop_gaussians(self,**kwargs):
		"""
		Method for choosing the positions of Gaussians.
		"""
		pos_spec = kwargs.get('curvature_positions',{})
		method = pos_spec.get('method',None)
		extent = kwargs.get('extents',{}).get('extent',{})
		if method=='protein_subselection':
			for sn in self.sns:
				selections = pos_spec.get('selections',None)
				if not selections: raise Exception('need selections in protein_subselection')
				#---determine the centers of the protein according to the selections
				#---...noting that the protein_abstractor points are stored by the residue, not bead/atom 
				points = np.array([np.transpose(self.data_prot[sn]['data']['points'],(1,0,2))[s] 
					for s in selections])
				points = points.mean(axis=1)[...,:2]
				#---save the points for later
				self.memory[(sn,'drop_gaussians_points')] = points
				ndrops = len(points)
				#---get data from the memory
				hqs = self.memory[(sn,'hqs')]
				self.nframes = len(hqs)
				mn = hqs.shape[1:]
				vecs = self.memory[(sn,'vecs')]
				vecs_mean = np.mean(vecs,axis=0)
				#---formulate the curvature request
				curvature_request = dict(curvature=1.0,mn=mn,sigma_a=extent,sigma_b=extent,theta=0.0)
				#---construct unity fields
				fields_unity = np.zeros((self.nframes,ndrops,mn[0],mn[1]))
				reindex,looper = zip(*[((fr,ndrop),
					dict(vecs=vecs[fr],centers=[points[ndrop][fr]/vecs[fr][:2]],**curvature_request)) 
					for fr in range(self.nframes) for ndrop in range(ndrops)])
				status('computing curvature fields for %s'%sn,tag='compute')
				incoming = basic_compute_loop(make_fields,looper=looper)
				#---! inelegant
				for ii,(fr,ndrop) in enumerate(reindex): fields_unity[fr][ndrop] = incoming[ii]
				self.memory[(sn,'fields_unity')] = fields_unity
		elif method=='pixel':
			#---recall that the loop over sns is pretty much redundant
			for sn in self.sns:
				#---construct a box-vector-scaled grid of points which we call "pixels"
				#---get data from the memory
				hqs = self.memory[(sn,'hqs')]
				self.nframes = len(hqs)
				mn = hqs.shape[1:]
				vecs = self.memory[(sn,'vecs')]
				vecs_mean = np.mean(vecs,axis=0)
				#---get the grid spacing from the metadata
				spacer = pos_spec.get('spacer',None)
				spacer_x,spacer_y = [pos_spec.get('spacer_%s'%i,spacer) for i in 'xy']
				npts = (vecs_mean[:2]/np.array([spacer_x,spacer_y])).astype(int)
				posts = np.array([[np.linspace(0,vecs[fr][d],npts[d]+1) 
					for d in range(2)] for fr in range(self.nframes)])
				fence = np.array([[(posts[fr][d][1:]+posts[fr][d][:-1])/2. for d in range(2)]
					for fr in range(self.nframes)])
				points = np.array([np.concatenate(np.transpose(np.meshgrid(*fence[fr]))) 
					for fr in range(self.nframes)])
				ndrops = len(points[0])
				#---formulate the curvature request
				curvature_request = dict(curvature=1.0,mn=mn,sigma_a=extent,sigma_b=extent,theta=0.0)
				#---construct unity fields
				fields_unity = np.zeros((self.nframes,ndrops,mn[0],mn[1]))
				reindex,looper = zip(*[((fr,ndrop),
					dict(vecs=vecs[fr],centers=[points[fr][ndrop]/vecs[fr][:2]],**curvature_request)) 
					for fr in range(self.nframes) for ndrop in range(ndrops)])
				status('computing curvature fields for %s'%sn,tag='compute')
				incoming = basic_compute_loop(make_fields,looper=looper)
				#---! inelegant
				for ii,(fr,ndrop) in enumerate(reindex): fields_unity[fr][ndrop] = incoming[ii]
				self.memory[(sn,'fields_unity')] = fields_unity
				self.memory[(sn,'drop_gaussians_points')] = points
		elif method=='neighborhood':
			#---extra distance defines a border around the average hull
			extra_distance = pos_spec['distance_cutoff']
			spacer = pos_spec['spacer']
			def rotate2d(pts,angle):
				x = pts[:,0]*np.cos(angle)-pts[:,1]*np.sin(angle)
				y = pts[:,1]*np.cos(angle)+pts[:,0]*np.sin(angle)
				return np.transpose((x,y))
			def arange_symmetric(a,b,c):
				return np.unique(np.concatenate((np.arange(a,b,c),-1*np.arange(a,b,c))))
			for sn in self.sns:
				###---!!! beware this code might have an indexing problem !!!
				#---! nope now that you fixed the index error each protein gets its own neighborhood 
				#---! ...presumably with an indeterminate position
				#---for each frame we compute the centroid and orientation
				points_all = self.data_prot[sn]['data']['points_all']
				cogs = points_all.mean(axis=2).mean(axis=1)[:,:2]
				#---get the average set of points
				average_pts = points_all.mean(axis=0).mean(axis=0)[:,:2]
				average_pts -= average_pts.mean(axis=0)
				average_axis = principal_axis(average_pts)
				angle = np.arccos(np.dot(vecnorm(average_axis),[1.0,0.0]))
				direction = 1.0-2.0*(np.cross(vecnorm(average_axis),[1.0,0.0])<0)
				rot = rotate2d(average_pts,direction*angle)
				#---get the span of the points plus the extra distance in each direction
				span_x,span_y = np.abs(rot).max(axis=0) + extra_distance
				ref_grid = np.concatenate(np.transpose(np.meshgrid(
					arange_symmetric(0,span_x,spacer),arange_symmetric(0,span_y,spacer))))
				#import matplotlib as mpl;import matplotlib.pyplot as plt;ax = plt.subplot(111);
				#ax.scatter(*average_pts.T);plt.show()
				#import ipdb;ipdb.set_trace()
				vecs = self.memory[(sn,'vecs')]
				vecs_mean = np.mean(vecs,axis=0)
				#---for each frame, map the ref_grid onto the principal axis
				self.nframes = len(points_all)
				points = np.zeros((len(ref_grid),self.nframes,2))
				for fr in range(self.nframes):
					pts = points_all[fr].mean(axis=0)[:,:2]
					offset = pts.mean(axis=0)
					average_axis = principal_axis(pts-offset)
					angle = np.arccos(np.dot(vecnorm(average_axis),[1.0,0.0]))
					direction = 1.0-2.0*(np.cross(vecnorm(average_axis),[1.0,0.0])<0)
					ref_grid_rot = rotate2d(ref_grid,direction*angle)+offset
					#---handle PBCs by putting everything back in the box
					ref_grid_rot_in_box = (ref_grid_rot + 
						(ref_grid_rot<0)*vecs[fr,:2] - (ref_grid_rot>=vecs[fr,:2])*vecs[fr,:2])
					points[:,fr] = ref_grid_rot_in_box
				#---debug with a plot if desired
				if False:
					import matplotlib.pyplot as plt
					plt.scatter(*ref_grid_rot_in_box.T)
					from base.store import picturesave
					fn = 'fig.DEBUG.curvature_undulation_coupling_neighborhood'
					picturesave(fn,self.work.plotdir)
					raise Exception('dropped debugging image for your review and deletion '
						'to %s. remove it and then turn this debugger off to continue'%fn)
				#---save the position of the curvature fields for later
				self.memory[(sn,'drop_gaussians_points')] = points
				ndrops = len(ref_grid)
				#---! ULTRA REPETITIVE WITH THE OTHER OPTIONS
				#---get data from the memory
				hqs = self.memory[(sn,'hqs')]
				mn = hqs.shape[1:]
				#---formulate the curvature request
				curvature_request = dict(curvature=1.0,mn=mn,sigma_a=extent,sigma_b=extent,theta=0.0)
				#---construct unity fields
				fields_unity = np.zeros((self.nframes,ndrops,mn[0],mn[1]))
				reindex,looper = zip(*[((fr,ndrop),
					dict(vecs=vecs[fr],centers=[points[ndrop][fr]/vecs[fr][:2]],**curvature_request)) 
					for fr in range(self.nframes) for ndrop in range(ndrops)])
				status('computing curvature fields for %s'%sn,tag='compute')
				incoming = basic_compute_loop(make_fields,looper=looper,run_parallel=True)
				#---! inelegant
				for ii,(fr,ndrop) in enumerate(reindex): fields_unity[fr][ndrop] = incoming[ii]
				self.memory[(sn,'fields_unity')] = fields_unity
				if False:
					import matplotlib as mpl;import matplotlib.pyplot as plt;
					plt.imshow(fields_unity[0][0].T);plt.show()
					import ipdb;ipdb.set_trace()
		#---one field per protein, for all proteins
		elif method in ['protein_dynamic_single','protein_dynamic_single_uniform',
			#---the catchall represents curvatures saved for both methods for drilldown later
			'protein_dynamic_single_catchall']:
			#---precomputed values or subsequent reanalysis precludes the construction of the fields
			if (set([(self.sns[0],k) for k in ['sampling','fields_unity','drop_gaussians_points']])
				<=set(self.memory.keys())): 
				self.nframes = len(self.memory[(self.sns[0],'fields_unity')])
				self.sampling = self.memory[(self.sns[0],'sampling')]
				return
			#---! the following code is very repetitive with the protein subselection method
			for sn in self.sns:
				#---points_all is nframes by proteins by beads/atoms by XYZ
				points = self.data_prot[sn]['data']['points_all'].mean(axis=2)[...,:2].transpose(1,0,2)
				#---save the points for later
				self.memory[(sn,'drop_gaussians_points')] = points
				ndrops = len(points)
				#---get data from the memory
				hqs = self.memory[(sn,'hqs')]
				if 'nframes' not in pos_spec:
					sampling = np.arange(0,len(hqs),pos_spec.get('frequency',1))
				#---nframes gets exactly the right number of frames
				else: sampling = np.linspace(0,len(hqs)-1,pos_spec['nframes']).astype(int)
				mn = hqs.shape[1:]
				vecs = self.memory[(sn,'vecs')]
				vecs_mean = np.mean(vecs,axis=0)
				#---formulate the curvature request
				curvature_request = dict(curvature=1.0,mn=mn,sigma_a=extent,sigma_b=extent,theta=0.0)
				#---construct unity fields
				fields_unity = np.zeros((len(sampling),ndrops,mn[0],mn[1]))
				reindex,looper = zip(*[((fr,ndrop),
					dict(vecs=vecs[fr],centers=[points[ndrop][fr]/vecs[fr][:2]],**curvature_request)) 
					for fr in sampling for ndrop in range(ndrops)])
				status('computing curvature fields for %s'%sn,tag='compute')
				incoming = basic_compute_loop(make_fields,looper=looper)
				#---! inelegant
				sampling_reindex = dict([(key,ii) for ii,key in enumerate(sampling)])
				for ii,(fr,ndrop) in enumerate(reindex): 
					fields_unity[sampling_reindex[fr]][ndrop] = incoming[ii]
				self.memory[(sn,'sampling')] = self.sampling = sampling
				self.memory[(sn,'fields_unity')] = fields_unity
				self.nframes = len(sampling)
		else: raise Exception('invalid selection method')

	def curvature_sum(self,cfs,curvatures,**kwargs):
		"""
		Curvature fields are maxed not summed.
		"""
		return curvature_sum_function(cfs,curvatures,**kwargs)

	def wilderness(self,**kwargs):
		"""
		Wandering on the error landscape to optimize the curvature coupling hypothesis.
		"""
		bundle,solutions = {},{}
		global Nfeval
		spec = self.design
		extents = spec.get('extents',{})
		extents_method = extents.get('method')
		curvature_sum_method = self.design['curvature_sum']
		#---special handling for pixel extents if None: extent is half the spacer
		#---! should this check of the curvature position method is pixel?
		if extents.get('extent',False)==None and spec['curvature_positions']['method']:
			extents['extent'] = spec['curvature_positions']['spacer']/2.
		#---flag for uniform or variable curvatures
		do_uniform_curvature = spec['curvature_positions']['method']=='protein_dynamic_single_uniform'
		#---prepare curvature fields
		if extents_method=='fixed_isotropic': self.drop_gaussians(**spec)
		elif extents_method=='protein_dynamic': pass
		else: raise Exception('invalid extents_method %s'%extents_method)
		#---in the previous version we manually saved the data (and checked if it was already computed here)
		#---...however this class has been ported directly into omnicalc now
		#---optimize over simulations, however only one simulation will ever be present in the current code
		for snum,sn in enumerate(self.sns):
			status('starting optimization for %s %d/%d'%(sn,snum+1,len(self.sns)),tag='optimize')
			#---! we could apply a second filter here for the fit, after the filter for the frames
			#---! ... frameslice = np.arange(0,self.nframes,spec.get(
			#---! ...     'curvature_positions',{}).get('frequency_fit',1))
			frameslice = self.sampling

			#---load the source data
			hqs = self.memory[(sn,'hqs')][frameslice]
			#---previously saved a composite curvature field but this condition was removed in favor of
			#---...saving the separate curvature fields, even for re-optimization
			needs_curvature_sum = True
			if (sn,'fields_unity') in self.memory: 
				cfs = self.memory[(sn,'fields_unity')]
			vecs = self.memory[(sn,'vecs')][frameslice]
			ndrops = cfs.shape[1]
			ndrops_uniform = 0 
			if do_uniform_curvature: 
				ndrops_uniform = int(ndrops)
				ndrops = 1

			#---formulate the wavevectors
			Lx,Ly = np.mean(vecs,axis=0)[:2]
			q_raw = formulate_wavevectors(lx=Lx,ly=Ly,dims=np.shape(hqs)[1:])['wavevectors']
			area = Lx * Ly

			tweak = self.fitting_parameters
			#---definition of inner sign is hard-coded here. used in one other place (below) which should
			#---...have the same value. this term is the same sign as height. the optimizer will find the 
			#---...right sign, so this parameter really only matters for the purposes of interpreting 
			#---...the curvature
			signterm = tweak.get('inner_sign',1.0)
			initial_kappa = tweak.get('initial_kappa',25.0)
			lowcut = kwargs.get('lowcut',tweak.get('low_cutoff',0.0))
			band = cctools.filter_wavevectors(q_raw,low=lowcut,high=tweak.get('high_cutoff',2.0))
			residual_form = kwargs.get('residual_form',tweak.get('residual_form','log'))
			binner_method = spec.get('binner','explicit')
			weighting_scheme = spec.get('weighting_scheme','standard')
			if binner_method=='explicit': q_raw_binned = q_raw
			elif binner_method=='perfect':
				q_raw_binned,_,_ = perfect_collapser(q_raw,q_raw)
			elif binner_method=='blurry':
				q_raw_binned,_,_ = blurry_binner(q_raw,q_raw)
			else: raise Exception('invalid binner_method %s'%binner_method)
			band = cctools.filter_wavevectors(q_raw_binned,low=lowcut,high=tweak.get('high_cutoff',2.0))
			#---weighting scheme
			if False:
				if weighting_scheme=='blurry_explicit' and binner_method!='explicit': 
					raise Exception('incompatible')
				if weighting_scheme=='blurry_explicit':
					q_raw_binned_temp,_,_,q_mapping = blurry_binner(q_raw,q_raw,return_mapping=True)
					weights = 1./np.array([len(i) for i in q_mapping])[band]

			if False:

				#---! remove this
				#---residual mode for the return format from the objective
				if residual_form == 'log':
					if weighting_scheme=='blurry_explicit':
						def residual(values):
							return np.sum(weights*np.log10(values.clip(min=machine_eps))**2)/float(len(values))
					else: 
						def residual(values):
							return np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
				#---deprecated residual method
				elif residual_form == 'linear':
					raise Exception('linear is not the right residual option') 
					def residual(values): 
						return np.sum((values-1.0)**2)/float(len(values))
				else: raise Exception('unclear residual form %s'%residual_form)

				def multipliers(x,y): 
					"""Multiplying complex matrices in the list of terms that contribute to the energy."""
					return x*np.conjugate(y)

				def callback(args):
					"""Watch the optimization."""
					global Nfeval
					name_groups = ['kappa','gamma','vibe']+['curve(%d)'%i for i in range(ndrops)]
					text = ' step = %d '%Nfeval+' '.join([name+' = '+dotplace(val)
						for name,val in zip(name_groups,args)+[('error',objective(args))]])
					status('searching! '+text,tag='optimize')
					Nfeval += 1

				#---! phasing this out so currently renamed
				#---! ...soon we will rerun with the at-large objective function from a decorator
				def objective_interior(args,mode='residual'):
					"""
					Fit parameters are defined in sequence for the optimizer.
					They are: kappa,gamma,vibe,*curvatures-per-dimple.
					"""
					kappa,gamma,vibe = args[:3]
					vibe = np.abs(vibe)
					#---uniform curvatures
					if ndrops_uniform!=0: curvatures = [args[3] for i in range(ndrops_uniform)]
					#---one curvature per field
					else: curvatures = args[3:]
					composite = self.curvature_sum(cfs,curvatures,method=curvature_sum_method)
					cqs = cctools.fft_field(composite)
					termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
					termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
					#---skipping assertion and dropping imaginary
					#---! explain logic behind real/imaginary here
					termlist = [np.real(k) for k in termlist]
					hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
						+signterm*termlist[2]*q_raw**2+termlist[3])
						+gamma*area*(termlist[0]*q_raw**2))
					#---! incorrect: ratio = hel/((vibe*q_raw+machine_eps)/(np.exp(vibe*q_raw)-1)+machine_eps)
					ratio = hel/((vibe*q_raw+machine_eps)/((np.exp(vibe*q_raw)-1)+machine_eps))
					#---note that the band is prepared in advance above
					if binner_method=='explicit': ratio = ratio
					elif binner_method=='perfect':
						q_binned,ratio,q_binned_inds = perfect_collapser(q_raw,ratio)
					elif binner_method=='blurry':
						q_binned,ratio,q_binned_inds = blurry_binner(q_raw,ratio)
					else: raise Exception('invalid binner_method %s'%binner_method)
					#----compute residuals and return
					if mode=='residual': return residual(ratio[band])
					elif mode=='ratio': return ratio
					else: raise Exception('invalid residual mode %s'%mode)


			def residual(hel,hosc): return np.mean((np.log10(hel)-np.log10(hosc))**2)

			residual_function = residual

			objective = prepare_objective(
				hqs=hqs,curvature_fields=cfs,
				wavevectors=q_raw,area=area,
				curvature_sum_function=curvature_sum_function,fft_function=cctools.fft_field,
				band=band,residual_function=residual_function,blurry_binner=blurry_binner,
				binner_method='explicit',positive_vibe=True,inner_sign=1.0,
				ndrops_uniform=ndrops_uniform)

			def callback(args):
				"""Watch the optimization."""
				global Nfeval
				name_groups = ['kappa','gamma','vibe']+['curve(%d)'%i for i in range(ndrops)]
				text = ' step = %d '%Nfeval+' '.join([name+' = '+dotplace(val)
					for name,val in zip(name_groups,args) ])#[('error',objective(args))]])
				status('searching! '+text,tag='optimize')
				Nfeval += 1

			Nfeval = 0
			optimize_method = self.design.get('optimize_method','Nelder-Mead')
			if optimize_method!='wait':
				#---either one distinct curvature per field or one curvature for all fields
				initial_conditions = [initial_kappa,0.0,0.001]+[0.0 for i in range(ndrops)]
				test_ans = objective(initial_conditions)
				if not isinstance(test_ans,np.floating): 
					raise Exception('objective_residual function must return a scalar')
				fit = scipy.optimize.minimize(objective,
					x0=tuple(initial_conditions),method=optimize_method,
					callback=callback)
				#---package the result
				bundle[sn] = dict(fit,success=str(fit.success))
				solutions['x'] = bundle[sn].pop('x')
				#---duplicate curvatures for uniform fields
				if ndrops_uniform!=0:
					solutions['x'] = np.concatenate((solutions['x'][:3],[solutions['x'][3] 
						for i in range(ndrops_uniform)]))
				#---cannot save final_simplex to attributes if it is available
				try: 
					final_simplex = bundle[sn].pop('final_simplex')
					solutions['final_simplex_0'] = final_simplex[0]
					solutions['final_simplex_1'] = final_simplex[1]
				except: pass
				#---not all integrators have the jacobean
				try: solutions['jac'] = bundle[sn].pop('jac')
				except: pass
			if needs_curvature_sum and optimize_method!='wait':
				solutions['cf'] = np.array(self.curvature_sum(cfs,fit.x[3:],
					method=curvature_sum_method).mean(axis=0))
				solutions['cf_first'] = np.array(self.curvature_sum(cfs,fit.x[3:],
					method=curvature_sum_method)[0])
				#---save explicit fields
				if spec.get('store_instantaneous_fields',False):
					solutions['cfs'] = np.array(self.curvature_sum(
						cfs,fit.x[3:],method=curvature_sum_method))
			if needs_curvature_sum:
				if spec.get('store_instantaneous_fields_explicit',False):
					solutions['fields_unity'] = self.memory[(sn,'fields_unity')]
			#---we also save the dropped gaussian points here
			if extents_method=='fixed_isotropic': 
				solutions['drop_gaussians_points'] = self.memory[(sn,'drop_gaussians_points')]
			else: raise Exception('invalid extents_method %s'%extents_method)
			if optimize_method!='wait': 
				solutions['ratios'] = objective(fit.x,mode='elastic')
			solutions['qs'] = q_raw
			solutions['qs_binned'] = q_raw_binned
			solutions['sampling'] = self.sampling
		#---we return contributions to result,attrs for the calculation
		#---! bundle is indexed by sn but solutions is not. this is an awkward mixture of bookkeeping
		return dict(result=solutions,attrs=dict(bundle=bundle,spec=spec))

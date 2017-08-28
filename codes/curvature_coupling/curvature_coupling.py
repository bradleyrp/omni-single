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
		self.memory = self.loader_func(data=dict(undulations=self.data))
		#---the style should be a function in this class
		if self.do_calculation: self.finding = getattr(self,self.style)()

	def gopher(self,spec,module_name,variable_name):
		"""Load an external module. Useful for changing the workflow without changing the code."""
		#---note that this function was folded into the omnicalc codebase in omni/base/tools.py
		mod = importlib.import_module(spec[module_name])
		target = mod.__dict__.get(spec[variable_name],None)
		#---the database design might need work so we always export it
		mod.work = self.work
		if not target: raise Exception('add %s and %s to the specs'%(module_name,variable_name))
		return target

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

				#---nope .. now that you fixed the index error each protein gets its own neighborhood presumably with an indeterminate position

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
				#import matplotlib as mpl;import matplotlib.pyplot as plt;ax = plt.subplot(111);ax.scatter(*average_pts.T);plt.show()
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
				incoming = basic_compute_loop(make_fields,looper=looper)
				#---! inelegant
				for ii,(fr,ndrop) in enumerate(reindex): fields_unity[fr][ndrop] = incoming[ii]
				self.memory[(sn,'fields_unity')] = fields_unity
		else: raise Exception('invalid selection method')

	def curvature_sum(self,cfs,curvatures,**kwargs):
		"""
		Curvature fields are maxed not summed.
		"""
		method = kwargs.get('method')
		if method=='mean':
			combo = np.array([np.transpose(cfs,(1,0,2,3))[cc]*c 
				for cc,c in enumerate(curvatures)]).mean(axis=0)
		else: raise Exception('invalid method %s'%method)
		return combo

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
		#---prepare curvature fields
		if extents_method=='fixed_isotropic': self.drop_gaussians(**spec)
		else: raise Exception('invalid extents_method %s'%extents_method)
		#---in the previous version we manually saved the data (and checked if it was already computed here)
		#---...however this class has been ported directly into omnicalc now
		#---optimize over simulations
		for snum,sn in enumerate(self.sns):
			status('starting optimization for %s %d/%d'%(sn,snum+1,len(self.sns)),tag='optimize')

			#---load the source data
			hqs = self.memory[(sn,'hqs')][:self.nframes]
			cfs = self.memory[(sn,'fields_unity')][:self.nframes]
			vecs = self.memory[(sn,'vecs')][:self.nframes]
			ndrops = cfs.shape[1]

			#---formulate the wavevectors
			lenscale = 1.0
			m,n = mn = np.shape(hqs)[1:]
			Lx,Ly = np.mean(vecs,axis=0)[:2]
			q2d = lenscale*np.array([[np.sqrt(
				((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
				((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
				for j in range(0,n)] for i in range(0,m)])
			q_raw = np.reshape(q2d,-1)[1:]
			area = (Lx*Ly/lenscale**2)

			tweak = self.fitting_parameters
			signterm = tweak.get('inner_sign',-1.0)
			initial_kappa = tweak.get('initial_kappa',25.0)
			lowcut = kwargs.get('lowcut',tweak.get('low_cutoff',0.0))
			band = cctools.filter_wavevectors(q_raw,low=lowcut,high=tweak.get('high_cutoff',1.0))
			residual_form = kwargs.get('residual_form',tweak.get('residual_form','log'))
			if residual_form == 'log':
				def residual(values): 
					return np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
			elif residual_form == 'linear': 
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

			def objective(args,mode='residual'):
				"""
				Fit parameters are defined in sequence for the optimizer.
				They are: kappa,gamma,vibe,*curvatures-per-dimple.
				"""
				kappa,gamma,vibe = args[:3]
				curvatures = args[3:]
				composite = self.curvature_sum(cfs,curvatures,method=curvature_sum_method)
				cqs = cctools.fft_field(composite)
				termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
				termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
				#---skipping assertion and dropping imaginary
				termlist = [np.real(k) for k in termlist]
				hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma*area*(termlist[0]*q_raw**2))
				ratio = hel/((vibe*q_raw+machine_eps)/(np.exp(vibe*q_raw)-1)+machine_eps)
				if mode=='residual': return residual(ratio[band])
				elif mode=='ratio': return ratio
				else: raise Exception('invalid mode %s'%mode)

			Nfeval = 0
			initial_conditions = [initial_kappa,0.0,0.01]+[0.0 for i in range(ndrops)]
			test_ans = objective(initial_conditions)
			if not isinstance(test_ans,np.floating): 
				raise Exception('objective_residual function must return a scalar')
			fit = scipy.optimize.minimize(objective,
				x0=tuple(initial_conditions),method='SLSQP',callback=callback)
			#---package the result
			bundle[sn] = dict(fit,success=str(fit.success))
			solutions['x'] = bundle[sn].pop('x')
			solutions['jac'] = bundle[sn].pop('jac')
			solutions['cf'] = np.array(self.curvature_sum(cfs,fit.x[3:],
				method=curvature_sum_method).mean(axis=0))
			solutions['cf_first'] = np.array(self.curvature_sum(cfs,fit.x[3:],
				method=curvature_sum_method)[0])
			#---we also save the dropped gaussian points here
			if extents_method=='fixed_isotropic': 
				solutions['drop_gaussians_points'] = self.memory[(sn,'drop_gaussians_points')]
			else: raise Exception('invalid extents_method %s'%extents_method)
			solutions['ratios'] = objective(fit.x,mode='ratio')
			solutions['qs'] = q_raw
		#---we return contributions to result,attrs for the calculation
		return dict(result=solutions,attrs=dict(bundle=bundle,spec=spec))

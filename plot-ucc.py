#!/usr/bin/env python

"""
UCC
Undulation Curvature Coupling (v5, ca 2017.12.06)
"""

import time,copy,pprint
import scipy
import scipy.optimize
import scipy.interpolate

machine_eps = eps = np.finfo(float).eps

def basic_undulation_fitter(qs,hqs,area,**kwargs):
	"""
	Fit the undulations with no curvature.
	"""
	#! this code is legacy and needs checked carefully
	#!   it is used later on in the coupling function
	#!   one the ucc class works and matches undulate.py, this should be retired
	# directions
	fit_tension = kwargs.pop('fit_tension',True)
	fit_correction = kwargs.pop('fit_correction',True)
	fit_protrusion = kwargs.pop('fit_protrusion',True)
	print('[STATUS calling basic_undulation_fitter: %s'%(dict(
		fit_tension=fit_tension,fit_correction=fit_correction,
		fit_protrusion=fit_protrusion)))
	# settings
	q_min = kwargs.pop('q_min',0.0)
	q_cut = kwargs.pop('q_cut',4.0)
	binner_method = kwargs.pop('binner_method','explicit')
	bin_width = kwargs.pop('bin_width',0.05)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	# using complex here. see also build_hqhq 
	hqhq = (np.abs(hqs).reshape((len(hqs),-1)).mean(axis=0)[1:])**2
	band = np.all((qs>=q_min,qs<=q_cut),axis=0)
	def model(q_raw,sigma,kappa,gamma_p,vibe,
		protrusion=True,correct=True,tension=True,only_tension=False):
		kappa,gamma_p = np.abs(kappa),np.abs(gamma_p)
		if not tension: sigma = 0.0
		if protrusion==False:
			pure = ((1.0)/(area))*(
				(1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps))
		else: 
			pure = ((1.0)/(area))*(
				(1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps)+
				(1.)/(gamma_p*q_raw**2+machine_eps))
		# we include epsilon here to avoid divide by zero errors
		#   however note that this may affect the optimizer (see grease below)
		if correct: osc = (vibe*q_raw)*(1./(np.exp(vibe*q_raw)-1+machine_eps))
		else: osc = 1.0
		return pure * osc
	def residual(a,b): 
		# added epsilon here by default
		return ((np.log10(a/(b+machine_eps)))**2).mean()
	binner = Binner(wavevectors=qs[band],mode=binner_method,bin_width=bin_width)
	def objective(args,
		protrusion=fit_protrusion,
		correct=fit_correction,
		tension=fit_tension):
		(sigma,kappa,gamma_p,vibe) = args
		heights = hqhq[band]
		heights_model = model(qs[band],sigma,kappa,gamma_p,vibe,
			protrusion=protrusion,correct=correct,tension=tension)
		return residual(binner.bin(heights),binner.bin(heights_model))
	global stepno
	stepno = 1
	def callback(args,silent=False):
		"""Watch the optimization."""
		global stepno
		args_string = ' '.join([stringer(a) for a in args])
		output = ((u'\r' if stepno>0 else '')+
			'[OPTIMIZE] step %s: %s'%(stepno,args_string))
		if not silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		stepno += 1
	initial_conditions = (0.,20.,0.,0.)
	fit = scipy.optimize.minimize(objective,
		x0=initial_conditions,method='Nelder-Mead',callback=callback)
	fitted = model(qs,*fit.x)
	print('\n[OPTIMIZE] completed with error %s'%fit.fun)
	return dict(sigma=fit.x[0] if fit_tension else 0.0,band=band,hqhq=hqhq,
		kappa=np.abs(fit.x[1]),vibe=fit.x[3] if fit_correction else 0.0,
		gamma=fit.x[2] if fit_protrusion else 0.0,
		model=model)

@autoload(plotrun)
def loader():
	"""
	Load undulation and protein data into globals.
	"""
	global seepspace
	seepspace = 'plotspecs,data,calc'.split(',')
	if any([s not in globals() for s in seepspace]):
		# the following flexible names are meant to allow inputs from 
		#   other sources, namely mesoscale simulations
		plotspecs = work.plotspec.specs
		protein_abstractor_name = plotspecs.get('protein_abstractor_name','protein_abstractor')
		undulations_name = plotspecs.get('undulations_name','undulations')
		data,calc = plotload('ucc')
		# export to globals after loading
		for key in seepspace: globals()[key] = locals()[key]

# log utility functions
# logs use epsilon to prevent division by zero by default
def logdiff(x,y): return 10.**(np.log10(x+machine_eps)-np.log10(y+machine_eps))
def logsum(x,y): return 10.**(np.log10(x+machine_eps)+np.log10(y+machine_eps))

def formulate_wavevectors(dims,vecs=None,lx=None,ly=None,lenscale=1.0):
	"""
	Generate wavevectors from box dimensions and the number of grid points.
	"""
	if type(vecs)!=type(None):
		Lx,Ly = np.mean(vecs,axis=0)[:2]
	elif type(lx)==type(None) or type(ly)==type(None):
		raise Exception('send box dimensions (lx and ly) or framewise box vectors (vecs)')
	else: Lx,Ly = lx,ly
	if len(dims)!=2:
		raise Exception('dims must be the number of grid points in each of two directions')
	# formulate the wavevectors
	lenscale = 1.0
	m,n = mn = dims
	q2d = lenscale*np.array([[np.sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	q_raw = np.reshape(q2d,-1)[1:]
	area = (Lx*Ly/lenscale**2)
	return dict(wavevectors=q_raw,area=area,wavevectors_2d=q2d)

class Binner:
	"""
	Manage all questions of binning your wavevectors (prefer explicit for theory reasons).
	"""
	def __init__(self,mode,wavevectors,**kwargs):
		self.mode = mode
		self.bin_width = kwargs.pop('bin_width',0.05)
		self.blurry_focus = kwargs.pop('blurry_focus','average')
		if kwargs: raise Exception
		self.wavevectors = self.qs = wavevectors
		print('[STATUS] binner is mode: %s'%mode)
		# trivial case where all bins are unique
		if mode=='explicit': 
			self.indexer = np.transpose([np.arange(len(self.qs))])
			self.binned_independent = self.bin(self.qs)
		elif mode=='blurry':
			# collect abstract values
			self.abstract = (self.qs/self.bin_width).astype(int)
			self.uniques_abs,self.indexer,self.counts = np.unique(
				self.abstract,return_inverse=True,return_counts=True)
			# save the bins within each wavevector
			if self.blurry_focus=='average':
				# run the binner on the wavevectors to get the average wavevector for each bin
				# bins will not be perfectly spaced so you would add another method to make it even
				self.binned_independent = np.bincount(self.indexer,weights=self.qs)/self.counts
			else: raise Exception
		# the perfect binner averages anything with exactly the same magnitude
		elif mode=='perfect':
			self.uniques,self.indexer,self.counts = np.unique(
				self.qs,return_inverse=True,return_counts=True)
			self.binned_independent = self.bin(self.qs)
		else: raise Exception
		#! should we save weights for the blurry-explicit method
	def bin(self,values): 
		"""Fast binner method."""
		# note the speed and usability of bincount over manually vectorizing things
		if self.mode=='explicit': 
			return values[self.indexer].mean(axis=1)
		else: 
			try: out = np.bincount(self.indexer,weights=values)/self.counts
			except:
				print('[ERROR] lengths of indexer, values, counts: %d, %d, %d'%(
					len(self.indexer),len(values),len(self.counts)))
				import ipdb;ipdb.set_trace()
				raise
			return out

def fft_field(dat):
	"""
	Take the FFT of a regular grid.
	Assumes multiple frames.
	"""
	#---take the norm of the FFT
	factor = np.product(dat.shape[-2:])
	raw = np.fft.fft2(np.array(dat))
	#---hard coding the right method here!
	fft_flag = 'complex'
	if fft_flag == 'real': hqs = np.real(raw)
	elif fft_flag == 'absolute': hqs = np.absolute(raw)
	elif fft_flag in ['complex','complex_real']: hqs = raw
	else: raise
	return hqs/factor

def fftwrap(dat): 
	"""
	Simple wrapper for the 2D FFT algorithm in scipy.
	"""
	return np.fft.fft2(np.array(dat))

def multipliers(x,y): 
	"""Multiplying complex matrices in the list of terms that contribute to the energy."""
	return x*np.conjugate(y)

def make_isotropic_fields(foci,vecs,ngrid,magnitude=1.0,extent=1.0):
	"""
	Turn a set of positions, box vectors, and number of grid points into a set of isotropic fields.
	Twentyfold speedup over previous method due to maximal vectorizing.
	"""
	# incoming data should be spots by frames by beads by coordinate
	if len(foci.shape)==4: positions = foci
	elif len(foci.shape)==3: positions = np.tile(foci,(1,1,1,1))
	else: raise Exception
	# spots are XY positions averaged over the second axis which is the beads
	spots = positions[...,:2].mean(axis=2)
	# nprots is really the number of foci or spots
	nprots,nframes = spots.shape[:2]
	#! more elegant way to do this trick? itertools? np.unravel_index(1000,(nprots,nframes))?
	#! previously used the concatenate,transpose,meshgrid trick 
	#!   for indices, deprecated by the new trick
	inds = np.concatenate(np.transpose(np.meshgrid(np.arange(nprots),np.arange(nframes))))
	# compute the grid points 
	prop_pts = np.transpose(np.meshgrid(range(0,ngrid[0]),range(0,ngrid[1])))/ngrid.astype(float)
	# new durable trick for making a grid inside of a simulation
	xypts = (np.tile(np.reshape(prop_pts,(ngrid[0],ngrid[1],1,2)),(nframes,1))*vecs[...,:2])
	# switch to 1D indexing
	xypts = np.transpose(xypts,(2,0,1,3))
	#! check indexing
	xypts = xypts.reshape((nframes,-1,2))
	# compute all euclidean distances from each protein foci to grid points for all frames
	distances = np.array([np.linalg.norm(xypts-
		np.tile(spots[i].reshape((nframes,1,2)),(1,xypts.shape[1],1)),axis=2) 
		for i in range(len(spots))])
	# only a square, division, and exp stand between us and the field
	fields = np.exp(-distances**2/extent**2)
	# return to a square shape
	return fields.reshape((nprots,nframes,ngrid[0],ngrid[1]))

def residual_standard(hel,hosc): 
	# we use epsilon by default
	return np.mean((
		np.log10(hel+machine_eps)-
		np.log10(hosc+machine_eps))**2)

def residual_linear(hel,hosc): 
	# testing purposes only
	return np.mean((hel-hosc)**2)

def residual_modelfix(hel,hosc): 
	# testing purposes only
	return np.mean((
		np.log10(hel+machine_eps)+
		np.log10(hosc+machine_eps))**2)

def stringer(x,p=5):
	"""A nice way to print numbers."""
	return (' '*(1*(x>0)))+('%.1e'%x 
		if (abs(x)<10**(-1*p) or abs(x)>10**p) else ('{:>%df}'%(p+2+1*(x<0))).format(x))

class QuickOpt:
	def __init__(self,objective,init,**kwargs):
		self.optimize_method = kwargs.pop('optimize_method','Nelder-Mead')
		self.scipy_optimize_function = kwargs.pop('scipy_optimize_function','minimize')
		self.silent = kwargs.pop('silent',False)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		self.stepno = 0
		# get the objective function from globals
		self.objective = objective
		if self.scipy_optimize_function=='minimize':
			self.fit = scipy.optimize.minimize(self.objective,x0=tuple(init),
				callback=self.callback,method=self.optimize_method)
			self.fitted = self.fit.x
			self.error = self.fit.fun
		elif self.scipy_optimize_function=='fmin':
			self.fit = scipy.optimize.fmin(self.objective,x0=tuple(init),disp=True)
			self.error = self.objective(self.fit)
			self.fitted = self.fit
		else: raise Exception('unclear optimize procedure %s'%self.scipy_optimize_function)
		print('\n[RESULT] error: %s, solution: %s'%(self.error,str(self.fitted)))
	def callback(self,args):
		"""Watch the optimization."""
		args_string = ' '.join([stringer(a) for a in args])
		output = ((u'\r' if self.stepno>0 else '')+
			'[OPTIMIZE] step %s: %s'%(self.stepno,args_string))
		if not self.silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		self.stepno += 1

def log_reverse(raw):
	return 10.**(-1.0*np.log10(raw))

###
### Manage Fields
###
class UCCFields:
	"""
	Manage a set of fields.
	"""
	def __init__(self,parent=None,**kwargs):
		if not parent: raise Exception('dev')
		else: self.get_from_parent(parent)
	def get_from_parent(self,parent):
		"""Gather relevant information."""
		self.vecs = parent.vecs
		self.ngrid = parent.ngrid
		self.sn = parent.sn
		self.foci = parent.foci
		self.nframes = parent.nframes
		if not len(self.foci)==self.nframes: 
			raise Exception('foci must be frames by foci by beads by coordinate')
		self.nfoci = self.foci.shape[1]
	def curvature_multiplex(self,fields,*curvatures,**kwargs):
		"""Multiplex and sum the curvature fields."""
		if len(fields.shape)!=4: raise Exception
		if len(curvatures)==0: raise Exception('no curvatures')
		if len(curvatures)==1: curvatures = [curvatures[0] for i in fields]
		elif len(curvatures)!=len(fields): raise Exception('not enough curvatures nor only one')
		else: raise Exception('not enough curvatures')
		method = kwargs.pop('method','mean')
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		if method=='mean':
			return (ucc.curvature.fields.transpose(1,2,3,0)*curvatures).mean(axis=-1)
		else: raise Exception('invalid method %s'%method)
	def fields_from_hypo(self,**hypo):
		"""
		Generate many fields from a hypothesis. Useful for building curvature fields.
		"""
		extent = hypo.get('extent',1.0)
		if hypo['mode']=='protein_dynamic' and self.__dict__.get('extent_cursor',None)==extent:
			status('repeating with the same curvature field',tag='note')
		elif hypo['mode']=='protein_dynamic' and self.__dict__.get('extent_cursor',None)!=extent:
			self.fields = make_isotropic_fields(
				vecs=self.vecs,ngrid=self.ngrid,
				#! use foci_dims to figure out the right way to transpose?
				foci=self.foci.transpose(1,0,2,3),
				magnitude=1.0,extent=extent)
			# remember the extent to avoid recomputing things
			self.extent_cursor = extent
		else: raise Exception
	def fields_from_heights(self,heights):
		# get the heights and center them
		# self.fields = np.reshape(heights.reshape(self.nframes,-1).T - 
		#	heights.reshape(self.nframes,-1).mean(axis=1),
		#   (self.nframes,self.ngrid[0],self.ngrid[1]))
		self.fields = (heights.transpose(1,2,0)-
			np.array([np.mean(i) for i in heights])).transpose(2,0,1)
		# the height transform only happens once
		self.fields_q = fft_field(self.fields)

###
### Parent Class for Undulation Curvature Coupling
###
class UndulationCurvatureCoupling:
	"""
	Supervise the undulation-curvature coupling algorithm.
	"""
	model_style_valid = ['bending','tension','tension_protrusion','protrusion']
	def __init__(self,**kwargs):
		"""Load data."""
		self.sn = kwargs.pop('sn')
		self.grid_spacing = kwargs.pop('grid_spacing')
		# optional parameters
		self.couple_sign = kwargs.pop('couple_sign',1.0)
		self.model_style = kwargs.pop('model_style','bending')
		self.q_min = kwargs.pop('q_min',0.0)
		self.q_cut = kwargs.pop('q_cut',1.0)
		self.equipartition_model = kwargs.pop('equipartition_model','default')
		self.residual_name = kwargs.pop('residual_name','standard')
		self.positive_vibe = kwargs.pop('positive_vibe',True)
		self.oscillator_reverse = kwargs.pop('oscillator_reverse',False)
		# we use epsilon by default however previous notes suggest this
		#   should be investigated further because it may affect the optimizer
		self.grease = kwargs.pop('grease',True)
		self.binner_method = kwargs.pop('binner_method','explicit')
		self.bin_width = kwargs.pop('bin_width',0.05)
		self.subtract_protrusions = kwargs.pop('subtract_protrusions',False)
		#! settings for the mapping
		self.nstrengths = kwargs.pop('nstrengths',1)
		if kwargs: raise ValueError('kwargs remain %s'%kwargs)
		# valid model styles
		if self.model_style not in self.model_style_valid:
			raise Exception('model_style not in %s'%self.model_style_valid)
		# load heights from data
		global data
		if data==None: raise Exception
		self.build() #! future routing here

	def build(self):
		"""Prepare the calculation."""
		# vectors from the data
		#! get the specific names for later
		self.vecs = data['undulations'][self.sn]['data']['vecs']
		#! should there by a GMX switch that checks this?
		self.vecs_prot = data['protein_abstractor'][self.sn]['data']['vecs']
		if not np.all(self.vecs==self.vecs_prot):
			raise Exception('box vectors from focii and heights do not match')
		self.nframes = len(self.vecs)
		# save the "foci" or positions of supposed curvature fields
		# points_all from protein_abstractor has dimensions frames, protein, beads/atoms, xyz
		#! this will fail if the proteins have different 
		#!   numbers of beads. this case needs to be considered
		self.foci = data['protein_abstractor'][self.sn]['data']['points_all']
		self.foci_dims = ['frames','focus','beads','coordinate']
		self.ngrid = self.compute_grid_counts(self.vecs,self.grid_spacing)
		# prepare heights
		status('preparing heights',tag='compute')
		self.zs = UCCFields(self)
		self.heights = data['undulations'][self.sn]['data']['mesh'].mean(axis=0)
		self.zs.fields_from_heights(self.heights)
		# prepare wavevectors
		form_wavevectors = formulate_wavevectors(self.ngrid,vecs=self.vecs)
		self.qs_raw = form_wavevectors['wavevectors']
		self.qs_2d = form_wavevectors['wavevectors_2d']
		self.area = form_wavevectors['area']
		# get the residual function with a prefixed name from globals
		self.residual = globals()['residual_%s'%self.residual_name]
		# prepare the binner
		self.binner = Binner(
			wavevectors=self.qs_raw,
			mode=self.binner_method,
			bin_width=self.bin_width)
		self.qs = self.binner.binned_independent
		# prepare the band according to the reduced wavevectors from the binner
		self.band = np.where(np.all((self.qs>=self.q_min,self.qs<=self.q_cut),axis=0))[0]
		self.band_raw = np.where(np.all((
			self.qs_raw>=self.q_min,
			self.qs_raw<=self.q_cut),axis=0))[0]

		#!!??
		if 1:
			#! previously used this only when subtracting protrusions because later the
			#!   optimizer should set these values more carefully. later we might want to use
			#!   this to subtract the protrusions, but it would be better to fit everything together
			result = basic_undulation_fitter(self.qs_raw,self.zs.fields_q,
				self.area,q_cut=1.0,fit_tension=False,fit_correction=False,
				binner_method=self.binner_method,fit_protrusion=self.subtract_protrusions)
			# microscopic surface tension must be positive and is enforced in the fitter above
			# if we fit protrusions then we store the fits here for reference later
			self.gamma = np.abs(result['gamma'])
			self.kappa = np.abs(result['kappa'])
			self.vibe = np.abs(result['vibe'])
			self.hqhq = result['hqhq']

	def compute_grid_counts(self,vecs,grid_spacing):
		"""Compute the number of grid points to get as close as possible to a particular spacing."""
		return np.array([round(i) for i in np.mean(vecs,axis=0)/grid_spacing])[:2].astype(int)

	def build_curvature_fields(self,**hypo):
		status('preparing curvature fields',tag='compute')
		# make a fields object if absent
		if 'curvature' not in self.__dict__: self.curvature = UCCFields(self)
		# populate the fields object with curvature fields according to the hypothesis
		#! reduce reexecution for the same parameters
		self.curvature.fields_from_hypo(**hypo)

	def coupling(self,*curvatures):
		"""Compute the terms in equation 23."""
		hqs = self.zs.fields_q
		# protrusions are subtracted from the explicit data 
		#! discarded method: remove protrusions from the spectra before the fit. this method
		#!   was discarded for the obvious reason that it would be better to simply include
		#!   protrusions in the model. comments retained here for posterity
		#!     if self.subtract_protrusions: hqs = (hqs - 1./(self.qs_2d**2*self.gamma*self.area))
		# curvature sum function
		self.curvature.fields_sum = self.curvature.curvature_multiplex(
			self.curvature.fields,*curvatures)
		# we FFT after generating the final field
		cqs = self.curvature.fields_q = fft_field(self.curvature.fields_sum)
		termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
		#! discarded method
		if self.subtract_protrusions or False:
			self.termlist_alt = np.abs(np.reshape(np.mean(termlist[0],axis=0),-1)[1:])
			kappa = self.kappa
			#! why is kappa necessary here? it makes a difference for the answer
			m = 1./(self.area*kappa/2.*self.qs_2d**4+machine_eps)
			m2 = (
				1./(self.area*kappa/2.*self.qs_2d**4+machine_eps) + 
				1./(self.area*self.gamma*self.qs_2d**2+machine_eps))
			termlist[0] = logdiff(termlist[0],logdiff(m2,m))
			#! other attempts to properly filter things out
			#! termlist[0] = termlist[0]-1./(self.area*self.gamma*self.qs_2d**2)*self.qs_2d**4
			#! termlist[0] = logdiff(termlist[0]/self.qs_2d**4,
			#!   1./(self.area*self.gamma*self.qs_2d**2))*self.qs_2d**4
		# reshape the terms into one-dimensional lists, dropping the zeroth mode and converting to real
		self.termlist = np.array([
			np.abs(np.reshape(np.mean(k,axis=0),-1)[1:]) for k in termlist])
		# same exact result if you take the magnitude and average or average and then take the magnitude
		if 0: self.termlist_alt = np.array([
			(np.abs(k).reshape((len(k),-1)).mean(axis=0)[1:]) for k in termlist])
		# compare this to the previous function, the basic_undulation_fitter, which squared the mean
		#   of the hqs. this gives a slightly different result than self.termlist[0] which comes from
		#   np.abs(np.reshape(np.mean(multipliers(hqs,hqs),axis=0),-1)[1:]) which itself does the 
		#   matrix multiplication over complex numbers. we continue with the solution in termlist because
		#   both the wavevectors and FFT fields (hqs or cqs) should be complex
		# we hold self.hqs and self.hqhq for reference only
		self.hqs = self.zs.fields_q
		self.hqhq = (np.abs(hqs).reshape((len(hqs),-1)).mean(axis=0)[1:])**2	

		#!!?? 
		#! unclear what happens below but we have turned it off for now
		if self.subtract_protrusions:
			raise Exception('dev')
			model_q4 = 1./(self.area*self.kappa/2.*self.qs**4)
			# plot the standard model
			kappa,gamma_p,vibe = self.kappa,self.gamma,self.vibe
			sigma = 0.0
			q_raw = self.qs
			area = self.area
			model_basic = (
				((1.0)/(area))*((1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps)+
				(1.)/(gamma_p*q_raw**2+machine_eps))*
				((vibe*q_raw)*(1./(np.exp(vibe*q_raw)-1))))
			#! dev: enter debugging mode for alternate binning methods because size mismatch
			if self.binner_method!='explicit':
				import ipdb;ipdb.set_trace()
			self.termlist[0] = logsum(model_q4,logdiff(self.termlist[0],model_basic))

	def route_args_hel(self,*args):
		if self.model_style=='tension_protrusion': arglist = ['sigma','kappa','gamma']
		elif self.model_style=='tension': arglist = ['sigma','kappa']
		elif self.model_style=='bending': arglist = ['kappa']
		elif self.model_style=='protrusion': arglist = ['kappa','gamma']
		else: raise Exception
		if len(args)!=len(arglist)+(self.nstrengths if not self.fix_curvatures else 0): 
			raise Exception('incorrect arguments to the Hamiltonian for style %s: %s'%(
				self.model_style,args))
		out = dict(zip(arglist,args))
		if not self.fix_curvatures: out.update(curvatures=args[len(arglist):])
		return out

	def build_elastic_hamiltonian(self,*curvatures):
		"""Construct the H_{el} function, the elastic energy function."""
		# locals from self
		couple_sign = self.couple_sign
		# use qs_raw because binning happens downstream
		qs = self.qs_raw
		# note that we do not jitter the area here and we use NPT simulations so there
		#   should be slight variations in the area which are inconsequential
		area = self.area
		if len(curvatures)>0: 
			self.coupling(*curvatures)
			self.fix_curvatures = True
		else: self.fix_curvatures = False
		def hamiltonian(*args,debug=False):
			params = self.route_args_hel(*args)
			kappa = params.get('kappa')
			sigma = params.get('sigma',0.0)
			gamma = params.get('gamma',0.0)
			# if curvatures are free parameters we recompute them (slow)
			if not self.fix_curvatures: 
				curvatures = params.get('curvatures',())
				self.coupling(*curvatures)
			termlist = self.termlist
			if debug:
				import ipdb;ipdb.set_trace()
			return (area/2 * ( kappa * (
				termlist[0] * qs**4 
				+ couple_sign*termlist[1] * qs**2
				+ couple_sign*termlist[2] * qs**2 
				+ termlist[3])
				+ gamma * (termlist[0] * qs**2 )))
		# the elastic hamiltonian stays with the instance
		self.hel = hamiltonian

	def build_equipartition(self):
		"""Construct our working definition of equipartition."""
		if self.equipartition_model=='default':
			# if no oscillator correction we use equipartition here
			self.energy_per_mode = 1./2
			def equipartition_default(): 
				# we must use the binner qs to ensure we operate in explicit 
				#   wavevectors since binning will occur downstream
				return self.energy_per_mode * np.ones(self.binner.qs.shape[0])
			#! replace this function with free variables if oscillator
			self.equipartition = equipartition_default
		elif self.equipartition_model=='harmonic_oscillators':
			# for the oscillator correction we fit to 1 kBT 
			self.energy_per_mode = 1.
			self.build_oscillator_function()
			self.equipartition = self.oscillator
		else: raise Exception('invalid equipartition_model: %s'%equipartition_model)

	def build_oscillator_function(self):
		"""Build an harmonic oscillators function."""
		def oscillator_function(vibe):
			"""Model a series of harmonic oscillators."""
			# we use qs_raw because the binner will happen downstream
			# vibration must be nonnegative
			if self.positive_vibe: vibe = np.abs(vibe)
			# no grease here
			if not self.grease: 
				raise Exception('we have added epsilon to all opportunities for '
					'division by zero in this code so we recommend keeping the "grease" on '
					'in this location for consistency')
				raw = (vibe*self.qs_raw)*(0.+1./(np.exp(vibe*self.qs_raw)-1))
			# standard grease is outside the exponential
			else: 
				raw = (vibe*self.qs_raw)*(
					0.+1./(np.exp(vibe*self.qs_raw)-1+machine_eps))
			if not self.oscillator_reverse: return raw
			else: return log_reverse(raw)
		self.oscillator = oscillator_function

	def build_initial_conditions(self,**kwargs):
		"""Initial conditions must have the correct order. This must match route_args_objective."""
		init_kappa = kwargs.pop('kappa',20.)
		init_sigma = kwargs.pop('sigma',0.)
		init_gamma = kwargs.pop('gamma',0.)
		init_vibe = kwargs.pop('vibe',0.)
		init_curvature = kwargs.pop('curvature',0.)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		init = []
		if self.model_style=='tension_protrusion': init += [init_sigma,init_kappa,init_gamma]
		elif self.model_style=='bending': init += [init_kappa]
		elif self.model_style=='tension': init += [init_sigma,init_kappa]
		elif self.model_style=='protrusion': init += [init_kappa,init_gamma]
		else: raise Exception
		if self.equipartition_model=='default': init += []
		elif self.equipartition_model=='harmonic_oscillators': init += [init_vibe]
		else: raise Exception
		if self.fix_curvatures: pass
		else: init += [init_curvature for i in range(self.nstrengths)]
		return init

	def route_args_objective(self,*args):
		"""Route arguments incoming to the objective function to customers."""
		#! dev: this mimics the router for hel. avoid repetition here
		if self.model_style=='tension_protrusion': arglist_hel = ['sigma','kappa','gamma']
		elif self.model_style=='bending': arglist_hel = ['kappa']
		elif self.model_style=='tension': arglist_hel = ['sigma','kappa']
		elif self.model_style=='protrusion': arglist_hel = ['kappa','gamma']
		else: raise Exception
		# arguments must be in a list: hel then hosc then curvatures
		if self.equipartition_model=='default': arglist_equipartition = []
		elif self.equipartition_model=='harmonic_oscillators': arglist_equipartition = ['vibe']
		else: raise Exception
		if self.fix_curvatures: n_curvatures = 0
		else: n_curvatures = self.nstrengths
		n_expected_args = (n_curvatures+len(arglist_hel)+len(arglist_equipartition))
		if len(args)!=(n_curvatures+len(arglist_hel)+len(arglist_equipartition)):
			raise Exception('incorrect number of arguments to the objective. expected %d and got %d'%(
				n_expected_args,len(args)))
		hel_i = len(arglist_hel)
		heq_i = hel_i + len(arglist_equipartition)
		out = dict(args_equipartition=args[hel_i:heq_i],
			args_hel=tuple(list(args[0:hel_i])+([] if self.fix_curvatures else list(args)[heq_i:])),)
		return out

	def build_objective(self):
		"""
		Build an objective function from energies, 
		an equipartition model, the residual, and the band.
		"""
		def objective(args):
			"""Route the parameters."""
			router = self.route_args_objective(*args)
			args_hel = router['args_hel']
			args_equipartition = router['args_equipartition']
			# it is useful to have this in one place in the code to avoid repetition
			# methodology: we bin and then filter (but you could do this the other way around)
			error = self.residual(
				self.binner.bin(self.equipartition(*args_equipartition))[self.band],
				self.binner.bin(self.hel(*args_hel))[self.band])
			return error
		# keep the objective function with the instance
		self.objective = objective

	def optimize(self):
		"""Run the optimizer."""
		init = self.build_initial_conditions()
		self.opt = QuickOpt(self.objective,init=init)

if __name__=='__main__':

	# use replot() to rerun main block when adjusting the plots
	# note that we removed __replotting__ because absent from omnicalc.py 

	# switches
	do_demo = 0
	do_survey = 1
	do_spectra_survey_debug = 1
	do_compute_landscape = 0

	# constants
	halfsweep = np.concatenate((
		np.arange(0.01,0.1+0.01,0.01),
		np.arange(0.2,1.0+0.1,0.1),
		np.arange(2.,10.+1.,1.)))
	curvatures_extreme = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))
	halfsweep = np.array([0.0,0.005,0.01,0.014,0.018,0.02,0.024,0.028,0.032,0.04,0.05])
	halfsweep_bigger = np.array([0.0,0.005,0.01,0.02,0.04,0.05,0.1,0.2,0.3,0.5,1.,5,10.])
	curvatures_legacy_symmetric = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))
	curvatures_legacy_symmetric_bigger = np.concatenate((
		halfsweep_bigger[::-1]*-1,[0.],halfsweep_bigger))

	if do_demo:
	
		sn = work.sns()[0]
		self = ucc = UndulationCurvatureCoupling(sn=sn,grid_spacing=0.5,
			model_style='bending',equipartition_model='default',binner_method='perfect')
		ucc.build_curvature_fields(
			mode='protein_dynamic',isotropy_mode='isotropic',extent=2.0)
		ucc.build_elastic_hamiltonian(0.0)
		ucc.build_equipartition()
		ucc.build_objective()
		ucc.optimize()

	#! compute landscape code moved here while testing survey at the end
	if do_survey:

		# compute the landscape
		if do_compute_landscape:

			if 'jobs' not in globals():
				ucc = UndulationCurvatureCoupling(sn=sn,grid_spacing=0.5,q_cut=10.0,
					model_style='bending',equipartition_model='harmonic_oscillators',
					oscillator_reverse=oscillator_reverse,positive_vibe=positive_vibe,
					binner_method=binner_method,subtract_protrusions=subtract_protrusions)
				jobs,fitted = {},{}
				start = time.time()
				for hnum,hypo in enumerate(hypos[:]):
					status('searching %d/%d'%(hnum,len(hypos)),i=hnum,
						looplen=len(hypos),tag='optimize',start=start,refresh=False)
					extent,curvature = hypo['extent'],hypo['curvature']
					ucc.build_curvature_fields(mode='protein_dynamic',
						isotropy_mode='isotropic',extent=extent)
					ucc.build_elastic_hamiltonian(curvature)
					ucc.build_equipartition()
					ucc.build_objective()
					ucc.optimize()
					# collect the results
					jobs[(extent,curvature)] = ucc.opt.fit.fun
					fitted[(extent,curvature)] = ucc.opt.fit.x

			"""
			Make a figure for the error landscapes.
			"""
			figsize = (12,10)
			contour_interp_pts = 100
			contour_line_skip = 4
			contour_nlevels = 100
			under_color = 'm'
			figname = 'fig.ucc.review.%s.pv%d.ro%d.cs%d.b%d.sp%d'%(sn,
				positive_vibe,oscillator_reverse,curvature_sweep_number,
				binners.index(binner_method),subtract_protrusions)

			figsize = (12,8)
			axes,fig = square_tiles(2,figsize,favor_rows=True)
			ax = axes[0]
			raw = np.array([[jobs[(e,c)] for e in extents] for c in curvatures])
			kwargs = dict(extent=[min(curvatures),max(curvatures),
				min(extents),max(extents)],aspect=(curvatures.ptp()/extents.ptp()))
			kwargs = dict(extent=[0,len(curvatures),0,len(extents)],
				aspect=(float(len(curvatures))/len(extents)))
			ax.imshow(raw.T,origin='lower',interpolation='nearest',**kwargs)
			ax.set_xticks(np.arange(len(curvatures))+0.5)
			ax.set_yticks(np.arange(len(extents))+0.5)
			ax.set_xticklabels(['%.3f'%i for i in curvatures],rotation=90)
			ax.set_yticklabels(['%.1f'%i for i in extents])
			ax = axes[1]
			error_min = raw.min()
			error_max = raw.ptp()/2.+raw.min()
			contour_line_max = raw.ptp()/4.+raw.min()
			curvature_extent_error = np.array([(c,e,raw[cc,ee]) 
				for cc,c in enumerate(curvatures) for ee,e in enumerate(extents)])
			c0,c1 = min(curvatures),max(curvatures)
			e0,e1 = min(extents),max(extents)
			finex = np.linspace(c0,c1,contour_interp_pts)
			finey = np.linspace(e0,e1,contour_interp_pts)
			grid_x,grid_y = np.meshgrid(finex,finey)
			errormap = scipy.interpolate.griddata(
				curvature_extent_error[:,:2],
				curvature_extent_error[:,2],
				(grid_x,grid_y),method='cubic')
			levels = np.linspace(error_min,error_max,contour_nlevels)
			cs = ax.contourf(grid_x,grid_y,
				errormap,levels=levels,vmax=error_max,vmin=error_min,
				extend='both',origin='lower',lw=2,zorder=3,cmap=mpl.cm.jet)
			cs.cmap.set_over('w')
			if under_color: cs.cmap.set_under(under_color)
			levels_contour = levels[np.where(levels<=contour_line_max)][::contour_line_skip]
			if False: cs_lines = ax.contour(grid_x,grid_y,errormap,vmax=error_max,
				vmin=error_min,levels=levels_contour,
				extend='both',origin='lower',linewidths=0.5,colors='k',zorder=4)
			ax.set_aspect(curvatures.ptp()/extents.ptp())
			picturesave(figname,work.plotdir,backup=False,version=False,meta={})

	if do_survey:

		# settings
		positive_vibe = False
		oscillator_reverse = False
		#! dev: the following should be retired
		subtract_protrusions = False
		curvature_sweep_number = 2
		equipartition_model = ['default','harmonic_oscillators'][1]
		# the default preferred method is the perfect binner
		binner_method = ['explicit','perfect','blurry'][1]
		# decide which physical parameters to fit
		model_style = ['bending','tension','tension_protrusion','protrusion'][1]
		# the default method ("standard") uses log residuals
		residual_name = ['standard','linear','modelfix'][2]
		# the default equipartition model simply fits to 1/2 kBT and is useful for checking kappa 
		#   without any curvature to ensure that this fitting method matches the typical q^{-4} fit
		# the binner method called "perfect" will average valuse that occur at each unique wavevector
		#   which creates fewer wavevectors in our plot and also somewhat evens out the distribution of
		#   wavevectors along the x-axis. the perfect method treats each wavevector magnitude as a unique
		#   mode, and for this reason we think it is the appropriate measure. the blurry method bins the
		#   wavevectors into finite bins, thus collapsing the number of wavevectors even further. for
		#   CGMD simulations with the upward "swoosh" shape, the blurry binner produces a higher kappa
		#   as expected because it counteracts the extra weight from the higher density of wavevectors at
		#   higher values, which pulls the kappa down because the swoosh goes up
		# note on tension above. in the v650 (4xENTH) reference simulation, we see that adding tension
		#   causes kappa to drop slightly from 20.4 to 19.0 kBT, presumably because the intensities at 
		#   lower wavevectors are somewhat lower than expected, indicating a slight surface tension,
		#   however the inclusion of sigma in these calculations does not affect the results very
		#   much because we are well below the crossover
		# residuals should be standard, which takes the residuals in the log space
		# note that the linear method would have a more profound effect if we were 
		#   fitting hqhq to q^{-4} directly (see undulate.py) because then the linear residuals
		#   would shrink with higher wavevectors and hence lower hqhq. in the energy space, switching
		#   to linear residuals has much less effect because the intensities are all flanking a
		#   mostly uniform value (both the oscillators and equipartition are still somewhat constant)
		#   but as we noted above, the log residuals are more accurate because the fit should respect
		#   the intensities across several orders of magnitude in hqhq	

		sn = work.sns()[0]
		from codes.hypothesizer import hypothesizer
		# various curvature sweeps
		curvatures = [
			curvatures_extreme,
			curvatures_legacy_symmetric,
			curvatures_legacy_symmetric_bigger,
			][curvature_sweep_number]
		binners = ['explicit','perfect','blurry']
		extents = np.array([0.25,0.5,1.0,2.0,3.0,4.0,8.0])
		# hypotheses are built in argument-order so pick extent first 
		#   since that is slowest and gets recomputed the least at the front
		hypos = hypothesizer(*({'route':['extent'],'values':extents},
			{'route':['curvature'],'values':curvatures}))

		# new spectra figure for debugging
		if do_spectra_survey_debug:

			# test the no-curvature case
			extent,curvature = 1.0,0.0
			if 'ucc' not in globals():
				self = ucc = UndulationCurvatureCoupling(
					sn=sn,grid_spacing=0.5,q_cut=1.,
					model_style=model_style,
					equipartition_model=equipartition_model,
					oscillator_reverse=oscillator_reverse,
					positive_vibe=positive_vibe,
					binner_method=binner_method,
					subtract_protrusions=subtract_protrusions,
					residual_name=residual_name)
				ucc.build_curvature_fields(
					mode='protein_dynamic',
					isotropy_mode='isotropic',
					extent=extent)
				ucc.build_elastic_hamiltonian(curvature)
				ucc.build_equipartition()
				ucc.build_objective()
				ucc.optimize()

			# PLOT: survey the fitting procedure
			# in the deprecated validation procedure we check that changing from the q4 spectrum
			#   to energy in a fit without tension, protrusions, oscillator, or curvature, gives us
			#   the same spectrum that we originally fit things to. this is useful for checking that 
			#   the kappa fit for the UCC method matches a typical q4 spectrum fit
			do_survey_validate = False
			if do_survey_validate:
				fig = plt.figure(figsize=(8,12))
				axes = [fig.add_subplot(311),fig.add_subplot(312),fig.add_subplot(313)]
			else: 
				fig = plt.figure(figsize=(8,10))
				axes = [fig.add_subplot(211),fig.add_subplot(212)]
			args_parse = ucc.route_args_objective(*ucc.opt.fit.x)
			args_hel = ucc.route_args_hel(*args_parse['args_hel'])
			for key in ['kappa','sigma','gamma']:
				globals()[key] = args_hel.get(key,0.0)
			hel = ucc.hel(*args_parse['args_hel'])
			hosc = ucc.equipartition(*args_parse['args_equipartition'])
		
			# left panel shows energy spectra
			ax = axes[0]
			ax.set_ylim((0.1,10))
			ax.set_xlim((0.05,2))
			ax.axhline(ucc.energy_per_mode,c='k',lw=1)
			ax.set_ylabel('energy ($k_BT$)')
			ax.set_xlabel('wavevector (${nm}^{-1}$')
			# plot the explicit energy spectrum
			ax.plot(ucc.qs_raw,hel,'.',c='b',lw=0,
				label='$H_{el}$ ($\kappa=%.1f k_BT)$'%kappa)
			# plot the perfect binner results
			if residual_name=='modelfix':
				energy_apparent = ucc.binner.bin(hel)*ucc.binner.bin(hosc)
				ax.plot(ucc.qs,energy_apparent,'-',c='k',zorder=2,lw=4)
			# standard method follows
			else: 
				ax.plot(ucc.qs,ucc.binner.bin(hel),'-',c='k',zorder=2,lw=4)
			ends_check = np.argsort(ucc.qs_raw)[:2] 
			print('[STATUS] the left two intensities are: %s'%str(hel[ends_check]))
			label_hosc = {
				'default':'equipartition',
				'harmonic_oscillators':'oscillators'}[equipartition_model]
			ax.plot(ucc.qs_raw,hosc,'.',c='r',lw=0,label=label_hosc)
			if 0: ax.plot(binner.bin(self.qs),binner.bin(
				logsum(1.0,logdiff(self.hqhq,model_basic))),
					'-',c='m',lw=1,zorder=10,label='corrected')
			if 0: ax.plot(binner.bin(self.qs),binner.bin(model_basic)*q_raw**4*kappa*area,
				'.-',c='k',lw=1,label='basic')

			# right panel shows h2q
			ax = axes[1]
			ax.set_ylabel('$\|{h_q}{h_q}\|$ (${nm}^{-2}$)')
			ax.set_xlabel('wavevector (${nm}^{-1})$')
			model_q4 = 1./(self.area * kappa * self.qs**4)
			# plot the binned undulation spectrum
			# if you use the perfect binner, the ucc.qs are the reduced wavevectors
			# note that we are using ucc.termlist[0] which is hqhq computed from the complex fields
			ax.plot(ucc.qs,ucc.binner.bin(ucc.termlist[0]),'.',c='b',lw=1,label='observed')
			# plot the fitted line
			ax.plot(ucc.qs[ucc.band],model_q4[ucc.band],'-',c='r',lw=2,zorder=5,
				label='fit ($\kappa=%.1f{k_{B}T}$)'%kappa)

			if do_survey_validate:
				ax = axes[2]
				# plot the first coupled term (i.e. hqhq) in equation 23
				# use termlist not hqhq so we retire
				#   energy_this = kappa * ucc.area * ucc.hqhq * ucc.qs_raw**4
				energy_this = kappa * ucc.area * ucc.termlist[0] * ucc.qs_raw**4
				print('[CHECK] two left values are: %s'%str(energy_this[ends_check]))
				ax.plot(ucc.qs_raw,energy_this,'.',c='b',lw=1,
					label='observed ($\kappa=%.1f$)'%fitted[0])
				ax.plot(ucc.qs,ucc.binner.bin(energy_this),'-',c='k',zorder=2,lw=2)
				ax.axhline(1.0,c='k',lw=1)
				axes[2].legend()

			# formatting
			for ax in axes:
				ax.set_xscale('log')
				ax.set_yscale('log')
				ax.axvline(ucc.q_cut,c='k',lw=1)
			axes[0].legend()
			axes[1].legend()
			picturesave('fig.debug.v13',work.plotdir)


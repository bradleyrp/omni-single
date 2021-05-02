#!/usr/bin/env python

import os,sys,time
import numpy as np
from joblib import Parallel,delayed
from base.tools import status,framelooper
import multiprocessing as mp
machine_eps = eps = np.finfo(float).eps
import scipy
import scipy.optimize
from omnicalc import store,load

#---backwards compatibility
show_optimization_log = False
tweak = {}

###---ENERGY FUNCTIONS

def energy_form(energy_raw,q_raw,kappa,gamma,vibe): 
	"""Objective function which tests the basic energy function against the harmonic oscillators."""
	#---! excessive epsilons seem dangerous so the optimization is worth checking later
	return (energy_raw(kappa,gamma))/((abs(vibe)+eps)*q_raw)/(eps+1./(np.exp((abs(vibe)+eps)*q_raw)-1+eps))

def optimize_hypothesis(energy,residual,band,init_cond):
	"""Optimize the energy function within a particular range of wavevectors."""
	fit = scipy.optimize.fmin(residual(energy,band=band),x0=(tuple(init_cond),),disp=show_optimization_log)
	error = residual(energy,band=band)(fit)
	result = {'kappa':fit[0],'gamma':fit[1],'vibe':fit[2],'error':error}
	return result

###---TRANSFORMS

def perfect_collapser(xs,ys,trim=False):
	"""
	Return two arrays for a "perfect" collapse.
	"""
	xsort = np.array(np.sort(list(set(xs))))
	if trim: xsort = xsort[1:]
	inds = np.argmax(np.array([xs==xsort[i] for i in range(len(xsort))]),axis=0)
	if type(ys)==type(None): col = None
	else: col = np.array([np.mean(ys[np.where(inds==i)]) for i in range(len(xsort))])
	coly = np.array([np.mean(ys[np.where(inds==i)]) for i in range(len(xsort))])
	#---! previously returned xsort for some reason?
	return coly

def blurry_binner(xs,ys,bin_width=0.05,trim=False):
	"""
	Group wavevectors by bins.
	"""
	blurred = (xs/bin_width).astype(int)
	xsort = np.array(np.sort(list(set(blurred))))
	if trim: xsort = xsort[1:]
	inds = np.argmax(np.array([(xs/bin_width).astype(int)==xsort[i] for i in range(len(xsort))]),axis=0)
	if type(ys)==type(None): coly = None
	else: coly = np.array([np.mean(ys[np.where(inds==i)]) for i in range(len(xsort))])
	colx = np.array([np.mean(xs[np.where(inds==i)]) for i in range(len(xsort))])
	return colx

def explicit_binner(xs,ys):
	"""..."""
	xsort = np.array(np.sort(list(set(xs))))
	return ys[ysort]

def filter_wavevectors(q,low,high):
	"""
	Restrict attention to wavevectors in a particular band.
	"""
	filtered, = np.where(np.all((q<=high,q>=low),axis=0))
	return filtered

###---CURVATURE FIELDS

def fftwrap(dat): 
	"""
	Simple wrapper for the 2D FFT algorithm in scipy.
	"""
	return np.fft.fft2(np.array(dat))

def fft_field(dat):
	"""
	Take the FFT of a regular grid.
	"""
	#---take the norm of the FFT
	factor = np.product(dat.shape[1:])
	hqs_complex = np.array([fftwrap(i)/factor for i in dat])
	#---hard coding the right method here!
	fft_flag = 'complex'
	if fft_flag == 'real': hqs = np.real(hqs_complex)
	elif fft_flag == 'absolute': hqs = np.absolute(hqs_complex)
	elif fft_flag in ['complex','complex_real']: hqs = hqs_complex
	else: raise
	return hqs

def gauss2d(grid,**kwargs):
	"""
	Compute 2D Gaussian over a grid (faster than gauss2d).
	"""
	st = time.time()
	th = kwargs['theta']
	c0 = kwargs['curvature']
	sx = kwargs['sigma_a']
	sy = kwargs['sigma_b']
	z0 = kwargs.get('z0',0.0)
	x0,y0 = kwargs['x0'],kwargs['y0']
	x,y = np.transpose(np.reshape(grid,(np.product(np.shape(grid)[:2]),2)))
	if c0=='0': return zeros(np.shape(grid)[:2])
	t1 = -((x-x0)*np.cos(th)+(y-y0)*np.sin(th))**2/2./sx**2
	t2 = -(-(x-x0)*np.sin(th)+(y-y0)*np.cos(th))**2/2./sy**2
	field = np.reshape(z0+c0*np.exp(t1)*np.exp(t2),np.shape(grid)[:2])
	return field

def construct_curvature_field(**kwargs):
	"""
	Create a single curvature field over a regular grid over a fixed spatial extent.
	Note that we have a ceiling on curvature for added dimples.
	!!!!!!! PBCs under development.
	"""
	#---grid dimensions
	m,n = kwargs.pop('mn')
	#---default centers is in the ahem center otherwise a list of positions rescaled to box vectors
	centers = kwargs.get('centers',np.array([[0.5,0.5]]))
	curvature = kwargs.pop('curvature')
	vecs = kwargs.pop('vecs')
	#---lenscale is unity unless we wish to rescale other simulation types
	lenscale = kwargs.pop('lenscale',1.0)
	#---construct the grid over the box vectors under PBCs
	grid = np.array([[[i,j] for j in np.linspace(0,3*vecs[1]/lenscale,3*n)] 
		for i in np.linspace(0,3*vecs[0]/lenscale,3*m)])
	field = np.zeros(np.shape(grid)[:2])
	#---we apply PBCs here. this might slow things down
	#---previously PBCs were only implemented on the field (!) not the centers
	offsets = np.concatenate(np.transpose(np.meshgrid(np.arange(3),np.arange(3))))
	centers_pbc = np.array([c+offset for c in centers for offset in offsets])
	#---loop over dimples with replicates in PBCs
	for center in centers_pbc:
		field += gauss2d(grid,x0=vecs[0]*center[0],y0=vecs[1]*center[1],
			curvature=1.0,**kwargs)
	#---ceiling on curvature in case multiple inducers overlap
	#---! NOTE
	field[field>1.0] = 1.0
	#---remove duplicates under PBCs
	field_one = field[m:2*m,n:2*n]
	return field_one

def construct_curvature_fields_trajectory(**kwargs):
	"""
	Generate a curvature field for a whole trajectory. 
	This function applies the "dynamics mode" logic.
	"""
	nprocs = kwargs.get('nprocs',4)
	######## #---need the data
	#global data
	#data = kwargs.get('data')
	#---each field has the following attributes
	fields = {'scale_factor':0.0,'framewise':False,'fields':[]}
	#---incoming vectors for the full trajectory
	vecs = kwargs.pop('vecs')
	mapping = kwargs.pop('mapping')
	motion = kwargs.pop('motion')
	sn = kwargs['sn']
	nframes = int(data['undulations'][sn]['data']['nframes'])
	start = time.time()
	if mapping == 'single' and motion == 'static':
		#---take the mean box vectors and draw a single dimple in the center 
		#---the dimple is normalized to a peak of unity and can be rescaled
		vecs_mean = np.mean(vecs,axis=0)
		fields['fields'] = [construct_curvature_field(vecs=vecs_mean,**kwargs)]
	elif mapping == 'protein' and motion == 'static':
		vecs_mean = np.mean(vecs,axis=0)
		#---protein_abstractor saves points_all which is frames by proteins by atoms by xyz
		centers = data['protein_abstractor'][sn]['data']['points_all'].mean(axis=2).mean(axis=0)/vecs_mean
		fields['fields'] = [construct_curvature_field(vecs=vecs_mean,centers=centers,**kwargs)]	
	elif mapping == 'single' and motion == 'dynamic':
		collect = Parallel(n_jobs=nprocs,verbose=0)(
			delayed(construct_curvature_field)(vecs=vecs[fr],**kwargs)
			#for fr in framelooper(nframes,start=start))
			for fr in range(nframes))
		fields['fields'] = collect
	elif mapping == 'protein' and motion == 'dynamic':
		#---reworked the protein points using mean to get the right shape out of points_all
		collect = Parallel(n_jobs=nprocs,verbose=0)(
			delayed(construct_curvature_field)(vecs=vecs[fr],
				centers=data['protein_abstractor'][sn]['data']['points_all'][fr].mean(axis=1)/vecs[fr],
				**kwargs)
			for fr in range(nframes))
		fields['fields'] = collect
	else: raise Exception('unclear dynamics mode')
	return fields

###---TESTING HYPOTHESES

def residual(energy,*args,**kwargs):
	"""
	Decorate the energy function to create a residaul.
	"""
	fmin_format = kwargs.get('fmin_format',False)
	band = kwargs.get('band',slice(None,None))
	#---default form comes from "tweak" in globals 
	#---overrides are possible during debugging/comparison in pipeline_curvature_coupling_single.py
	residual_form = kwargs.get('residual_form',tweak.get('residual_form','log'))
	if residual_form == 'log':
		def residual_banded(arguments):
			energies = energy(*arguments)[band]
			return sum(np.log10(energies.clip(min=machine_eps))**2)/float(len(energies))
	elif residual_form == 'subtract':
		def residual_banded(arguments):
			energies = energy(*arguments)[band]
			return sum((energies)**2)/float(len(energies))
	elif residual_form == 'linear':
		def residual_banded(arguments):
			energies = energy(*arguments)[band]
			return sum((energies.clip(min=machine_eps)-1)**2)/float(len(energies))
	else: raise
	return residual_banded

def evaluate_hypothesis(hypo,sessions,full_output=False,preloaded_curvature=False,debug=False):
	"""
	Perform a single hypothesis test and return the result to the 
	callback which will save to the database.
	"""
	global Field,Hypothesis
	#---settings
	sn = hypo['sn']
	hypo_full = Hypothesis(**hypo)
	match = sessions['hypothesis'].query(Hypothesis).filter_by(**hypo_full.base()).one()
	#---get the heights and vectors from "memory"
	hqs = memory[(sn,'hqs')]
	nframes = len(hqs)
	#---get the curvature fields
	#---! later we will manage the memory externally?
	hypo_cf = Field(**hypo)
	match_cf = sessions['field'].query(Field).filter_by(**hypo_cf.dict()).one()
	fn_cf = namer_cf(match_cf.id)
	key_cf = ('curvature',match_cf.id)
	if key_cf not in memory and not preloaded_curvature: 
		memory[key_cf] = load(os.path.basename(fn_cf),rootdir_cf)
	curvature_fields = memory[key_cf]['fields']
	#---rescale the field to the right curvature
	cfs = np.array([hypo['curvature']*c for c in curvature_fields])
	#---for the "static" motion, we repeat the same curvature field over all frames
	if hypo['motion']=='static':
		if not len(cfs)==1: 
			raise Exception('motion is static but the pre-computed curvature fields have %s'%len(cfs))
		cfs = [cfs[0] for fr in range(nframes)]
	vecs = memory[(sn,'vecs')]
	mn = hqs.shape[1:]
	#---prepare the functions
	kwargs = dict(cfs=cfs) 
	if 'fallback' in hypo_full.dict() and hypo_full.dict()['fallback'] != 'none': 
		kwargs['fallback'] = True
	parts = couplecalc(hqs,vecs,**kwargs)
	binner = parts['binner']
	energy = parts['energy']
	q = parts['wavevectors']
	#---set the band
	#---! this needs added to the database
	band = filter_wavevectors(q,low=0.0,high=tweak.get('hicut',1.0))
	#---send to the optimizer
	init_cond = tweak.get('init',[20,0.0,10.0])
	result = optimize_hypothesis(energy,residual,band,init_cond)
	#---! on porting to the freshest omnicalc ca 2017.6.6 we find kappa not being saved
	#---! ...the following fix makes this whole function read kind of bad
	parts.update(**dict([(k,v) for k,v in result.items() if k in 'kappa gamma vibe error'.split()]))
	#---we save the parts for later using the record numbers
	parts['best_energy'] = np.array(energy(result['kappa'],result['gamma'],result['vibe']))
	result = dict([(key,val) for key,val in parts.items() if key not in ['binner','energy']])
	store(obj=result,name=os.path.basename(namer_cc(match.id)),path=rootdir_cc,
		attrs=hypo_full.base(),verbose=False)
	#---update the database
	for key,val in result.items(): match.__setattr__(key,val)
	sessions['hypothesis'].commit()
	if debug:
		import ipdb;ipdb.set_trace()

def couplecalc(hqs,vecs,**kwargs):
	"""
	Kernel of the curvature coupling calculation.
	"""
	#---default gets the functional form from tweak but could be overriden
	energy_functional = kwargs.get('energy_functional',energy_form)

	fallback = kwargs.get('fallback',False)
	lenscale = kwargs.get('lenscale',1.0)
	binner = kwargs.get('binner','explicit')
	#---if no incoming curvatures then we use the zero field
	#---identified a major error here: cqs = kwargs.get('cqs',fft_field(np.zeros(hqs.shape)))
	#---note the misleading naming
	cfs = kwargs.get('cfs',np.zeros(hqs.shape))
	cqs = fft_field(np.array(cfs))

	#---identify the wavevectors
	nframes = len(vecs)
	m,n = mn = np.shape(hqs)[1:]
	Lx,Ly = np.mean(vecs,axis=0)[:2]
	q2d = lenscale*np.array([[np.sqrt(
		((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
		((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
		for j in range(0,n)] for i in range(0,m)])
	q_raw = np.reshape(q2d,-1)[1:]
	area = (Lx*Ly/lenscale**2)
	binner_function = {
		'blurry':blurry_binner,
		'perfect':perfect_collapser,
		'explicit':lambda x:x,
		}[binner]

	#---on the fallbacks we might have different numbers of frames
	if len(hqs) != len(cqs):
		assert fallback
		cqs = cqs[np.array([i%len(cqs) for i in range(len(hqs))])]
		assert len(hqs)==len(cqs)

	#---crucial steps
	if tweak.get('fft_flag','complex')=='complex': 
		def multipliers(x,y):
			#---use the fact that all real-valued functions are symmetric under FFT
			return x*np.conjugate(y)
		termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
	else: raise
	#---take the ensemble average, convert the matrix into a vector, and drop the zeroth-order term
	termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
	#---we confirm that the incoming functions obey the FFT before casting to real below
	for t in [termlist[0],termlist[1]+termlist[2],termlist[3]]:
		assert np.all(np.absolute(np.imag(t))<machine_eps)
	termlist = [np.real(k) for k in termlist]
	
	#---equation 23
	def energy_raw(kappa,gamma):

		"""
		Define the energy.
		This function inherits the wavevector.
		Should the function be generic, or is it necessary to pack it up to optimize it later?
		"""

		signterm = tweak.get('inner_sign',-1.0)
		#---included the half-factor for hcv3 v28,29 and removed for v30
		curv = (kappa*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
			+signterm*termlist[2]*q_raw**2+termlist[3])
			#---removed square on ther first term in front of the tension term
			+gamma*area*(termlist[0]*q_raw**2))
		return curv

	def energy(kappa,gamma,vibe):

		"""
		Corrected energy function.
		"""

		spectrum = energy_functional(energy_raw,q_raw,kappa,gamma,0)
		return spectrum

	#---export
	out = {
		'area':area,
		'q2d':q2d,
		'wavevectors':q_raw,
		'energy':energy,
		'binner':binner_function,
		'grid':(m,n),
		#---HACKED then un-hacked because it was trying to send the function to store
		#'energy_raw':energy_raw,
		}

	#---return the raw energy function for debugging
	out_raw = kwargs.get('out_raw',False)
	if out_raw: out['energy_raw'] = energy_raw
	return out

###---MULTIPROCESSING functions for the "manual" hypothesis testing

def manyjob_worker(queue_hypothesis,session_constructors,kwargs={}):
	"""
	Individual worker which evaluates hypotheses (and creates necessary curvature fields).
	"""
	sessions = dict([(key,val()) for key,val in session_constructors.items()])
	while True:
		if queue_hypothesis.empty(): break
		hypo = queue_hypothesis.get(True)
		result = evaluate_hypothesis(hypo,sessions=sessions,**kwargs)
	for session in sessions.values(): session.close()

def manyjob(function,queue,objects,session_classes,kwargs=None,single=False):
	"""
	Run the computation many times over.
	"""
	kwargs = {} if not kwargs else kwargs
	#---single processor for testing
	njobs = len(objects)
	#---the single 
	if single:
		sessions = dict([(key,val()) for key,val in session_classes.items()])
		for hypo in objects[:]:
			status("solving in serial: %s"%str(hypo),tag='compute')
			#---remember that function is the worker that wraps evaluate_hypothesis
			evaluate_hypothesis(hypo,sessions,debug=True,**kwargs)
			print('debugging')
			import ipdb;ipdb.set_trace()
	#---multiprocessing
	else:
		interrupted = False
		try:
			for hypo in objects: queue.put(hypo)
			pool = mp.Pool(4,function,(queue,session_classes,kwargs))
			pool.close()
			start = time.time()
			while not queue.empty():
				status('hypothesizing',i=njobs-queue.qsize()+1,
					looplen=njobs,tag='compute',start=start)
				time.sleep(1)
			pool.join()
		except KeyboardInterrupt:
			print("[STATUS] interrupted!")
			pool.terminate()
			pool.join()
			interrupted = True
		if not interrupted:
			status('computations complete in %.1fmin'%((time.time()-start)/60.),tag='status')

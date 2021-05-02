#!/usr/bin/env python

"""
CURVATURE UNDULATION COUPLING
version 4 (redeveloped and audited version 3)
"""

import importlib,datetime,time
import scipy
import scipy.optimize
import scipy.interpolate
from codes.hypothesizer import hypothesizer

import matplotlib.patheffects as path_effects
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

machine_eps = eps = np.finfo(float).eps

# added here on 2021.04.27
# moved from plot-curvature_undulation_coupling_pixel.py on 2021.04.17
import copy
def collect_upstream_calculations_over_loop():
	"""
	Some plotting and analysis benefits from checking all calculations in an upstream loop (which is 
	contrary to the original design of )
	!!! This script is a candidate for inclusion in omnicalc.py
	"""
	global plotname,work
	#---! move to a separate function
	plotspecs = work.metadata.plots.get(plotname,work.metadata.calculations.get(plotname,{})).get('specs',{})
	calcname = plotspecs.get('calcname',plotname)
	#---load the canonical upstream data that would be the focus of a plot in standard omnicalc
	#---! load the upstream data according to the plot. note that this may fail in a loop hence needs DEV!
	try: data,calc = plotload(plotname)
	except:
		data,calc = None,None
		status('failed to load a single upstream calculation however this plot script has requested '
			'all of them so we will continue with a warning. if you have downstream problems consider '
			'adding a specific entry to plots to specify which item in an upstream loop you want',
			tag='warning')
	#---in case there is no plot entry in the metadata we copy it
	if plotname not in work.metadata.plots: work.metadata.plots[plotname] = copy.deepcopy(work.metadata.calculations[calcname])
	#---load other upstream data
	#---get all upstream curvature sweeps
	upstreams,upstreams_stubs = work.calcs.unroll_loops(work.metadata.calculations[calcname],return_stubs=True)
	#for u in upstreams_stubs: u['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for unum,upstream in enumerate(upstreams_stubs):
		#---use the whittle option to select a particular calculation
		dat,cal = plotload(calcname,whittle_calc={calcname:upstream['specs']})
		tag = upstreams_stubs[unum]['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[calcname][sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
	#---singluar means the typical "focus" of the upstream calculation, plural is everything else
	return dict(datas=datas,calcs=calcs,data=data,calc=calc)

###
### imports preserved for posterity from drilldown
###

def blurry_binner(xs,ys,bin_width=0.05,trim=True,return_mapping=False,mode='snapped'):
	"""
	Improved over original with sensible binning and interpolation.
	Removed the index return strategy in favor of handling things internally.
	"""
	#---get the range of the independent variable
	x0,x1 = xs.min(),xs.max()
	#---snap to the bin width
	x0b,x1b = bin_width*np.floor(x0/bin_width),bin_width*np.ceil(x1/bin_width)
	#---develop canonical bin markers
	bins = np.arange(x0b,x1b+bin_width,bin_width)
	#---bins are represented in the middle
	xmid = (bins[1:]+bins[:-1])/2.
	#---get the bin positions
	xbin = np.floor(xs/bin_width).astype(int) - np.floor(xs.min()/bin_width).astype(int)
	#---check that everything is in the right bin
	try: 
		if not (np.all(xs<=bins[1:][xbin]) and np.all(xs>=bins[:-1][xbin])):
			import ipdb;ipdb.set_trace()
			raise Exception('binning failure!')
	except:
		import ipdb;ipdb.set_trace()
	#---reindex each position into bins
	indices = np.unique(xbin)
	reindex = [np.where(xbin==i)[0] for i in indices]
	#---snapped mode just uses the middle of each bin and takes the average of the constituents
	if mode=='snapped':
		y_out = np.array([ys[r].mean() for r in reindex])
		x_out = xmid[indices]
	#---! recommend developing an interpolation method here to shift laterally if off-center constituents
	else: raise 
	#---for each point, return the indices of other points in the same bin
	#---...note that this is essential for developing weights in the blurry_explicit
	#---...where the weights end up: weights = 1./np.array([len(i) for i in q_mapping])[band]
	if return_mapping: 
		mapping = [np.where(xbin==i)[0] for i in np.sort(xbin)]
		return x_out,y_out,mapping
	else: return x_out,y_out

def log_reverse(raw):
	return 10.**(-1.0*np.log10(raw))

def prepare_oscillator_function(reverse=False,positive_vibe=True):
	def oscillator_function(vibe,qs):
		#---vibration must be nonnegative
		if positive_vibe: vibe = np.abs(vibe)
		#---! no grease here
		raw = (vibe*qs)*(0.+1./(np.exp(vibe*qs)-1))
		if not reverse: return raw
		else: return log_reverse(raw)
	return oscillator_function

def prepare_residual(mode='standard'):
	if mode in ['standard','comp','subtractor']:
		def residual(hel,hosc): return np.mean((np.log10(hel)-np.log10(hosc))**2)
	elif mode=='alt':
		def residual(hel,hosc): 
			energies = hel/hosc
			return sum(np.log10(energies.clip(min=machine_eps))**2)/float(len(energies))
	else: raise Exception('unclear residual mode')
	return residual

def prepare_objective(
	hqs,curvature_fields,wavevectors,area,
	curvature_sum_function,fft_function,band,
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
	imaginary_mode = kwargs.pop('imaginary_mode','complex')
	fix_curvature = kwargs.pop('fix_curvature',None)
	oscillator_function_local = kwargs.pop('oscillator_function',None)
	subtractor = kwargs.pop('subtractor',None)
	#---handle master modes
	master_mode = kwargs.pop('master_mode','standard')
	signterm = {'standard':1.0,'alt':-1.0,'comp':True}.get(master_mode,1.0)
	positive_vibe = {'standard':True,'alt':False,'comp':True}.get(master_mode,True)
	#---! adding this back in
	positive_vibe = kwargs.pop('positive_vibe',False)
	reverse_oscillator = {'standard':True,'alt':False,'comp':False,'subtractor':False}[master_mode]
	#---get the oscillator if not explicit
	if oscillator_function_local==None:
		oscillator_function_local = prepare_oscillator_function(reverse=reverse_oscillator,
			positive_vibe=positive_vibe)
	residual = kwargs.pop('residual_function',
		prepare_residual(mode={'comp':'standard'}.get(master_mode,master_mode)))
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
		q_raw_binned_temp,_,q_mapping = blurry_binner(q_raw,q_raw,return_mapping=True)
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
		###!!! special debugging for redev use
		if len(hqs)!=len(cqs) and fix_curvature==0:
			termlist = [multipliers(hqs,hqs)]+[multipliers(np.zeros(hqs.shape),np.zeros(hqs.shape)) 
				for jj in range(3)]
		else:
			#---construct the terms in eqn 23
			termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
		#---reshape the terms into one-dimensional lists, dropping the zeroth mode
		termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
		#---! explain logic behind real/imaginary here
		if imaginary_mode=='real': termlist = [np.real(k) for k in termlist]
		elif imaginary_mode=='complex': termlist = [np.abs(k) for k in termlist]
		else: raise Exception('invalid imaginary_mode %s'%imaginary_mode)
		def objective(args,mode='residual',return_spectrum=False,debug=False):
			"""
			Fit parameters are defined in sequence for the optimizer.
			They are: kappa,gamma,vibe,*curvatures-per-dimple.
			"""
			if master_mode in ['standard','alt']:
				if len(args)>3: raise Exception('only send three arguments')
				#---the first three arguments are always bending rigitidy, surface tension, vibration
				kappa,gamma,vibe = args[:3]
				#---CONSTRAINTS are enforced here
				if positive_vibe: vibe = np.abs(vibe)
			###
			### PROTRUSION IS UNDER DEVELOPMENT HERE SEE debug_resurvey in the drilldown
			###
			elif master_mode=='comp':
				kappa,gamma,gamma_p,vibe = args[:4]
				gamma_p = np.abs(gamma_p)
				#---! assume no macro tension for now
				gamma = 0.0
				if positive_vibe: vibe = np.abs(vibe)
			elif master_mode=='subtractor':
				kappa,gamma = args[:2]
				kappa,gamma = np.abs(kappa),np.abs(gamma)
			else: raise Exception('invalid master_mode %s'%master_mode)

			#---constructing the elastic Hamiltonian based on the wavevectors
			if master_mode in ['standard','alt']: 
				hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma*area*(termlist[0]*q_raw**2))
			elif master_mode=='comp':
				hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma_p*area*(termlist[0]*q_raw**2))
			elif master_mode=='subtractor':
				hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma*area*(termlist[0]*q_raw**2))
			else: raise Exception('invalid master_mode %s'%master_mode)			
			if master_mode=='subtractor': hosc = subtractor + 1.0
			else: hosc = oscillator_function_local(vibe,q_raw)
			#---note that the band is prepared in advance above
			if binner_method=='explicit': ratio = hel
			elif binner_method=='perfect':
				q_binned,ratio,q_binned_inds = perfect_collapser(q_raw,hel)
			elif binner_method=='blurry':
				q_binned,ratio = blurry_binner(q_raw,hel)
			else: raise Exception('invalid binner_method %s'%binner_method)
			#----compute residuals with relevant wavevectors (in the band) and return
			if mode=='residual': 
				if type(weights)!=type(None): value = residual(weights*hel[band],weights*hosc[band])
				else: value = residual(hel[band],hosc[band])
			elif mode=='elastic': value = hel
			else: raise Exception('invalid residual mode %s'%mode)
			if debug:
				print('!!!!!!!!!')
				import ipdb;ipdb.set_trace()
			return value
		return objective

	def objective(args,mode='residual',return_spectrum=False,debug=False):
		"""
		Fit parameters are defined in sequence for the optimizer.
		They are: kappa,gamma,vibe,*curvatures-per-dimple.
		"""
		if master_mode in ['standard','alt']:
			if len(args)>3: raise Exception('only send three arguments')
			#---the first three arguments are always bending rigitidy, surface tension, vibration
			kappa,gamma,vibe = args[:3]
			#---CONSTRAINTS are enforced here
			if positive_vibe: vibe = np.abs(vibe)
			#---uniform curvatures are multiplexed here
			if ndrops_uniform!=0: curvatures = [args[3] for i in range(ndrops_uniform)]
			#---one curvature per field
			else: curvatures = args[3:]

		###
		### PROTRUSION IS UNDER DEVELOPMENT HERE SEE debug_resurvey in the drilldown
		###
		elif master_mode=='comp':
			kappa,gamma,gamma_p,vibe = args[:4]
			gamma_p = np.abs(gamma_p)
			#---! assume no macro tension for now
			gamma = 0.0
			if positive_vibe: vibe = np.abs(vibe)
			#---uniform curvatures are multiplexed here
			if ndrops_uniform!=0: curvatures = [args[4] for i in range(ndrops_uniform)]
			#---one curvature per field
			else: curvatures = args[4:]
		else: raise Exception('invalid master_mode %s'%master_mode)
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
		if master_mode in ['standard','alt']: 
			hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
				+signterm*termlist[2]*q_raw**2+termlist[3])
				+gamma*area*(termlist[0]*q_raw**2))
		elif master_mode=='comp':
			hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
				+signterm*termlist[2]*q_raw**2+termlist[3])
				+gamma*area*(termlist[0]*q_raw**2))
		elif master_mode=='subtractor':
			hel = (kappa/2.0*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
				+signterm*termlist[2]*q_raw**2+termlist[3])
				+gamma*area*(termlist[0]*q_raw**2))
		else: raise Exception('invalid master_mode %s'%master_mode)
		#---apply the vibration correction
		if master_mode=='subtractor': hosc = subtractor+1.0
		else: hosc = oscillator_function_local(vibe,q_raw)
		#---note that the band is prepared in advance above
		if binner_method=='explicit': ratio = hel
		elif binner_method=='perfect':
			q_binned,ratio,q_binned_inds = perfect_collapser(q_raw,hel)
		elif binner_method=='blurry':
			q_binned,ratio = blurry_binner(q_raw,hel)
		else: raise Exception('invalid binner_method %s'%binner_method)
		#----compute residuals with relevant wavevectors (in the band) and return
		if mode=='residual':
			if type(weights)!=type(None): value = residual(weights*ratio[band],weights*hosc[band])
			else: value = residual(hel[band],hosc[band])
		elif mode=='elastic': value = ratio
		else: raise Exception('invalid residual mode %s'%mode)
		if debug:
			print('!!!!!!!!!')
			import ipdb;ipdb.set_trace()
		return value

	#---return the decorated function
	return objective

###
### redev begins here
###

def stringer(x,p=5):
	"""A nice way to print numbers."""
	return (' '*(1*(x>0)))+('%.1e'%x 
		if (abs(x)<10**(-1*p) or abs(x)>10**p) else ('{:>%df}'%(p+2+1*(x<0))).format(x))

def callback(args,silent=False):
	"""Watch the optimization."""
	global stepno
	args_string = ' '.join([stringer(a) for a in args])
	output = (u'\r' if stepno>0 else '\n')+'[OPTIMIZE] step %s: %s'%(stepno,args_string)
	if ~silent:
		sys.stdout.flush()
		sys.stdout.write(output)
	stepno += 1

class QuickOpt:
	def __init__(self,objective,init,**kwargs):
		self.optimize_method = kwargs.pop('optimize_method','Nelder-Mead')
		self.scipy_optimize_function = kwargs.pop('scipy_optimize_function','minimize')
		self.silent = kwargs.pop('silent',False)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		self.stepno = 0
		#---get the objective function from globals
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
		output = (u'\r' if self.stepno>0 else '\n')+'[OPTIMIZE] step %s: %s'%(self.stepno,args_string)
		if not self.silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		self.stepno += 1

@autoload(plotrun)
def loader():
	"""Load data."""
	#---only load once
	if 'data' not in globals():
		#---begin loading sequence
		plotname = 'curvature_undulation_coupling'
		if plotname not in work.metadata.plots: raise Exception('add %s to the plots metadata'%plotname)
		plotspecs = work.metadata.plots[plotname].get('specs',{})
		calcname = plotspecs.get('calcname',plotname)
		#---new method for getting all upstream calculations in the loop
		#! combodat = work.collect_upstream_calculations_over_loop(plotname)
		#! rescued the following function
		combodat = collect_upstream_calculations_over_loop()
		data,datas,calcs = combodat['data'],combodat['datas'],combodat['calcs']
		#---we expect no optimization
		#---repopulate curvature fields
		for tag in datas:
			for sn in datas[tag]:
				datas[tag][sn]['cf'] = datas[tag][sn]['fields_unity'].sum(axis=1).mean(axis=0)
				datas[tag][sn]['cf_first'] = datas[tag][sn]['fields_unity'].sum(axis=1)[0]
		protein_abstractor_name = plotspecs.get('protein_abstractor_name','protein_abstractor')
		undulations_name = plotspecs.get('undulations_name','undulations')
		#---check for alternate loaders
		alt_loader = plotspecs.get('loader',None)
		if alt_loader:
			from base.tools import gopher
			postdat = gopher(alt_loader)(data=data)
		#---compile the necessary (default) data into a dictionary
		else:
			postdat = dict([(sn,dict(
				vecs=data[undulations_name][sn]['data']['vecs'].mean(axis=0),
				points_protein=data[protein_abstractor_name][sn]['data']['points_all']))
				for sn in work.sns()])
		#---end loading sequence
		#---export to globals
		global seepspace
		seepspace = 'data,datas,postdat,undulations_name,calcs,protein_abstractor_name'.split(',')
		for key in seepspace: globals()[key] = locals()[key]

def gopher(spec,module_name,variable_name,work=None):
	"""Load an external module. Useful for changing the workflow without changing the code."""
	#---note that this function was folded into the omnicalc codebase in omni/base/tools.py
	mod = importlib.import_module(spec[module_name])
	target = mod.__dict__.get(spec[variable_name],None)
	#---the database design might need work so we always export it
	mod.work = work
	if not target: raise Exception('add %s and %s to the specs'%(module_name,variable_name))
	return target

#---! is the fft_field repetitive with gopher?

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

def log_reverse(raw): return 10.**(-1.0*np.log10(raw))

class Binner:
	trim = True
	def __init__(self,mode,wavevectors,**kwargs):
		self.mode = mode
		self.bin_width = kwargs.pop('bin_width',0.05)
		if kwargs: raise Exception
		self.wavevectors = self.qs = wavevectors
		#---trivial case where all bins are unique
		if mode=='explicit': self.indexer = np.transpose([np.arange(len(self.qs))])
		elif mode=='blurry':
			blurred = (self.qs/self.bin_width).astype(int)
			self.indexer = [np.where(blurred==i)[0] for i in np.unique(blurred)]
			#---! should also save weights here?
		elif mode=='perfect':
			uniques,inds = np.unique(self.qs,return_inverse=True)
			self.indexer = [np.where(inds==i)[0] for i in np.unique(inds)]
			#---! should also save weights here?
		else: raise Exception
	def bin(self,values): 
		if self.mode=='explicit': return values[self.indexer].mean(axis=1)
		#---! is this step slow?
		else: return np.array([values[i].mean() for i in self.indexer])

def basic_undulation_fitter(qs,hqs,area,**kwargs):
	"""
	Fit the undulations with no curvature.
	"""
	#---directions
	fit_tension = kwargs.pop('fit_tension',True)
	fit_correction = kwargs.pop('fit_correction',True)
	fit_protrusion = kwargs.pop('fit_protrusion',True)
	#---settings
	q_min = kwargs.pop('q_min',0.0)
	q_cut = kwargs.pop('q_cut',4.0)
	binner_method = kwargs.pop('binner_method','explicit')
	bin_width = kwargs.pop('bin_width',0.05)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	#---! check this. using complex here. this is also performed by the build_hqhq function in the class
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
				(1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps)+(1.)/(gamma_p*q_raw**2+machine_eps))
		#---! CRITICAL ERROR:
		#---! ... if correct: osc = (vibe*q_raw+machine_eps)*(1./(np.exp(vibe*q_raw)-1+machine_eps))
		#---! ... INCORRECT RESULTS DUE TO FUDGE FACTORS!
		#---! removed to be safer 
		#---! ... if correct: osc = (vibe*q_raw+machine_eps)*(1./(np.exp(vibe*q_raw)-1)+machine_eps)
		if correct: osc = (vibe*q_raw)*(1./(np.exp(vibe*q_raw)-1))
		else: osc = 1.0
		return pure * osc
	def residual(a,b): return ((np.log10(a/b))**2).mean()
	binner = Binner(wavevectors=qs[band],mode=binner_method,bin_width=bin_width)
	def objective(args,protrusion=fit_protrusion,correct=fit_correction,tension=fit_tension):
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
		output = (u'\r' if stepno>0 else '\n')+'[OPTIMIZE] step %s: %s'%(stepno,args_string)
		if not silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		stepno += 1
	initial_conditions = (0.,20.,0.,0.)
	fit = scipy.optimize.minimize(objective,
		x0=initial_conditions,method='Nelder-Mead',callback=callback)
	fitted = model(qs,*fit.x)
	print('\n[OPTIMIZE] completed with error %s'%fit.fun)
	return dict(sigma=fit.x[0] if fit_tension else 0.0,band=band,
		kappa=np.abs(fit.x[1]),vibe=fit.x[3] if fit_correction else 0.0,
		gamma=fit.x[2] if fit_protrusion else 0.0,
		model=model)

def build_elastic_hamiltonian(
	area,inner_sign,qs,hqs,curvature_sum_function,curvature_fields,
	curvature_sum_method,fft_function,ndrops_uniform,fit_tension):
	"""
	A decorator which prepares the Hamiltonian. 
	No keyword arguments. Incoming data proliferates to several internal functions.
	"""

	def multipliers(x,y): 
		"""Multiplying complex matrices in the list of terms that contribute to the energy."""
		return x*np.conjugate(y)

	def curvature_multiplexer(curvatures):
		"""
		Scale the curvature fields by the curvature hypothesis and generate the terms in equation 23.
		"""
		#---merge curvature fields
		composite = curvature_sum_function(curvature_fields,curvatures,method=curvature_sum_method)
		cqs = fft_function(composite)
		###!!! special debugging for redev use
		if len(hqs)!=len(cqs) and all([c==0.0 for c in curvatures]):
			termlist = [multipliers(hqs,hqs)]+[multipliers(np.zeros(hqs.shape),np.zeros(hqs.shape)) 
				for jj in range(3)]
		else:
			#---construct the terms in eqn 23
			termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
		#---reshape the terms into one-dimensional lists, dropping the zeroth mode
		#---also convert complex to magnitudes here
		termlist = [np.abs(np.reshape(np.mean(k,axis=0),-1)[1:]) for k in termlist]
		return termlist

	def arghandle(args):
		"""
		Manage the correct number of arguments for these fits.
		"""
		n_core_args = 2 if fit_tension else 1
		#---uniform curvatures are multiplexed here
		if ndrops_uniform!=0: curvatures = [args[n_core_args] for i in range(ndrops_uniform)]
		#---one curvature per field if they are distinct
		else: curvatures = args[n_core_args:]
		return dict(kappa=args[0],gamma=args[1] if fit_tension else 0.0,curvatures=curvatures)

	#---the energy calls on qs,termlist,fit_tension
	def energy_functional(args,debug=False):
		"""
		Scaffold for the Hamiltonian.
		"""
		argdict = arghandle(args)
		kappa,gamma,curvatures = argdict['kappa'],argdict['gamma'],argdict['curvatures']
		termlist = curvature_multiplexer(curvatures=curvatures)
		#---handle macroscopic tension
		if fit_tension: gamma_term = gamma*area*(termlist[0]*qs**2)
		else: gamma_term = 0.0
		#---construct the energy
		hel = (kappa/2.0*area*(termlist[0]*qs**4+inner_sign*termlist[1]*qs**2+inner_sign*termlist[2]*qs**2+termlist[3])+gamma_term)
		if debug:
			import ipdb;ipdb.set_trace()
		return hel

	return energy_functional

class InvestigateCurvatureMicro:
	"""
	"""
	version = 4.0
	def __init__(self,**kwargs):
		"""Formulate and validate the hypothesis."""
		global data
		#---required arguments
		self.mode = kwargs.pop('mode')
		self.sn = kwargs.pop('sn')
		self.curvature_fields = kwargs.pop('curvature_fields')
		self.hqs = kwargs.pop('hqs')
		#---optional arguments
		self.review = kwargs.pop('review',False)
		self.q_min = kwargs.pop('q_min',0.0)
		self.q_cut = kwargs.pop('q_cut',1.0)
		self.bin_width = kwargs.pop('bin_width',0.05)
		self.binner_method = kwargs.pop('binner_method','explicit')
		self.inner_sign = kwargs.pop('inner_sign',1.0)
		self.frameslice = kwargs.pop('frameslice',None)
		self.fit_tension = kwargs.pop('fit_tension',True)
		self.legacy_spec = kwargs.pop('legacy_spec',{})
		self.fit_protrusion = kwargs.pop('fit_protrusion',True)
		self.fit_correction = kwargs.pop('fit_correction',True)
		self.ndrops_uniform = kwargs.pop('ndrops_uniform',0)
		self.scipy_optimize_function = kwargs.pop('scipy_optimize_function','minimize')
		self.scipy_optimize_method = kwargs.pop('scipy_optimize_method','Nelder-Mead')
		self.midplane_method = kwargs.pop('midplane_method','flat')
		self.curvature_sum_method = kwargs.pop('curvature_sum_method','mean')
		self.name = kwargs.pop('name','icm_%s'%
			datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M'))
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		#---! get vectors from outside
		self.vecs = data['undulations'][self.sn]['data']['vecs']
		#---subsample frames if necessary
		if type(self.frameslice)==type(None): pass
		elif type(self.frameslice)==int: 
			self.hqs = self.hqs[np.linspace(0,len(self.hqs)-1,self.frameslice).astype(int)]
		else: raise Exception('cannot process frameslice %s'%frameslice)
		#---heights in Fourier space
		self.dims = self.hqs.shape[-2:]
		wavevectors_out = formulate_wavevectors(vecs=self.vecs,dims=self.dims)
		self.qs,self.area = wavevectors_out['wavevectors'],wavevectors_out['area']
		#---route the mode to the handler
		if self.mode=='subtract_protrusion': self.subtract_protrusion()
		elif self.mode=='legacy': self.legacy()
		else: raise Exception('invalid mode %s'%self.mode)
		#---review the energy
		if self.review: self.review_energy()

	def build_hqhq(self):
		#---! check this
		self.hqhq = (np.abs(self.hqs).reshape((len(self.hqs),-1)).mean(axis=0)[1:])**2 

	def residual(self,a,b): return ((np.log10(a)-np.log10(b))**2).mean()
	def residual_clip(self,a,b): return np.mean(np.log10((a/b).clip(min=machine_eps))**2)

	###
	### LEGACY MODE
	def legacy(self):
		"""
		Recapitulating the legacy mode.
		"""
		#---we cannot fit protrusion and correction from heights in legacy mode
		self.fit_protrusion,self.fit_correction = False,False
		self.build_hqhq()
		self.undulation_fit = basic_undulation_fitter(
			fit_protrusion=self.fit_protrusion,fit_correction=self.fit_correction,
			hqs=self.hqs,qs=self.qs,area=self.area,q_cut=self.q_cut,q_min=self.q_min,
			fit_tension=self.fit_tension,binner_method=self.binner_method,bin_width=self.bin_width)
		self.band = self.undulation_fit['band']
		#---prepare the hamiltonian
		self.hamiltonian = build_elastic_hamiltonian(
			area=self.area,inner_sign=self.inner_sign,qs=self.qs,hqs=self.hqs,
			#---! call out to the curvature sum function above but should it be passed through?
			curvature_sum_function=curvature_sum_function,curvature_fields=self.curvature_fields,
			curvature_sum_method=self.curvature_sum_method,fft_function=fft_field,
			ndrops_uniform=self.ndrops_uniform,fit_tension=self.fit_tension)
		#---create the binner
		self.binner = Binner(wavevectors=self.qs[self.band],
			mode=self.binner_method,bin_width=self.bin_width)
		#---remap the residual style
		if self.legacy_spec.get('residual_style','standard')=='clip': self.residual = self.residual_clip
		safe = self.legacy_spec.get('safe',True)
		#---constraints on vibration correction
		def vibe_correct(vibe):
			if safe=='negative': v = -1*np.abs(vibe)
			elif safe==True: v = np.abs(vibe)
			elif safe==False: v = vibe
			else: raise Exception
			return v
		#---initialize the vibration correction
		def vibe_init():
			if safe=='negative': v = -5.0
			elif safe==True: v = 5.0
			elif safe==False: v = -5.0
			else: raise Exception
			return v
		self.vibe_correct,self.vibe_init = vibe_correct,vibe_init
		def oscillator(vibe):
			v = self.vibe_correct(vibe)
			if self.legacy_spec.get('grease',False): 
				osc = (v*self.qs+machine_eps)*(1./(np.exp(v*self.qs)-1)+machine_eps)
			else: osc = (v*self.qs)*(1./(np.exp(v*self.qs)-1))
			if self.legacy_spec.get('reverse_oscillator',False): return log_reverse(osc)
			else: return log_reverse(osc)
		#---! DEVELOPMENT VERSION FOR USE WITH THE OTHER ELASTIC HAMILTONIAN FROM curvature_coupling.py
		def oscillator_qs(vibe,qs):
			v = self.vibe_correct(vibe)
			if self.legacy_spec.get('grease',False): 
				osc = (v*self.qs+machine_eps)*(1./(np.exp(v*qs)-1)+machine_eps)
			else: osc = (v*qs)*(1./(np.exp(v*qs)-1))
			if self.legacy_spec.get('reverse_oscillator',False): return log_reverse(osc)
			else: return log_reverse(osc)
		self.oscillator_qs = oscillator_qs
		def objective(args):
			"""Optimize over kappa,sigma,vibe,curvature."""
			if self.fit_tension: (kappa,gamma,vibe),curvatures = args[:3],args[3:]
			else: (kappa,vibe),gamma,curvatures = args[:2],0.0,args[2:]
			vibe = self.vibe_correct(vibe)
			if self.fit_tension: args_hamiltonian = tuple([kappa,gamma]+[curvatures])
			else: args_hamiltonian = tuple([kappa]+[curvatures])
			return self.residual(self.binner.bin(self.hamiltonian(args_hamiltonian)[self.band]),
				self.binner.bin(oscillator(vibe))[self.band])
		def build_objective_macros(curvatures):
			"""Prepare the objective function for a specific curvature."""
			def objective_macros(args,debug=False):
				"""Objective function with fixed curvature and free kappa and tension."""
				if self.fit_tension: kappa,sigma,vibe = args
				else: (kappa,vibe,),sigma = (args),0.0
				vibe = self.vibe_correct(vibe)
				if self.fit_tension: 
					#---! this was the site of two major errors where we sent vibe
					args_hamiltonian = tuple([np.abs(kappa),sigma]+list(curvatures))
				else: 
					args_hamiltonian = tuple([np.abs(kappa)]+list(curvatures))
				if debug:
					import ipdb;ipdb.set_trace()
				return self.residual(self.binner.bin(oscillator(vibe)[self.band]),
					self.binner.bin(self.hamiltonian(args_hamiltonian)[self.band]))
			return objective_macros
		#---save functions for later
		self.objective = objective
		self.oscillator = oscillator
		self.build_objective_macros = build_objective_macros
		self.optimize_curvature = self._optimize_curvature_legacy
		self.optimize = self._optimize_macros_legacy
	def _optimize_curvature_legacy(self,initial_curvature=0.0):
		"""Run the optimizer over the major parameters. Includes kappa,gamma,vibe,curvature."""
		global stepno
		stepno = 0
		#---! assumes uniform curvatures
		kappa,gamma = self.undulation_fit['kappa'],self.undulation_fit['gamma']
		vibe = self.vibe_init()
		if self.fit_tension: initial_conditions = [kappa,0.0,vibe,initial_curvature]
		else: initial_conditions = [kappa,vibe,initial_curvature]
		fit = scipy.optimize.minimize(self.objective,
			x0=tuple(initial_conditions),callback=callback,method='Nelder-Mead')
		print('\n[OPTIMIZE] solution: %s'%fit)
		self.opt_result = dict(error=fit.fun,kappa=fit.x[0],sigma=fit.x[1] if self.fit_tension else 0.0,
			gamma=0.0,vibe=fit.x[2],curvature=fit.x[3])
		return self.opt_result
	def _optimize_macros_legacy(self,curvature):
		"""Run the optimizer over the major parameters."""
		global stepno
		stepno = 0
		kappa,sigma = self.undulation_fit['kappa'],self.undulation_fit['sigma']
		vibe = self.vibe_init()
		initial_conditions = [kappa,0.0,vibe] if self.fit_tension else [kappa,vibe]
		#---! assume uniform curvatures
		curvatures = [curvature]
		if self.legacy_spec.get('previous',False):
			#---! this functions run much faster. will not be changed until questions are resolved
			#from codes.curvature_coupling.curvature_coupling import prepare_objective,prepare_residual
			#from codes.curvature_coupling.curvature_coupling import prepare_oscillator_function
			#from codes.undulate import blurry_binner
			#---use the previous oscillator
			if self.legacy_spec.get('previous_oscillator',True):
				oscillator_function = prepare_oscillator_function(
					reverse=self.legacy_spec.get('previous_oscillator_reverse',False),positive_vibe=False)
			#----! internal one gives 67kBT and works better for some reason. 74kBT for the previous one.
			else: oscillator_function = self.oscillator_qs
			objective_this = prepare_objective(
				hqs=self.hqs,curvature_fields=self.curvature_fields,
				wavevectors=self.qs,area=self.area,
				curvature_sum_function=curvature_sum_function,fft_function=fft_field,
				band=self.band,blurry_binner=blurry_binner,
				residual_function=self.residual_clip,
				binner_method=self.binner_method,weighting_scheme=None,
				positive_vibe=False,inner_sign=self.inner_sign,
				oscillator_function=oscillator_function,
				ndrops_uniform=8,fix_curvature=curvature)
		else: objective_this = self.build_objective_macros(curvatures=curvatures)
		opt = QuickOpt(objective=objective_this,
			scipy_optimize_function=self.scipy_optimize_function,
			optimize_method=self.scipy_optimize_method,
			silent=False,init=initial_conditions)
		self.opt_result = dict(error=opt.error,kappa=opt.fitted[0],
			sigma=opt.fitted[1] if self.fit_tension else 0.0,
			curvature=curvature,curvatures=curvatures,vibe=opt.fitted[2 if self.fit_tension else 1])
		return self.opt_result

	###
	### MINUS PROTRUSIONS
	def subtract_protrusion(self):
		"""
		Developing a method to remove protrusions and perform the fit.
		"""
		self.undulation_fit = basic_undulation_fitter(
			hqs=self.hqs,qs=self.qs,area=self.area,q_cut=self.q_cut,q_min=self.q_min,
			fit_tension=self.fit_tension,binner_method=self.binner_method,bin_width=self.bin_width)
		self.band = self.undulation_fit['band']
		self.build_hqhq()
		if self.review: self.review_undulation_fit()
		#---in this mode the model for the height fluctuations serves as the reference curve
		#---...for fitting curvatures in the energy space (but we must confirm this theory is sensible)
		self.reference_curve = self.undulation_fit['model'](
			self.qs,self.undulation_fit['sigma'],self.undulation_fit['kappa'],
			self.undulation_fit['gamma'],self.undulation_fit['vibe'])
		#---prepare the hamiltonian
		self.hamiltonian = build_elastic_hamiltonian(
			area=self.area,inner_sign=self.inner_sign,qs=self.qs,hqs=self.hqs,
			#---! call out to the curvature sum function above but should it be passed through?
			curvature_sum_function=curvature_sum_function,curvature_fields=self.curvature_fields,
			curvature_sum_method=self.curvature_sum_method,fft_function=fft_field,
			ndrops_uniform=self.ndrops_uniform,fit_tension=self.fit_tension)
		#---! modify the reference curve
		self.reference_curve_alt = (
			self.undulation_fit['kappa']*self.reference_curve*self.qs**4*self.area)
		#---create the binner
		self.binner = Binner(wavevectors=self.qs[self.band],
			mode=self.binner_method,bin_width=self.bin_width)
		def objective(args):
			"""Optimize curvatures, kappa, and possibly the tension sigma."""
			#---! validation step for incoming arguments?
			if self.fit_tension: 
				kappa,gamma,curvatures = args[0],args[1],args[2:]
				args_hamiltonian = tuple([kappa,gamma]+[curvatures])
			else: 
				kappa,gamma,curvatures = args[0],0.0,args[1:]
				args_hamiltonian = tuple([kappa]+[curvatures])
			return self.residual(self.binner.bin(self.reference_curve_alt[self.band]),
				self.binner.bin(self.hamiltonian(args_hamiltonian)[self.band]))
		def build_objective_macros(curvatures):
			"""Prepare the objective function for a specific curvature."""
			def objective_macros(args):
				"""Objective function with fixed curvature and free kappa and tension."""
				if self.fit_tension: kappa,sigma = args
				else: (kappa,),sigma = (args),0.0
				args_hamiltonian = tuple([np.abs(kappa),np.abs(sigma)]+list(curvatures))
				return self.residual(self.binner.bin(self.reference_curve_alt[self.band]),
					self.binner.bin(self.hamiltonian(args_hamiltonian)[self.band]))
			return objective_macros
		#---save the functions for later
		self.objective = objective
		self.build_objective_macros = build_objective_macros
		self.optimize_curvature = self._optimize_curvature_subtract_protrusion
		self.optimize = self._optimize_macros_subtract_protrusion
	def _optimize_macros_subtract_protrusion(self,curvature):
		"""Run the optimizer over the major parameters."""
		global stepno
		stepno = 0
		kappa,gamma = self.undulation_fit['kappa'],self.undulation_fit['gamma']
		initial_conditions = [kappa,0.0] if self.fit_tension else [kappa]
		#---! assume uniform curvatures
		curvatures = [curvature]
		objective_this = self.build_objective_macros(curvatures=curvatures)
		fit = scipy.optimize.minimize(objective_this,
			x0=tuple(initial_conditions),callback=callback,method='Nelder-Mead')
		print('\n[OPTIMIZE] solution: %s'%fit)
		self.opt_result = dict(error=fit.fun,kappa=fit.x[0],sigma=fit.x[1] if self.fit_tension else 0.0,
			curvature=curvature,curvatures=curvatures)
		return self.opt_result
	def _optimize_curvature_subtract_protrusion(self,initial_curvature=0.0):
		"""Run the optimizer over the major parameters. Includes kappa,gamma,curvature."""
		global stepno
		stepno = 0
		#---! assumes uniform curvatures
		kappa,gamma = self.undulation_fit['kappa'],self.undulation_fit['gamma']
		if self.fit_tension: initial_conditions = [kappa,0.0,initial_curvature]
		else: initial_conditions = [kappa,initial_curvature]
		fit = scipy.optimize.minimize(self.objective,
			x0=tuple(initial_conditions),callback=callback,method='Nelder-Mead')
		print('\n[OPTIMIZE] solution: %s'%fit)
		self.opt_result = dict(error=fit.fun,kappa=fit.x[0],sigma=fit.x[1] if self.fit_tension else 0.0,
			curvature=fit.x[2])
		return self.opt_result

	def review_comprehensive(self,name,**kwargs):
		"""
		Review a single, particular solution to the fit.
		"""
		legacy_vibe = kwargs.pop('legacy_vibe',False)
		solution = kwargs.pop('solution',{})
		legacy_energy = kwargs.pop('legacy_energy',False)
		frozen_weight_backwards = kwargs.pop('frozen_weight_backwards',False)
		extraspectra = kwargs.pop('extraspectra',{})
		height_energy_on_energy = kwargs.pop('height_energy_on_energy',False)
		if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
		axes,fig = square_tiles(2,figsize=(18,10),favor_rows=True)
		#---binner for simplicity
		binner = Binner(mode='perfect',wavevectors=self.qs)
		#---height undulations on the left
		ax = axes[0]
		sigma,kappa,gamma,vibe = (self.undulation_fit['sigma'],self.undulation_fit['kappa'],
			self.undulation_fit['gamma'],self.undulation_fit['vibe'])
		fitted = self.undulation_fit['model'](self.qs,sigma,kappa,gamma,vibe,
			protrusion=self.fit_protrusion,tension=self.fit_tension,correct=self.fit_correction)
		qsort = np.argsort(self.qs)
		#---plot the spectrum
		ax.plot(self.qs[qsort],self.hqhq[qsort],'.',color='m',zorder=1,label='data')
		#---plot the model
		ax.plot(self.qs[qsort],fitted[qsort],'-',color='c',lw=3,zorder=3,label='model')
		ax.axvline(self.q_cut,c='k',lw=1)
		ax.set_xscale('log')
		ax.set_yscale('log')
		#---plot heights in the inset on the energy scale for reference
		axins = inset_axes(ax,width="35%",height="35%",loc=1)
		heights_on_energy = self.hqhq*self.qs**4*kappa*self.area
		xs,ys = binner.bin(self.qs),binner.bin(heights_on_energy) 
		xsort = np.argsort(xs)
		axins.plot(xs[xsort],ys[xsort],'-',color='m',lw=2,zorder=3,label='model',alpha=0.5)
		axins.plot(self.qs[qsort],fitted[qsort]*self.qs[qsort]**4*kappa*self.area,
			'-',color='c',lw=3,zorder=3,label='model')
		axins.axvline(self.q_cut,c='k',lw=1)
		axins.axhline(1.0,c='k',lw=1)
		axins.set_xscale('log')
		axins.set_yscale('log')
		#---energy spectrum on the right
		ax = axes[1]
		if height_energy_on_energy:
			ax.plot(self.qs[qsort],fitted[qsort]*self.qs[qsort]**4*kappa*self.area,
				'-',color='c',lw=3,zorder=3,label='model')
		kappa = self.undulation_fit['kappa']
		#---plot the elastic hamiltonian on the energy scale
		spectrum = self.hamiltonian((self.opt_result['kappa'],ic.opt_result['sigma'],
			#---! critical error here where you sent vibe instead of curvature and it wasted a ton of time
			ic.opt_result['curvature']))
		xs,ys = binner.bin(self.qs),binner.bin(spectrum)
		xsort = np.argsort(xs)
		ax.plot(xs,ys,'-',c='k',lw=1,alpha=0.6,zorder=1)
		ax.axvline(self.q_cut,c='k',lw=1)
		ax.axhline(1.0,c='k',lw=1)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.set_ylim((0.1,10.0))
		#---plot the vibrational correction
		if legacy_vibe:
			vibe = self.oscillator(vibe=legacy_vibe)
			if frozen_weight_backwards: ys = log_reverse(vibe[qsort])
			else: ys = vibe[qsort]
			ax.plot(self.qs[qsort],ys,'-',color='r',lw=2,zorder=4,label='oscillator')
		if legacy_energy:
			if self.fit_tension: 
				args_hamiltonian = tuple([solution['kappa'],solution['sigma']]+[solution['curvature']])
			else: args_hamiltonian = tuple([solution['kappa'],]+[solution['curvature']])
			energies = self.hamiltonian(args_hamiltonian)
			frozen = self.oscillator(vibe=legacy_vibe)
			if frozen_weight_backwards: xs,ys = binner.bin(self.qs),binner.bin(energies*frozen)
			else: xs,ys = binner.bin(self.qs),binner.bin(energies/frozen)
			ax.plot(xs,ys,lw=2,zorder=10,color='cyan',solid_capstyle='round',
				path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])
		#---additional energy spectra
		for ename,spec in extraspectra.items():
			kwargs = dict(color=spec.get('color','k'))
			if 'zorder' in spec: kwargs.update(zorder=spec['zorder'])
			if spec.get('outline',True): kwargs.update(
				path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])
			qsort = np.argsort(spec['qs'])
			ax.plot(spec['qs'][qsort],spec['energy'][qsort],**kwargs)
		picturesave('fig.%s.%s'%(name,self.sn),
			work.plotdir,backup=False,version=False,meta={},extras=[])

def refresh_memory(**kwargs):
	"""
	"""
	global midplane_method,memory
	#---refreh memory if the midplane_method is different
	if midplane_method!=kwargs['midplane_method'] or 'memory' not in globals() or memory==None:
		#---call to the external transformer
		loader_spec = {'module':'codes.curvature_coupling_loader',
			'function':'curvature_coupling_loader_membrane'}
		loader_func = gopher(loader_spec,module_name='module',variable_name='function',work=work)
		#---load and transform the heights
		memory = loader_func(midplane_method=midplane_method,
			data=dict(undulations=data['undulations']))
		midplane_method = kwargs['midplane_method']

def prepare_curvatures(curvatures_sweep):
	"""
	Prepare a list of curvatures.
	"""
	if curvatures_sweep=='extreme':
		halfsweep = np.concatenate((np.arange(0.01,0.1+0.01,0.01),np.arange(0.2,1.0+0.1,0.1),
			np.arange(2.,10.+1.,1.)))
		curvatures = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))
	elif curvatures_sweep=='small':
		halfsweep = np.concatenate((
			np.arange(0.001,0.01+0.001,0.001),
			np.arange(0.02,0.1+0.01,0.01)))
		curvatures = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))
	elif curvatures_sweep=='reasonable':
		halfsweep = np.concatenate((
			np.arange(0.001,0.01+0.001,0.001),
			np.arange(0.02,0.1+0.01,0.01),
			np.arange(0.2,1.0+0.2,0.2)))
		curvatures = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))
	elif curvatures_sweep=='legacy':
		curvatures = np.array([0.0,0.005,0.01,0.014,0.018,0.02,0.024,0.028,0.032,0.04,0.05])
	else: raise Exception('unclear curvatures_sweep %s'%curvatures_sweep)
	return curvatures

def survey_error_landscape_singleton(name,**kwargs):
	"""
	Brute force search of the error landscape.
	"""
	q_cut = kwargs['q_cut']
	ndrops_uniform = kwargs['ndrops_uniform']
	frameslice = kwargs['frameslice']
	extent_keys = kwargs['extent_keys']
	binner_method = kwargs['binner_method']
	mode = kwargs['mode']
	curvatures_sweep = kwargs['curvatures_sweep']
	fit_tension = kwargs['fit_tension']
	sn = kwargs['sn']
	legacy_spec = kwargs.pop('legacy_spec',{})
	inner_sign = kwargs.pop('inner_sign',1.0)

	curvatures = prepare_curvatures(curvatures_sweep)
	refresh_memory(midplane_method=kwargs['midplane_method'])

	global job_toc
	job_toc = {}
	extents = np.array([datas[k][sn]['spec']['extents']['extent'] for k in extent_keys])

	start = time.time()
	#---prepare objective functions for each extent
	#---! cannot pickle function objects so cannot do this in parallel
	for ii,(design_name,extent) in enumerate(zip(extent_keys,extents)):
		hqs = memory[(sn,'hqs')]
		curvature_fields = datas[design_name][sn]['fields_unity']
		ic = self = InvestigateCurvatureMicro(sn=sn,hqs=hqs,name=name,frameslice=frameslice,
			curvature_fields=curvature_fields,fit_tension=fit_tension,ndrops_uniform=ndrops_uniform,
			mode=mode,midplane_method=midplane_method,q_cut=q_cut,legacy_spec=legacy_spec)
		kappa,sigma = ic.undulation_fit['kappa'],ic.undulation_fit['sigma']
		for jj,curvature in enumerate(curvatures):
			print('\n')
			status('preparing objective function',i=ii*len(curvatures)+jj,
				looplen=len(extents)*len(curvatures),start=start,tag='build')
			error = ic.optimize(curvature=curvature)['error']
			job_toc[(extent,curvature)] = error

	###---PLOT ERROR LANDSCAPES
	figname = 'fig.coupling.error_landscape.%s.%s'%(sn,name)
	figsize = (12,10)
	contour_interp_pts = 100
	contour_line_skip = 4
	contour_nlevels = 100
	under_color = 'm'
	
	axes,fig = square_tiles(2,figsize,favor_rows=True)
	ax = axes[0]
	raw = np.array([[job_toc[(e,c)] for e in extents] for c in curvatures])
	#---literal
	kwargs = dict(extent=[min(curvatures),max(curvatures),
		min(extents),max(extents)],aspect=(curvatures.ptp()/extents.ptp()))
	#---figurative
	kwargs = dict(extent=[0,len(curvatures),0,len(extents)],aspect=(float(len(curvatures))/len(extents)))
	ax.imshow(raw.T,origin='lower',interpolation='nearest',**kwargs)
	ax.set_xticks(np.arange(len(curvatures))+0.5)
	ax.set_yticks(np.arange(len(extents))+0.5)
	ax.set_xticklabels(['%.3f'%i for i in curvatures],rotation=90)
	ax.set_yticklabels(['%.1f'%i for i in extents])
	#---contour
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
	errormap = scipy.interpolate.griddata(curvature_extent_error[:,:2],curvature_extent_error[:,2],
		(grid_x,grid_y),method='cubic')
	levels = np.linspace(error_min,error_max,contour_nlevels)
	cs = ax.contourf(grid_x,grid_y,errormap,levels=levels,vmax=error_max,vmin=error_min,
		extend='both',origin='lower',lw=2,zorder=3,cmap=mpl.cm.jet)
	cs.cmap.set_over('w')
	if under_color: cs.cmap.set_under(under_color)
	levels_contour = levels[np.where(levels<=contour_line_max)][::contour_line_skip]
	if False: cs_lines = ax.contour(grid_x,grid_y,errormap,vmax=error_max,
		vmin=error_min,levels=levels_contour,
		extend='both',origin='lower',linewidths=0.5,colors='k',zorder=4)
	ax.set_aspect(curvatures.ptp()/extents.ptp())
	picturesave(figname,work.plotdir,backup=False,version=False,meta={})

def survey_error_landscapes(trials):
	"""
	Loop over trials to generate error landscapes.
	"""
	#---midplane method is handled with a global cursor
	global midplane_method,memory
	for tnum,trial in enumerate(trials):
		midplane_method = trial['midplane_method']
		trial['inner_sign_note'] = 'p' if trial.get('inner_sign',1)>0 else 'n'
		name = ('REVIEW.%(mode)s.mm.%(midplane_method)s.bin.%(binner_method)s.'+
			'qcut.%(q_cut)s.csweep.%(curvatures_sweep)s.i%(inner_sign_note)s')%trial
		status('generating %s'%name,tag='SWEEP')
		survey_error_landscape_singleton(name=name,**trial)

#---iterative development
if __name__=='__replotting__':

	#---direction
	devmode = ['survey','survey_single','live'][-1]
	#---protect memory
	if len(work.sns())>1: raise Exception

	#---choose simulation
	mode = ['legacy','subtract_protrusion'][0]
	sn = ['membrane-v1005','membrane-v650-enthx4-dev',
		'membrane-v651-enthx8',][-1]
	legacy_spec_name = ['legacy_original','legacy_original_reverse',
		'legacy_standard','legacy_unstuck'][-1]
	frameslice_name = ['limited','everything'][0]
	extents_keys_name = ['v1_high','v1','v2_all'][1]
	curvatures_sweep = ['extreme','small','reasonable','legacy'][-1]
	legacy_spec = {
		'legacy_original':{
			'grease':True,'safe':'negative','previous':True,'previous_oscillator':True,
			'previous_oscillator_reverse':False},
		'legacy_original_reverse':{
			'grease':True,'safe':'negative','previous':True,'previous_oscillator':True,
			'previous_oscillator_reverse':True},
		'legacy_standard':{
			'grease':True,'safe':True,'previous':True,'previous_oscillator':False,
			'previous_oscillator_reverse':True},
		# slow but it works
		'legacy_unstuck':{
			'grease':False,'safe':True,'previous':False,'previous_oscillator':False,
			'previous_oscillator_reverse':True,'reverse_oscillator':False},
		}[legacy_spec_name]
	frameslice = {'everything':None,'limited':100}[frameslice_name]
	extent_keys = {
		'v1':['v10_fixed_extent_0.25','v3_fixed_extent_0.5',
			'v1_fixed_extent_1','v2_fixed_extent_2','v5_fixed_extent_3','v4_fixed_extent_4'][:6],
		'v1_high':['v10_fixed_extent_0.25','v3_fixed_extent_0.5',
			'v1_fixed_extent_1','v2_fixed_extent_2',
			'v5_fixed_extent_3','v4_fixed_extent_4','v6_fixed_extent_5','v7_fixed_extent_6',
			'v8_fixed_extent_8','v9_fixed_extent_10'][:],
		'v2_all':['v3_fixed_extent_all_frames_0.5','v1_fixed_extent_all_frames_1',
			'v2_fixed_extent_all_frames_2','v5_fixed_extent_all_frames_3',
			'v4_fixed_extent_all_frames_4'],}[extents_keys_name]

	#---default parameters
	params = dict(sn=sn,mode=mode,fit_tension=True,inner_sign=-1.0,
		ndrops_uniform=work.meta[sn]['nprots'],curvatures_sweep=curvatures_sweep,
		legacy_spec=legacy_spec,extent_keys=extent_keys,frameslice=frameslice)

	trials = hypothesizer(*(
		{'route':['midplane_method'],'values':['flat','average'][:1]},
		{'route':['q_cut'],'values':[1.0,2.0,4.0][:1]},
		{'route':['binner_method'],'values':['explicit','perfect','blurry'][:1]},),
		default=params)

	#---make the survey
	if devmode=='survey': survey_error_landscapes(trials)
	#---bleeding edge debugging
	elif devmode=='live': 

		#---directions
		legacy_compare = 0 #---include a spectra from the original calculation
		always_redo = 1 #---remake the ic each time during development
		use_precomp = 0 #---use previous optimization hardcoded below
		free_curvature = 0 #---fit with free or fixed curvature
		render = 1 #---make the figure
	
		#---fold these above
		midplane_method = 'flat'
		frameslice = 100
		frameslice = None
		fit_tension = True
		q_cut = 1.0
		inner_sign = -1
		frozen_weight_backwards = True
		binner_method = ['explicit','perfect','blurry'][0]
		ndrops_uniform = work.meta[sn]['nprots']
		refresh_memory(midplane_method=midplane_method)
		hqs = memory[(sn,'hqs')]
		#---! extent is irrelevant for now because we are just investigating the no-curvature case
		design_name = ['v10_fixed_extent_0.25',
			'v3_fixed_extent_0.5','v2_fixed_extent_2'][-1]
		curvature_fields = datas[design_name][sn]['fields_unity']
		#---name for the figure
		name = ['legacy_compare_original','legacy_compare','subtract_protrusion',
			'legacy_compare.unstuck'][-1]
		#---pull previous spectra for comparison
		if 'original' not in globals() and legacy_compare:
			import pickle
			original = pickle.load(open('/home/rpb/spectra.pkl'))
			binner = Binner(mode='perfect',wavevectors=original['wavevectors'])
			extraspectra = {'original':
				{'qs':binner.bin(original['wavevectors']),'energy':binner.bin(original['best_energy']),
				'color':'#db4d76','zorder':10,'outline':True}}
		elif not legacy_compare: extraspectra = {}
		if 'ic' not in globals() or always_redo:
			ic = InvestigateCurvatureMicro(sn=sn,hqs=hqs,name=name,frameslice=frameslice,
				curvature_fields=curvature_fields,fit_tension=fit_tension,ndrops_uniform=ndrops_uniform,
				mode=mode,midplane_method=midplane_method,q_cut=q_cut,inner_sign=inner_sign,
				binner_method=binner_method,legacy_spec=legacy_spec)
			if free_curvature:
				if not use_precomp: ic.optimize_curvature()
				else: raise Exception
			else:
				if not use_precomp: opt_result = ic.optimize(0.0)
				else: raise Exception
		if render: 
			if mode=='subtract_protrusion':
				es = dict(model=dict(energy=ic.hamiltonian((ic.opt_result['kappa'],
					ic.opt_result['sigma'],0.0)),qs=ic.qs,color='r',outline=True))
				ic.review_comprehensive(name,legacy_vibe=ic.opt_result.get('vibe',False),solution=ic.opt_result,legacy_energy=False,frozen_weight_backwards=frozen_weight_backwards,extraspectra=es,
					height_energy_on_energy=True)
			elif mode=='legacy':
				ic.review_comprehensive(name,legacy_vibe=ic.opt_result.get('vibe',False),solution=ic.opt_result,legacy_energy=True,frozen_weight_backwards=frozen_weight_backwards,
					extraspectra=extraspectra)

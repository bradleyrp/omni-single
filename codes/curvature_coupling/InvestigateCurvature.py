#!/usr/bin/env python

"""
CURVATURE-UNDULATION COUPLING (version 3)
!!! Currently abandoned for version 3b!!!
"""

import copy,importlib,time,os,re
from codes.hypothesizer import hypothesizer
from base.tools import status
from omnicalc import store,load
import numpy as np
import scipy
import scipy.optimize
import multiprocessing as mp
machine_eps = eps = np.finfo(float).eps
import codes.curvature_coupling.tools as cctools
from codes.curvature_coupling.tools import manyjob,manyjob_worker
from codes.curvature_coupling.tools import construct_curvature_fields_trajectory
global Nfeval

###---UTILITIES

def make_fields(**kwargs):
	return cctools.construct_curvature_field(**kwargs)
def framelooper(total,start=None,text='frame'):
	"""
	When performing parallel calculations with joblib we pass a generator to count the number of 
	tasks and report the time.
	"""
	for fr in range(total):
		status(text,i=fr,looplen=total,tag='parallel',start=start)
		yield fr
from joblib import Parallel,delayed
from joblib.pool import has_shareable_memory
def basic_compute_loop(compute_function,looper,run_parallel=True,debug=False):
	"""
	Canonical form of the basic compute loop.
	"""
	start = time.time()
	if run_parallel:
		incoming = Parallel(n_jobs=8,verbose=10 if debug else 0)(
			delayed(compute_function,has_shareable_memory)(**looper[ll]) 
			for ll in framelooper(len(looper),start=start))
	else: 
		incoming = []
		for ll in framelooper(len(looper)):
			incoming.append(compute_function(**looper[ll]))
	return incoming

dotplace = lambda n : re.compile(r'(\d)0+$').sub(r'\1',"%3.5f"%float(n)).rjust(8)

class SpecialClass(object):
	"""
	Decorator class for mapping our "database" specification onto SQLAlchemy.
	"""
	def __init__(self,**kwargs):
		self.test_spec = kwargs['test_spec']
		self.change_spec = kwargs['change_spec']
		self.order = kwargs['order']
		self.table = kwargs['table']
		self.table_name = kwargs['table_name']
		self.anticipate = kwargs['anticipate']
	def __call__(self,cls):
		class Wrapped(cls):
			test_spec = self.test_spec
			change_spec = self.change_spec
			order = self.order
			anticipate = self.anticipate
			__table__ = self.table
			__tablename__ = self.table_name
			def __init__(self,**kwargs):
				for key in self.order: self.__dict__[key] = None
				for key in self.anticipate: self.__dict__[key] = -1
				for key,val in kwargs.items(): self.__dict__[key] = val
				#---test and validate
				for tt,test in enumerate(self.test_spec): 
					if not test(self.__dict__):
						raise Exception('failed database test number %d (from zero)'%(tt))
				for change in self.change_spec: self.__dict__.update(**change(self.__dict__))
			def base(self):
				return {key:self.__dict__[key] for key in self.order if key not in self.anticipate}
			def dict(self):
				return {key:self.__dict__[key] for key in self.order}
			def present(self):
				#---! recent modification here !!!!!!! 2017.07.19
				return {key:self.__dict__[key] for key in self.order if self.__dict__[key]!=None}
		return Wrapped

if False:
	def callback(args):
		"""
		Watch the optimization.
		"""
		global Nfeval,name_groups,objective
		text = ' step = %d '%Nfeval+' '.join([name+' = '+dotplace(val)
			for name,val in zip(name_groups,args)+[('error',objective(args))]])
		status('searching! '+text,tag='optimize')
		Nfeval += 1

###---CONTROLLER

class InvestigateCurvature:

	"""
	Manage the entire curvature coupling calculation.
	"""

	def __init__(self,**kwargs):
		"""Assemble the necessary modules and execute."""
		self.plotloader = kwargs.pop('plotloader',plotload)
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		#---get the data according to instructions in the metadata
		self.fetch_data_and_specs()
		self.memory = self.loader_func(self.data)
		self.attempts = work.plots.get('curvature',{}).get('specs',{}).get('attempts',{})
		#---compute and plot modes
		self.mode = kwargs.pop('mode','plot')
		if self.mode=='plot': 
			curvature_specs = work.plots.get('curvature',{}).get('specs',{})
			plot_cursor = curvature_specs['plot_cursor']
			plot_specs = curvature_specs['plot_specs'].get(plot_cursor,
				curvature_specs['attempts'][plot_cursor])
			#!!! check if dat is here for opt...
			self.signifier = plot_cursor
			self.attempt(signifier=plot_cursor,**plot_specs)
		elif self.mode=='compute': self.main()
		else: raise Exception('unclear mode %s'%self.mode)
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)

	def main(self):
		"""..."""
		#---master loop over the attempts
		for signifier,attempt in self.attempts.items():
			self.signifier = signifier
			self.attempt(signifier,**attempt)

	def attempt(self,signifier,**attempt):
		"""..."""
		#---cursors
		self.sig,self.spec = signifier,attempt
		#---select the hypotheses
		self.style = attempt.get('style',None)
		if self.style=='manual': 
			self.prepare_rootdir('curvature_%s'%signifier)
			self.prepare_database()
			self.manual()
		elif self.style=='gridded': self.gridded()
		elif self.style=='wilderness': self.wilderness()
		else: raise Exception('unclear attempt style: %s'%self.style)

	def gopher(self,spec,module_name,variable_name):
		"""..."""
		mod = importlib.import_module(spec[module_name])
		target = mod.__dict__.get(spec[variable_name],None)
		#---the database design might need work so we always export it
		global work
		mod.work = work
		if not target: raise Exception('add %s and %s to the specs'%(module_name,variable_name))
		return target

	def fetch_data_and_specs(self):
		"""..."""
		#---load the upstream data
		self.data,self.calc = self.plotloader(plotname)
		self.meta_spec = work.plots.get(plotname,{}).get('specs',{})
		#---load the module responsible for packing the data
		loader_spec = self.meta_spec.get('loader',None)
		if not loader_spec: raise Exception('you must add specs,loader to the plot entry for curvature')
		self.loader_func = self.gopher(loader_spec,module_name='module',variable_name='function')
		#---load the module responsible for packing the data
		database_spec = self.meta_spec.get('design',None)
		if not database_spec: raise Exception('you must add specs,design to the plot entry for curvature')
		self.database = self.gopher(database_spec,module_name='module',variable_name='def')

	def prepare_rootdir(self,dropname):
		"""This calculation is off-pathway so we make a folder in the post directory for it."""
		#---root directory
		self.rootdir = os.path.join(work.paths['post_data_spot'],dropname,'')
		if not os.path.isdir(self.rootdir): 
			os.mkdir(self.rootdir)
			#---make subdirectories for hypotheses and curvature fields
			for name,sub_dn in [('rootdir_cc','hypotheses'),('rootdir_cf','curvature_fields')]:
				os.mkdir(os.path.join(work.paths['post_data_spot'],dropname,sub_dn))
		#### else: raise Exception('refusing to write into preexisting directory: %s'%self.rootdir)
		else: status('data already exists',tag='note')
		for name,sub_dn in [('rootdir_cc','hypotheses'),('rootdir_cf','curvature_fields')]:
			self.__dict__[name] = os.path.join(work.paths['post_data_spot'],dropname,sub_dn)
		#---name the data files
		self.namer_cf = lambda pk : os.path.join(self.rootdir_cf,'curvature_field.%d.dat'%pk)
		self.namer_cc = lambda pk : os.path.join(self.rootdir_cc,'hypothesis.%d.dat'%pk)

	def magic_data(self,database,engine,Session):
		"""
		Automatically generate declaratives and mappings for SQL Alchemy to handle the data with restrictions
		defined by the SpecialClass above.
		"""
		from sqlalchemy import MetaData
		from sqlalchemy.ext.declarative import declarative_base
		from sqlalchemy import Table,Column,Float,String,Integer
		type_convert = {float:Float,str:String,int:Integer}
		outgoing,self.magic_tables = [],[]
		metadata = MetaData()
		Base = declarative_base()
		for dataname,dataspec in database.items():
			basename = 'table_%s'%dataspec['table_name']
			if basename in self.__dict__.keys(): raise Exception('%s in self already'%basename)
			self.__dict__[basename] = Table(
				dataspec['table_name'],Base.metadata,
				*[Column('id',Integer,primary_key=True)]+[Column(key,type_convert[val]) 
				for key,val in dataspec['data_map']])
		for dataname,dataspec in database.items():
			class_name = dataspec['table_name'].capitalize()
			#---mimic decorator "@SpecialClass(**kwargs)"" over "class Hypothesis(object): pass"
			if class_name+'_inner' in self.__dict__.keys(): 
				raise Exception('%s in self already'%(class_name+'_inner'))
			#---make note of the tables we have created
			self.magic_tables.append(class_name)
			self.__dict__[class_name+'_inner'] = SpecialClass(
				test_spec=dataspec['test_spec'],
				change_spec=dataspec['change_spec'],
				order=zip(*dataspec['data_map'])[0],
				table=Base.metadata.tables[dataspec['table_name']],
				table_name=dataspec['table_name'],
				anticipate=dataspec['anticipate'],
				)(object)
			#---note that we replaced globals with self here to contain this hack
			exec('class %s(self.__dict__[\'%s_inner\'],Base): pass'%(class_name,class_name))
			outgoing.append(eval(class_name))
		Base.metadata.create_all(engine)
		session = Session()
		session.commit()
		return session,outgoing

	def prepare_database(self):
		"""..."""
		#---database interface (now separated by table)
		self.sessions = {}
		#---create engines and sessions
		from sqlalchemy import Table,Column,Float,String,Integer,create_engine
		from sqlalchemy.orm import sessionmaker,mapper
		engine_cc = create_engine('sqlite:///%s/hypotheses.sqlite'%self.rootdir_cc)
		Session_cc = sessionmaker(bind=engine_cc)
		engine_cf = create_engine('sqlite:///%s/curvature_fields.sqlite'%self.rootdir_cf)
		Session_cf = sessionmaker(bind=engine_cf)
		#---attach engines and sessions
		self.database['hypothesis']['engine'] = engine_cc
		self.database['hypothesis']['session_maker'] = Session_cc
		self.database['field']['engine'] = engine_cf
		self.database['field']['session_maker'] = Session_cf
		#---???
		for table_name in self.database:
			self.sessions[table_name],outgoing = self.magic_data({table_name:self.database[table_name]},
				self.database[table_name]['engine'],self.database[table_name]['session_maker'])
			#---note that we replaced globals with self here to contain this hack
			for key in outgoing: self.__dict__[key.__name__] = key
		self.session_makers = dict([(key,self.database[key]['session_maker']) 
			for key in ['hypothesis','field']])

	def manual_prepare_compute(self):
		"""
		"""
		#---load the database
		start = time.time()
		session = self.sessions['hypothesis']
		for hh,hypo in enumerate(self.hypotheses):
			status('populating hypos',tag='load',i=hh,looplen=len(self.hypotheses),start=start)
			#---reduce step before checking database
			hypo_full = self.Hypothesis(**hypo)
			matches = session.query(self.Hypothesis).filter_by(**hypo_full.base()).all()
			if not any(matches): session.add(hypo_full)
		session.commit()
		session = self.sessions['field']
		for hh,hypo in enumerate(self.hypotheses):
			status('populating fields',tag='load',i=hh,looplen=len(self.hypotheses),start=start)
			#---reduce step before checking database
			hypo_full = self.Field(**hypo)
			matches = session.query(self.Field).filter_by(**hypo_full.base()).all()
			if not any(matches): 
				session.add(hypo_full)
		session.commit()
		#---integrity checks on database rows
		hypotheses_reduced = [i.dict() for i in self.sessions['hypothesis'].query(self.Hypothesis).all()]
		fields_reduced = [i.dict() for i in self.sessions['field'].query(self.Field).all()]
		assert not [i for i in hypotheses_reduced if i['mapping']=='protein' and i['curvature']==0.0]
		assert not [i for i in hypotheses_reduced if i['curvature']==0.0 
			and not (i['sigma_a']==1.0 and i['isotropy']==1.0 and i['sigma_b']==1.0)]

	def manual_populate_fields(self):
		"""
		"""
		cctools.data = self.data
		#---compute pending fields according to populated rows
		fns = [(i.id,self.namer_cf(i.id)) for i in self.sessions['field'].query(self.Field).all()]
		pending = [(pk,fn) for pk,fn in fns if not os.path.isfile(fn)]
		if pending:
			#---loop over absent files
			start = time.time()
			for ii,(pk,fn) in enumerate(pending):
				status('computing curvature field',tag='compute',i=ii,looplen=len(pending),start=start)
				hypo = self.sessions['field'].query(self.Field).filter_by(id=pk).one().dict()
				sn = hypo['sn']
				dat = self.data['undulations'][sn]['data']
				vecs = dat['vecs']
				mn = np.shape(dat['mesh'])[2:]
				fields = construct_curvature_fields_trajectory(vecs=vecs,mn=mn,**hypo)
				store({'fields':np.array(fields['fields'])},os.path.basename(fn),self.rootdir_cf,
					attrs={key:val for key,val in fields.items()+hypo.items() 
					if key!='fields'},verbose=False)

	def manual_evaluate_hypotheses(self):
		"""
		"""
		#---manual execution requires export of the data tables to the tools
		#---! prefer this to be systematic, but exporting is already offbeat
		cctools.namer_cf = self.namer_cf
		cctools.namer_cc = self.namer_cc
		cctools.Field = self.Field
		cctools.Hypothesis = self.Hypothesis
		cctools.memory = self.memory
		cctools.rootdir_cf = self.rootdir_cf
		cctools.rootdir_cc = self.rootdir_cc
		#---solve the hypotheses
		#---for memory efficiency we queue up hypotheses according to which curvature field they require
		#---note that we had a simpler, memory-hogging loop in a previous iteration of this code
		fns = [(i.id,self.namer_cc(i.id)) for i in self.sessions['hypothesis'].query(self.Hypothesis).all()]
		pending = [(pk,fn) for pk,fn in fns if not os.path.isfile(fn)]
		if pending:
			self.hypotheses = [self.sessions['hypothesis'].query(
				self.Hypothesis).filter_by(id=pk).one().dict() 
				for pk in zip(*pending)[0]]
			fields_required = [self.sessions['field'].query(self.Field).filter_by(**f.dict()).one() 
				for f in [self.Field(**h) for h in self.hypotheses]]
			field_ids_by_hypothesis = np.array([f.id for f in fields_required])
			unique_field_ids = np.unique(field_ids_by_hypothesis)
			#---compute the curvatures in batches
			for uu,ufid in enumerate(unique_field_ids):
				status('computing all hypotheses for field %d/%d'%(uu,len(unique_field_ids)),tag='compute')
				hypo_subset = [self.hypotheses[j] for j in np.where(field_ids_by_hypothesis==ufid)[0]]
				key_cf = ('curvature',ufid)
				self.memory[key_cf] = load(os.path.basename(self.namer_cf(ufid)),
					cwd=os.path.dirname(self.namer_cf(ufid)))
				#---queue for each part of the computation
				queue_hypothesis = mp.Queue()
				#---solve
				manyjob(single=False,
					function=manyjob_worker,
					queue=queue_hypothesis,
					session_classes=self.session_makers,
					objects=hypo_subset,
					kwargs={'preloaded_curvature':True})
				#---clear that hypothesis from memory
				del self.memory[key_cf]
			status('done all batches',tag='compute')

	def manual(self):
		"""
		Original "manual" optimization method in which we sweep over many hypotheses and then later 
		choose the best one.
		"""
		#---unpack the parameters from the attempt spec
		sweep_sigma_a = self.spec.get('sigma_a',[1.0])
		sweep_c0 = self.spec.get('C_0',[0.1])
		sweep_isotropy = self.spec.get('isotropy',[1.0])
		sweep_theta = np.array(self.spec.get('theta',[0.0])).astype(float)
		#---theta in the YAML is in degrees but we convert to radians here
		sweep_theta = sweep_theta/180.0*np.pi
		motion = self.spec.get('motion')
		mapping = self.spec.get('mapping')

		global work
		sns = work.sns()
		sns_unfree = [s for s in sns if work.meta[s].get('nprots',1)>0]
		#---currently we store a timecode in the database however it only depends on the name in the meta
		#---! this needs improved
		timecode = self.calc['calcs']['undulations']['specs']['slice_name']

		#---prepare hypotheses independent of the simulations
		hypotheses_base = hypothesizer(
			dict(route=['curvature'],values=sweep_c0),
			dict(route=['sigma_a'],values=sweep_sigma_a),
			dict(route=['isotropy'],values=sweep_isotropy),
			dict(route=['theta'],values=sweep_theta),
			default={'mapping':mapping,'fallback':'none','motion':motion,'theta':0.0})
		hypotheses = []
		#---cross the hypotheses with 
		for hypo in hypotheses_base:
			for sn in sns:
				if sn not in sns_unfree:
					for sn_alt in sns_unfree:
						hypo_new = dict(hypo)
						hypo_new.update(sn=sn,fallback=sn_alt,timecode=timecode)
						hypotheses.append(hypo_new)
				else:
					hypo_new = dict(hypo)
					hypo_new.update(sn=sn)
					hypo_new.update(sn=sn,timecode=timecode)
					hypotheses.append(hypo_new)
		self.hypotheses = hypotheses
		self.manual_prepare_compute()
		self.manual_populate_fields()
		self.manual_evaluate_hypotheses()

	###---PLOTTING UTILITIES

	def generate_standard_manual_landscapes(self,**kwargs):
		"""Generate the error landscape over curvature and extent."""
		global work
		sns = work.sns()
		hypos = self.hypotheses
		#---first we filter the hypotheses by kwargs
		hypos_sub = [h for h in hypos if all([h[key]==val for key,val in kwargs.items()])]
		#---get the sorted extents and curvatures
		curvatures = np.unique([i['curvature'] for i in hypos_sub])
		extents = np.unique([i['sigma_a'] for i in hypos_sub])
		curvatures_inds,extents_inds = [dict([(i,ii) 
			for ii,i in enumerate(v)]) for v in [curvatures,extents]]
		#---we repackage the hypotheses by simulation then curvatures by extent
		landscapes = dict([(sn,0*np.ones((len(curvatures),len(extents)))) for sn in sns])
		#---ensure that we have one result per hypothesis otherwise we need to whittle things down
		#---perform this check by sweeping over simulations, curvatures, extents checking for redundancy
		for sn in sns:
			for c0 in curvatures:
				for sigma_a in extents:
					candidates = [h for h in hypos_sub if h['sn']==sn and 
						h['curvature']==c0 and h['sigma_a']==sigma_a]
					if len(candidates)!=1: raise Exception('non-unique result: %s'%candidates)
		#---collect the results
		for hypo in hypos_sub:
			#---each hypothesis has a single row, which is why we perform the check above
			row = self.sessions['hypothesis'].query(self.Hypothesis).filter_by(
				**self.Hypothesis(**hypo).base()).one()
			landscapes[row.sn][curvatures_inds[hypo['curvature']],extents_inds[hypo['sigma_a']]] = row.error
		return landscapes

	###---FIRST (FAILED!) GRIDDER ATTEMPT

	def master_decorator(self,sn,extent,lenscale=1.0):
		"""
		Given the wavevectors, Fourier-transformed heights, and core curvature fields, we return
		the master objective function.
		"""
		global nprots
		#---prepare
		hqs = self.memory[(sn,'hqs')]
		vecs = self.memory[(sn,'vecs')]
		vecs_mean = np.mean(vecs,axis=0)
		#---! HACKING and hardcoding!
		npts = (vecs_mean/self.spec['spacer']).round()[:2].astype(int)
		np.meshgrid(*[np.linspace(0,i,npts[ii]) for ii,i in enumerate(vecs_mean[:2])])
		centers = np.concatenate(np.transpose(np.meshgrid(*[np.linspace(0,i,npts[ii])[1:-1] 
			for ii,i in enumerate(vecs_mean[:2])])))
		nprots = len(centers)
		#---! previously: centers = np.mean(data_prot[sn]['data']['trajectory_com'],axis=0)/vecs_mean
		kwargs = dict(curvature=1.0,mn=hqs.shape[1:],theta=0.0,sigma_a=extent,sigma_b=extent)
		core_cfs = np.array([cctools.construct_curvature_field(vecs=vecs_mean,
			centers=[c],**kwargs) for c in centers])
		nframes = len(vecs)
		m,n = mn = np.shape(hqs)[1:]
		Lx,Ly = np.mean(vecs,axis=0)[:2]
		q2d = lenscale*np.array([[np.sqrt(
			((i-m*(i>m/2))/((Lx)/1.)*2*np.pi)**2+
			((j-n*(j>n/2))/((Ly)/1.)*2*np.pi)**2)
			for j in range(0,n)] for i in range(0,m)])
		q_raw = np.reshape(q2d,-1)[1:]
		area = (Lx*Ly/lenscale**2)

		#---! remove this eventually
		tweak = {}
		signterm = tweak.get('inner_sign',-1.0)
		lowcut = kwargs.get('lowcut',tweak.get('lowcut',0.0))
		band = cctools.filter_wavevectors(q_raw,low=lowcut,high=tweak.get('hicut',1.0))
		residual_form = kwargs.get('residual_form',tweak.get('residual_form','log'))
		if residual_form == 'log':
			def residual(values): return sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
		elif residual_form == 'linear': 
			def residual(values): return sum((values-1.0)**2)/float(len(values))
		else: raise

		def multipliers(x,y): return x*np.conjugate(y)

		def master(args):
			"""
			Fit parameters are defined in sequence for the optimizer.
			They are: kappa,gamma,vibe,*curvatures-per-dimple.
			"""
			kappa,gamma,vibe = args[:3]
			curvatures = args[3:]
			composite = np.sum(core_cfs.T*np.array(curvatures),axis=2)
			cqs = np.tile(cctools.fft_field(np.array([composite])),(len(hqs),1,1))
			termlist = [multipliers(x,y) for x,y in [(hqs,hqs),(hqs,cqs),(cqs,hqs),(cqs,cqs)]]
			termlist = [np.reshape(np.mean(k,axis=0),-1)[1:] for k in termlist]
			#---skipping assertion and dropping imaginary
			termlist = [np.real(k) for k in termlist]
			hel = (kappa*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
				+signterm*termlist[2]*q_raw**2+termlist[3])
				+gamma*area*(termlist[0]*q_raw**2))
			ratio = hel/((vibe*q_raw+machine_eps)/(np.exp(vibe*q_raw)-1)+machine_eps)
			return residual(ratio[band])

		fit = scipy.optimize.minimize(objective,x0=tuple(initial_conditions),method=minmethod,
			callback=callback)

		return {'master':master,'core_cfs':core_cfs}

	def gridded(self):
		"""Probe curvature in a regular grid."""
		global work,callback,Nfeval,name_groups,objective,nprots
		sns = work.sns()
		#---! still settling on the formation. previously swept extents here
		extent = self.spec['spacer']/2.0
		package = dict([(sn,{'data':self.master_decorator(sn,extent)}) for sn in work.sns()])
		#---optimize for each simulation
		for sn in sns:
			###### nprots = work.meta[sn].get('nprots',1)
			method = ['fmin','minimize','basinhopping','differential_evolution','brute'][1]
			minmethod = ['L-BFGS-B','CG','COBYLA','dogleg','Nelder-Mead','SLSQP'][-1]
			initial_conditions = [25.0*2,0.0,-0.1]+[0.0 for i in range(nprots)]
			bounds_dict = {'kappa':(0.0,100.0),'gamma':(-10,10),'vibe':(-10,10),'curvature':(-0.05,0.05)}

			#---run
			Nfeval = 1
			start_time = time.time()
			name_groups = ['kappa','gamma','vibe']+['curve(%d)'%i for i in range(nprots)]
			bounds = [bounds_dict[k] for k in ['kappa','gamma','vibe']]
			bounds += [bounds_dict['curvature'] for k in range(nprots)]
			objective = package[sn]['data']['master']

			status('starting to optimize %s'%sn,tag='compute')
			if method == 'fmin':
				fit = scipy.optimize.fmin(objective,x0=tuple(initial_conditions),
					disp=True,callback=callback,full_output=True,)
			elif method == 'minimize':
				fit = scipy.optimize.minimize(objective,x0=tuple(initial_conditions),method=minmethod,
					callback=callback)
				status('finished\n'+str(fit),tag='result')
			elif method == 'basinhopping':
				def callback(args): 
					print "step = %d, args = %s"%(Nfeval,str(args))
					Nfeval += 1
				fit = scipy.optimize.basinhopping(objective,x0=tuple(initial_conditions),
					disp=True,callback=callback)
			elif method == 'differential_evolution':
				fit = scipy.optimize.differential_evolution(objective,bounds=bounds,disp=True,
					callback=callback)
			elif method == 'brute':
				fit = scipy.optimize.brute(objective,ranges=bounds,Ns=5,disp=True)
			else: raise Exception()
			status('%.1fmin elapsed'%((time.time()-start_time)/60.),tag='time')

			#---!!!!!!!!!!!!!!!!!!!!!!!1
			import ipdb;ipdb.set_trace()

	###---OPTIMIZED SOLUTIONS

	def drop_gaussians(self,**kwargs):
		"""
		Method for choosing the positions of Gaussians.
		"""
		pos_spec = kwargs.get('curvature_positions',{})
		method = pos_spec.get('method',None)
		extent = kwargs.get('extents',{}).get('extent',{})
		if not method: raise Exception('need a method for setting the curvature fields')
		elif method=='protein_subselection':
			self.data_prot,_ = plotload('protein_abstractor')
			for sn in work.sns():
				selections = pos_spec.get('selections',None)
				if not selections: raise Exception('need selections in protein_subselection')
				#---determine the centers of the protein according to the selections
				#---...noting that the protein_abstractor points are stored by the residue, not bead/atom 
				points = np.array([np.transpose(self.data_prot[sn]['data']['points'],(1,0,2))[s] 
					for s in selections])
				#points = np.transpose(self.data_prot[sn]['data']['points'],(1,0,2))[selections]
				points = points.mean(axis=1)[...,:2]
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
				status('computing curvature fields for %s'%sn)
				incoming = basic_compute_loop(make_fields,looper=looper)
				#---! inelegant
				for ii,(fr,ndrop) in enumerate(reindex): fields_unity[fr][ndrop] = incoming[ii]
				self.memory[(sn,'fields_unity')] = fields_unity

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
		#(1.0-2.0*(xxx<0.0))
		#import ipdb;ipdb.set_trace()
		#return np.array([np.transpose(cfs,(1,0,2,3))[cc]*c 
		#	for cc,c in enumerate(curvatures)]).max(axis=0)

	def wilderness(self,**kwargs):
		"""
		Wandering on the error landscape to optimize the curvature coupling hypothesis.
		"""
		bundle,solutions = {},{}
		global Nfeval
		spec = self.spec
		extents = spec.get('extents',{})
		extents_method = extents.get('method')
		curvature_sum_method = self.spec['curvature_sum']
		#---prepare curvature fields
		if not extents_method: raise Exception('set extents method')
		elif extents_method=='fixed_isotropic': self.drop_gaussians(**spec)
		#---if the output file exists then the optimization is done and we just keep the curvature fields
		out_fn_name = 'curvature_%s.dat'%self.signifier
		if os.path.isfile(os.path.join(work.postdir,out_fn_name)): return

		#---optimize over simulations
		for snum,sn in enumerate(work.sns()):
			status('starting optimization for %s %d/%d'%(sn,snum+1,len(work.sns())),tag='optimize')

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

			tweak = {}
			#---! remove this eventually
			signterm = tweak.get('inner_sign',-1.0)
			lowcut = kwargs.get('lowcut',tweak.get('lowcut',0.0))
			band = cctools.filter_wavevectors(q_raw,low=lowcut,high=tweak.get('hicut',1.0))
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

			def objective(args):
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
				hel = (kappa*area*(termlist[0]*q_raw**4+signterm*termlist[1]*q_raw**2
					+signterm*termlist[2]*q_raw**2+termlist[3])
					+gamma*area*(termlist[0]*q_raw**2))
				ratio = hel/((vibe*q_raw+machine_eps)/(np.exp(vibe*q_raw)-1)+machine_eps)
				return residual(ratio[band])

			Nfeval = 0
			initial_conditions = [25.0*2,0.0,-0.1]+[0.0 for i in range(ndrops)]
			fit = scipy.optimize.minimize(objective,
				x0=tuple(initial_conditions),method='SLSQP',callback=callback)
			#---package the result
			bundle[sn] = dict(fit,success=str(fit.success))
			solutions['%s_%s'%(sn,'x')] = bundle[sn].pop('x')
			solutions['%s_%s'%(sn,'jac')] = bundle[sn].pop('jac')
			solutions['%s_%s'%(sn,'cf')] = np.array(self.curvature_sum(cfs,fit.x[3:],
				method=curvature_sum_method).mean(axis=0))
		try: store(obj=solutions,name=out_fn_name,path=work.postdir,attrs=dict(bundle=bundle,spec=self.spec))
		except: 
			import ipdb;ipdb.set_trace()

#import matplotlib as mpl;import matplotlib.pyplot as plt;plt.imshow(np.array(self.curvature_sum(cfs,fit.x[3:]).mean(axis=0)).T,interpolation='nearest',origin='lower');plt.show()

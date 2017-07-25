#!/usr/bin/env python

"""
CURVATURE-UNDULATION COUPLING (version 3)
"""

import copy,importlib,time,os
from codes.hypothesizer import hypothesizer
from base.tools import status
from omnicalc import store,load
import numpy as np
import multiprocessing as mp
import codes.curvature_coupling.tools as cctools
from codes.curvature_coupling.tools import manyjob,manyjob_worker
from codes.curvature_coupling.tools import construct_curvature_fields_trajectory

###---CONTROLLER

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

class InvestigateCurvature:

	"""
	Manage the entire curvature coupling calculation.
	"""

	def __init__(self,**kwargs):
		"""Assemble the necessary modules and execute."""
		if kwargs: raise Exception('unprocessed kwargs: %s'%kwargs)
		#---get the data according to instructions in the metadata
		self.fetch_data_and_specs()
		self.memory = self.loader_func(self.data)
		#---master loop over the attempts
		global work
		for signifier,attempt in work.plots.get('curvature',{}).get('specs',{}).get('attempts',{}).items():
			self.prepare_rootdir('curvature_%s'%signifier)
			self.prepare_database()
			#---cursors
			self.sig,self.spec = signifier,attempt
			#---select the hypotheses
			self.style = attempt.get('style',None)
			if self.style=='manual': self.manual()
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
		self.data,self.calc = plotload('curvature')
		self.meta_spec = work.plots.get('curvature',{}).get('specs',{})
		#---load the module responsible for packing the data
		loader_spec = self.meta_spec.get('loader',None)
		if not loader_spec: raise Exception('you must add specs,loader to the plot entry for cuvature')
		self.loader_func = self.gopher(loader_spec,module_name='module',variable_name='function')
		#---load the module responsible for packing the data
		database_spec = self.meta_spec.get('design',None)
		if not database_spec: raise Exception('you must add specs,loader to the plot entry for cuvature')
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

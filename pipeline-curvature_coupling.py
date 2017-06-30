#!/usr/bin/env python

import time
from codes import curvature_coupling as couplecalc
#---we manually store/load data instead of using omnicalc
from omnicalc import store,load
import multiprocessing as mp
import scipy

#---switches
do_single = False
do_strict = False
show_optimization_log = False
do_next = 'plotting'
eps = np.finfo(np.float32).eps
huge_number = 1./eps

#---everything must be perfect (avoid division-by-zero errors)
#---! this is too strict. fails on sphinx.
if do_strict:
	import warnings
	warnings.filterwarnings('error')

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

###

#---hard-coded version controls
roundsig = 'hcv3'
tweak_menu = []
#---wish fulfillment? (note that 40 is the original 29)
#---turned out to be correct (thank the optimizer)
tweak_menu += [('v41',{
	'residual_form':'log',
	'sweep_code':'sweep-v1',
	'energy_form':energy_form,
	'curvature_field_spec':'cfv5',
	'inner_sign':-1.0,'fft_flag':'complex',
	'init':[20,0.0,2.0],'hicut':2.0,'opt':optimize_hypothesis,
	'cut_curve':(0.0,0.032),'error_ceiling':0.0012,
	'kappa_gamma_prefactor':{'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}})]
tweak_menu += [('v42',{
	'residual_form':'log',
	'sweep_code':'sweep-v2',
	'energy_form':energy_form,
	'curvature_field_spec':'cfv5',
	'inner_sign':-1.0,'fft_flag':'complex',
	'init':[20,0.0,2.0],'hicut':1.0,'opt':optimize_hypothesis,
	'cut_curve':(0.0,0.032),'error_ceiling':0.0014,
	'kappa_gamma_prefactor':{'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}})]
#---! THIS IS THE CURRENT CALCULATION FOR THE PNAS DRAFT
tweak_menu += [('v40',)]

#---select a hard-coded version
re_hcv = '^hcv=(.+)$'
if False: hcv_flags = [i for i in sys.argv if re.match(re_hcv,i)]
#---reworking the code with the "good" version
else: hcv_flags = ['hcv=v40']
hcv_index = zip(*tweak_menu)[0].index(re.findall(re_hcv,hcv_flags[0])[0])
hcv,tweak = tweak_menu[hcv_index]

#---curvature field dataset is specified by the tweak
cfsig = tweak['curvature_field_spec']
tweak_cf = {
	#---the cfv4 was first generated with hcv2 and linked in to the following codes because it is 32GB
	'cfv4':[{'motion':j,'mapping':i} for i in ['single','protein'] for j in ['static','dynamic']],
	#---ignore static for now despite the cost of dynamic
	'cfv5':[{'motion':j,'mapping':i} for i in ['single','protein'] for j in ['static','dynamic'][1:]],
}[cfsig]

#---set the motion of the curvature fields
tweak_cf = [{'motion':j,'mapping':i} 
	for i in ['single','protein'] for j in ['static','dynamic'][1:]]

#---specify the energy function, prefactors, and default behaviors
tweak = {
	'residual_form':'log',
	'sweep_code':'sweep-v1',
	'energy_form':energy_form,
	'init':[20,0.0,10.0],
	'hicut':1.0,'opt':optimize_hypothesis,
	'kappa_gamma_prefactor':{
		'kappa':1/2.0,'gamma':1/2.0,'vibe':1.0}}

#---scripting sequence logic
goto = 'plotting'
load_hqs = False
quit_on_complete = False
if 'noplot' in sys.argv: goto = None
elif 'single' in sys.argv: goto = 'single'
elif 'single' in sys.argv: goto = 'single'
if 'qoc' in sys.argv: quit_on_complete = True
#---plot while calculating
skip_calculation = False
#---keep midplanes in memory
debug_carefully = True
if goto == 'single':
	skip_calculation,load_hqs,debug_carefully = True,True,True

#---DATABASE
#-------------------------------------------------------------------------------------------------------------

absolute_start_time = time.time()
#---! see above for: from codes import couplecalc
couplecalc.tweak = tweak

#---load upstream data
rootdir_cc = 'curvature_coupling-%s-%s'%(roundsig,hcv)
rootdir_cf = 'curvature_coupling-%s-%s'%(roundsig,cfsig)
data,calc = plotload('curvature_coupling',work)

#---root directory
rootdir = os.path.join(work.paths['post_data_spot'],rootdir_cc,'')
if not os.path.isdir(rootdir): os.mkdir(rootdir)
rootdir_cf = os.path.join(work.paths['post_data_spot'],rootdir_cf,'')
if not os.path.isdir(rootdir_cf): os.mkdir(rootdir_cf)

#---sqlalchemy
from sqlalchemy import Table,Column,Float,String,Integer,create_engine
from sqlalchemy.orm import sessionmaker,mapper
engine_cc = create_engine('sqlite:///%s/hypotheses.sqlite'%rootdir)
Session_cc = sessionmaker(bind=engine_cc)
engine_cf = create_engine('sqlite:///%s/curvature_fields.sqlite'%rootdir_cf)
Session_cf = sessionmaker(bind=engine_cf)
namer_cf = lambda pk : os.path.join(rootdir_cf,'curvature_field.%d.dat'%pk)
namer_cc = lambda pk : os.path.join(rootdir,'hypothesis.%d.dat'%pk)
couplecalc.namer_cc = namer_cc
couplecalc.namer_cf = namer_cf

#---DEFINITIONS
#-------------------------------------------------------------------------------------------------------------

database = {}
database['hypothesis'] = dict(
	table_name = 'hypothesis',
	data_map = (
		('sn',str),
		('timecode',str),		
		('curvature',float),
		('sigma_a',float),
		('sigma_b',float),
		('isotropy',float),
		('mapping',str),
		('motion',str),
		('theta',float),
		('vibe',float),
		('gamma',float),
		('kappa',float),
		('error',float),
		('fallback',str),
		),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein','multipoint'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],
		],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0 or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},
		],
	anticipate = ['kappa','gamma','vibe','error'],
	)

database['field'] = dict(
	table_name = 'field',
	data_map = (
		('sn',str),
		('timecode',str),		
		('curvature',float),
		('sigma_a',float),
		('sigma_b',float),
		('isotropy',float),
		('mapping',str),
		('motion',str),
		('theta',float),
		),
	test_spec = [
		lambda x:'sigma_a' in x,
		lambda x:x['mapping'] in ['single','protein'],
		lambda x:x['motion'] in ['static','dynamic'],
		lambda x:x['isotropy'] or x['sigma_b'],
		],
	change_spec = [
		lambda x:{'isotropy':x['sigma_a']/x['sigma_b']} if not x['isotropy'] else {},
		lambda x:{'sigma_b':x['sigma_a']/x['isotropy']} if not x['sigma_b'] else {},
		lambda x:{'sigma_a':1.0,'sigma_b':1.0,'isotropy':1.0,'mapping':'single','curvature':0.0} 
			if (x['curvature']==0.0  or (x['sigma_a']==0.0 and x['sigma_b']==0.0)) else {},
		lambda x:{'curvature':1.0} if x['curvature']!=0.0 else {},
		lambda x:{'sn':x['fallback'],'fallback':'none'} if work.meta[x['sn']]['nprots']==0 
			and x['curvature'] != 0.0 and x['mapping']=='protein' else {},
		#---test specs throw errors when the fail but change specs simply "reduce" like hypotheses
		#---incorrect: lambda x:{'mapping':'single'} if work.meta[x['sn']]['nprots']==0 else {},
		],
	anticipate = ['kappa','gamma','vibe','error'],
	)

#---database interface (now separated by table)
sessions = {}
database['hypothesis']['engine'],database['hypothesis']['session_maker'] = engine_cc,Session_cc
database['field']['engine'],database['field']['session_maker'] = engine_cf,Session_cf
for table_name in database:
	sessions[table_name],outgoing = couplecalc.magic_data({table_name:database[table_name]},
		database[table_name]['engine'],database[table_name]['session_maker'])
	for key in outgoing: globals()[key.__name__] = key
session_makers = dict([(key,database[key]['session_maker']) for key in ['hypothesis','field']])

#---PREPARE CALCULATIONS
#-------------------------------------------------------------------------------------------------------------

if tweak['sweep_code'] == 'sweep-v1':

	#---define sweeps
	hypotheses = []
	#---HACKING FOR A NEGATIVE CURVATURE INDUCER
	base_range_sweep = [0,1,2,3,4,5,6,7,8,9,10,12,18,24,32,64]
	base_range = [0.0,0.005,0.01,0.014,0.018,0.02,0.024,0.028,0.032,0.04,0.05,0.06,0.07,0.08,0.09,0.1]
	curvature_sweep = np.sort(np.concatenate((-1.0*np.array(base_range),1.0*np.array(base_range),)))
	if 'extended_curvatures' in tweak and tweak['extended_curvatures'] == 'v1':
		curvature_sweep += [0.06,0.07,0.08,0.09,0.1,0.2]+[0.4,0.6,0.8,1.0]
	elif 'extended_curvatures' in tweak and tweak['extended_curvatures'] == 'v2':
		curvature_sweep += [0.06,0.07,0.08,0.09,0.1,0.2]
	elif 'extended_curvatures' in tweak and tweak['extended_curvatures'] == 'v3':
		curvature_sweep += [0.1]
	sweep_curvatures_single = dict(route=['curvature'],values=curvature_sweep)
	sweep_extents_single = dict(route=['sigma_a'],values=[float(j) for j in base_range_sweep])
	sweep_curvatures_pd = dict(route=['curvature'],values=curvature_sweep)
	#---! beefed up the dynamic since it should match the single
	sweep_extents_pd = dict(route=['sigma_a'],values=[float(j) for j in base_range_sweep])

	#---complete the hypotheses with further sweep specification
	#---outer loop over motions and inner if-loop over mappings
	#sns = data.keys()
	sns = work.sns()
	if len(set(v['specs']['slice_name'] for k,v in calc.items()))!=1:
		raise Exception('inconsistent slice names in upstream data from: %s'%calc.keys())
	slice_name = [calc['undulations']['specs']['slice_name']]
	sns_unfree = [s for s in sns if work.meta[s]['nprots']>0]
	#---also save some of these hypotheses in groups according to the tweak_cf
	hypo_groups = {i:[] for i in range(len(tweak_cf))}
	for tweak_cf_subset in tweak_cf:
		group_key = tweak_cf.index(tweak_cf_subset)
		motion = tweak_cf_subset['motion']
		if tweak_cf_subset['mapping'] == 'single':
			#---the hypothesizer does the product/right-prism hypotheses in the phase space
			hypotheses_single = couplecalc.hypothesizer(sweep_curvatures_single,sweep_extents_single,
				default={'mapping':'single','isotropy':1.0,'fallback':'none','motion':motion,'theta':0.0})
			#---loop single fields over all simulations
			for hypo in hypotheses_single:
				for sn in sns:
					hypo_new = dict(hypo)
					hypo_new.update(sn=sn,timecode='>'.join(slice_name))
					hypotheses.append(hypo_new)
					hypo_groups[group_key].append(hypo_new)
		if tweak_cf_subset['mapping'] == 'protein':
			#---loop over simulations with "fallback" for free bilayer
			hypotheses_pd = couplecalc.hypothesizer(sweep_curvatures_pd,sweep_extents_pd,
				default={'mapping':'protein','isotropy':1.0,'fallback':'none','motion':motion,'theta':0.0})
			for hypo in hypotheses_pd:
				for sn in sns:
					if sn not in sns_unfree:
						for sn_alt in sns_unfree:
							hypo_new = dict(hypo)
							hypo_new.update(sn=sn,fallback=sn_alt,timecode='>'.join(slice_name))
							hypotheses.append(hypo_new)
							hypo_groups[group_key].append(hypo_new)
					else:
						hypo_new = dict(hypo)
						hypo_new.update(sn=sn)
						hypo_new.update(sn=sn,timecode='>'.join(slice_name))
						hypotheses.append(hypo_new)
						hypo_groups[group_key].append(hypo_new)

else: raise Exception('unclear sweep code in tweak')

#---PREPARE COMPUTE
#-------------------------------------------------------------------------------------------------------------

#---load the database
start = time.time()
session = sessions['hypothesis']
for hh,hypo in enumerate(hypotheses):
	status('populating hypos',tag='load',i=hh,looplen=len(hypotheses),start=start)
	#---reduce step before checking database
	hypo_full = Hypothesis(**hypo)
	matches = session.query(Hypothesis).filter_by(**hypo_full.base()).all()
	if not any(matches): session.add(hypo_full)
session.commit()
session = sessions['field']
for hh,hypo in enumerate(hypotheses):
	status('populating fields',tag='load',i=hh,looplen=len(hypotheses),start=start)
	#---reduce step before checking database
	hypo_full = Field(**hypo)
	matches = session.query(Field).filter_by(**hypo_full.base()).all()
	if not any(matches): 
		session.add(hypo_full)
session.commit()

#---integrity checks on database rows
hypotheses_reduced = [i.dict() for i in sessions['hypothesis'].query(Hypothesis).all()]
fields_reduced = [i.dict() for i in sessions['field'].query(Field).all()]
assert not [i for i in hypotheses_reduced if i['mapping']=='protein' and i['curvature']==0.0]
assert not [i for i in hypotheses_reduced if i['curvature']==0.0 
	and not (i['sigma_a']==1.0 and i['isotropy']==1.0 and i['sigma_b']==1.0)]

#---point heights into "memory"
status('populating memory',tag='load')
memory = {}
for sn in sns:
	if (sn,'hqs') not in memory and (not skip_calculation or load_hqs):
		dat = data['undulations'][sn]['data']
		vecs = dat['vecs']
		mesh = dat['mesh']
		#---development test using previous dataset
		if tweak.get('old_data',False)=='v650_only' and sn=='membrane-v650-enthx4-dev':
			status('loading previous data',tag='backwards')
			v650fns = ['postproc.structures.membrane-v650.s06.%d-%d-160.grid_spacing-0.5.dat'%(
				i*10**6,(i+1)*10**6) for i in range(1,5+1)]
			from base.store import load
			old_dir = '/home/rpb/worker/repo-postproc/project-curve/'
			olds = [load(fn,path=old_dir) for fn in v650fns]
			oldcat = np.concatenate([np.array([[o['%d.%d'%(fr,mn)] for mn in range(2)] 
				for fr in range(o['nframes'])]) for o in olds])
			vecs = np.concatenate([o['vecs'] for o in olds])
			mesh = np.transpose(oldcat,(1,0,2,3))
		midplane = mesh.mean(axis=0)
		zmeans = midplane.reshape((midplane.shape[0],-1)).mean(axis=1)
		midplane = np.array([i-zmeans[ii] for ii,i in enumerate(midplane)])
		hqs = couplecalc.fft_field(midplane)
		#---! note that this is temporary -- might be wasteful
		if debug_carefully: memory[(sn,'hs')] = midplane
		memory[(sn,'hqs')] = hqs
		memory[(sn,'vecs')] = vecs

#---export to the module
couplecalc.Hypothesis = Hypothesis
couplecalc.Field = Field
couplecalc.rootdir = rootdir
couplecalc.rootdir_cf = rootdir_cf
couplecalc.memory = memory
couplecalc.calc = calc
couplecalc.data = data
couplecalc.work = work

#---POOL and MULTIPROCESS
#-------------------------------------------------------------------------------------------------------------

from codes.curvature_coupling import construct_curvature_fields_trajectory
#---compute pending fields according to populated rows
fns = [(i.id,namer_cf(i.id)) for i in sessions['field'].query(Field).all()]
pending = [(pk,fn) for pk,fn in fns if not os.path.isfile(fn)]
if pending and not skip_calculation:
	#---loop over absent files
	start = time.time()
	for ii,(pk,fn) in enumerate(pending):
		status('computing curvature field',tag='compute',i=ii,looplen=len(pending),start=start)
		hypo = sessions['field'].query(Field).filter_by(id=pk).one().dict()
		sn = hypo['sn']
		dat = data['undulations'][sn]['data']
		vecs = dat['vecs']
		mn = np.shape(dat['mesh'])[2:]
		fields = construct_curvature_fields_trajectory(vecs=vecs,mn=mn,**hypo)
		store({'fields':np.array(fields['fields'])},os.path.basename(fn),rootdir_cf,
			attrs={key:val for key,val in fields.items()+hypo.items() if key!='fields'},verbose=False)

#---solve the hypotheses
#---for memory efficiency we queue up hypotheses according to which curvature field they require
#---note that we had a simpler, memory-hogging loop in a previous iteration of this code
fns = [(i.id,namer_cc(i.id)) for i in sessions['hypothesis'].query(Hypothesis).all()]
pending = [(pk,fn) for pk,fn in fns if not os.path.isfile(fn)]
if pending and not skip_calculation:
	hypotheses = [sessions['hypothesis'].query(Hypothesis).filter_by(id=pk).one().dict() 
		for pk in zip(*pending)[0]]
	fields_required = [sessions['field'].query(Field).filter_by(**f.dict()).one() 
		for f in [Field(**h) for h in hypotheses]]
	field_ids_by_hypothesis = np.array([f.id for f in fields_required])
	unique_field_ids = np.unique(field_ids_by_hypothesis)
	#---compute the curvatures in batches
	for uu,ufid in enumerate(unique_field_ids):
		status('computing all hypotheses for field %d/%d'%(uu,len(unique_field_ids)),tag='compute')
		hypo_subset = [hypotheses[j] for j in np.where(field_ids_by_hypothesis==ufid)[0]]
		key_cf = ('curvature',ufid)
		memory[key_cf] = load(os.path.basename(namer_cf(ufid)),cwd=os.path.dirname(namer_cf(ufid)))
		#---queue for each part of the computation
		queue_hypothesis = mp.Queue()
		#---solve
		couplecalc.manyjob(single=do_single,
			function=couplecalc.worker,
			queue=queue_hypothesis,
			session_classes=session_makers,
			objects=hypo_subset,
			kwargs={'preloaded_curvature':True})
		#---clear that hypothesis from memory
		del memory[key_cf]
	status('done all batches',tag='compute')

#---order the simulations
sns_all,sns_unfree = sns,[]
status('%.1fmin elapsed'%((time.time()-absolute_start_time)/60.),tag='time')
if do_next=='plotting': next_script = ns = 'calcs/pipeline-curvature_coupling_plot.py'
elif do_next=='single': raise Exception('dev. go rescue this code from legacy')
elif not do_next: sys.exit()
status('running "%s"'%next_script,tag='goto')
exec(open(next_script).read())

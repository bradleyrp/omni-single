#!/usr/bin/env python

"""
DRILLDOWN to investigate curvature coupling results
"""

from codes.undulate import calculate_undulations
from codes.undulate_plot import undulation_panel,add_undulation_labels,add_axgrid,add_std_legend
from codes.curvature_coupling.curvature_coupling_plots import individual_reviews_plotter
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature,prepare_oscillator_function
from codes.hypothesizer import hypothesizer
from codes.undulate import perfect_collapser,blurry_binner
from datapack import asciitree
from omnicalc import store
import scipy
import scipy.optimize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
machine_eps = eps = np.finfo(float).eps
import copy,json,glob
import h5py

#---where to catalog the completed optimizations
coordinator = 'curvature_coupling_drilldown.json'

"""
history:
	started with protein_dynamic_single_catchall and had to manually fix the coordinator
	added protein_dynamic_single_uniform
	added extents (v3,v4) and selected only protein_dynamic_single_uniform to save time
		the protein_dynamic_single runs for thousands of steps
	added v5,v6 to curvature.yaml and generated fields via `make compute` for all frames
	current batch of six tests per simulation is worth making a summary plot
		includes two extents for protein_dynamic_single and four for protein_dynamic_single_uniform
		otherwise lots of other combinations to try
		but currently working on the summary plots
"""

#---see history notes above. the round name selects different optimmization tests
round_name = ['r1','r2','r3','r4'][0]

#---all tags in the design loop in curvature.yaml
master_tags = {
	#---several rounds use blurry_explicit binner
	#---first round uses protein_dynamic_single
	'r1':['v1_fixed_extent_1','v2_fixed_extent_2'],
	#---second round uses protein_dynamic_single_uniform
	'r2':['v1_fixed_extent_1','v2_fixed_extent_2','v3_fixed_extent_0.5','v4_fixed_extent_4'],
	'r3':['v5_fixed_extent_all_frames_1','v6_fixed_extent_all_frames_2'],
	'r4':['v1_fixed_extent_1','v2_fixed_extent_2','v3_fixed_extent_0.5','v4_fixed_extent_4'],}

@autoload(plotrun)
def loader():
	"""Load data."""
	#---only load once
	if 'data' not in globals():
		#---begin loading sequence
		plotname = 'curvature_undulation_coupling'
		if plotname not in work.plots: raise Exception('add %s to the plots metadata'%plotname)
		plotspecs = work.plots[plotname].get('specs',{})
		calcname = plotspecs.get('calcname',plotname)
		#---new method for getting all upstream calculations in the loop
		combodat = work.collect_upstream_calculations_over_loop(plotname)
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

def individual_reviews_drilldown_height_summary():
	"""
	Plot height profiles once per curvature_undulation_coupling_calculation.
	"""
	plotspec = {
		'coupling_review_drilldown_heights':{
			'viewnames':['average_height','average_height_pbc','neighborhood_static',
				'neighborhood_dynamic','average_field','example_field','example_field_pbc',
				'spectrum','spectrum_zoom'][:4]},}
	global seepspace
	seep = dict([(key,globals()[key]) for key in seepspace])
	for out_fn,details in plotspec.items(): 
		#---! this currently only plots four tiles none of which have anything to do with the field
		#---! ...hence there appear to be repeats since datas is populated with different extents
		individual_reviews_plotter(out_fn=out_fn,seep=seep,**details)

def prepare_hypotheses():
	"""
	Generate a list of hypotheses for testing the optimizer.
	We generate the hypotheses in several "rounds" during development of the drilldown.
	"""
	#---start making hypotheses here with the blurry_explicit method
	if round_name in ['r1','r2','r3']:
	    #---default parameters
		opt_params = dict(design={
			'optimize_method':'Nelder-Mead',
			'store_instantaneous_fields':False,
			'store_instantaneous_fields_explicit':False},
			fitting={'initial_kappa':20.0})
		opt_params['design']['binner'] = 'explicit'
		opt_params['design']['weighting_scheme'] = 'blurry_explicit'
		if round_name=='r1':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single']},
				#---! include all elements of curvature_positions or it gets overwritten
				{'route':('design','curvature_positions','nframes'),'values':[100]},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r2':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('design','curvature_positions','nframes'),'values':[100]},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r3':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		elif round_name=='r4':
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
		else: raise Exception('fail')
	elif round_name in ['r4']:
		if round_name=='r4':
		    #---default parameters
			opt_params = dict(design={
				'optimize_method':'Nelder-Mead',
				'store_instantaneous_fields':False,
				'store_instantaneous_fields_explicit':False},
				fitting={'initial_kappa':20.0})
			opt_params['design']['binner'] = 'explicit'
			args = [{'route':('design','curvature_positions','method'),
				'values':['protein_dynamic_single_uniform']},
				{'route':('fitting','high_cutoff'),'values':[1.0]},]
			hypos = hypothesizer(*args,default=opt_params)
	else: raise Exception('fail')
	return hypos

def new_coordinator(fn):
	"""Make a new table."""
	with open(os.path.join(work.postdir,fn),'w') as fp: fp.write(json.dumps({}))

def fetch_hypothesis_table(fn):
	"""Load a master table of completed hypotheses from disk."""
	with open(os.path.join(work.postdir,fn),'r') as fp: return json.loads(fp.read())

def write_hypothesis_table(fn,table):
	"""Load a master table of completed hypotheses from disk."""
	with open(os.path.join(work.postdir,fn),'w') as fp: return fp.write(json.dumps(table))

def prepare_calcs(tags,hypos,sns):
	"""Cross hypotheses with simulations and curvature fields."""
	table = []
	for tag in tags:
		for sn in sns:
			for hnum,hypo in enumerate(hypos):
				status('computing %s,%s with the following hypothesis'%(tag,sn),tag='compute')
				asciitree(hypo)
				if False: ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d.%H%M')
				#---add design back to the hypothesis
				#---note that at least one of the following three deepcopy are necessary
				hypo_more = dict(design=copy.deepcopy(calcs[tag][sn]['calcs']['specs']['design']))
				hypo_more['design'].update(**copy.deepcopy(hypo['design']))
				hypo_more['fitting'] = copy.deepcopy(hypo['fitting'])
				table.append(dict(hypo=hypo_more,sn=sn,tag=tag))
	return table

def calc_namer(**item):
	"""Name a calculation."""
	tag = item['tag']
	sn = item['sn']
	hypo = item['hypo']
	fn = 'curvature_coupling_drilldown.%s.%s'%(tag,work.namer.short_namer(sn,spot='sims'))
	return fn

def run_optimizer(fn,**calc):
	"""
	Run the optimizer once and save the results.
	"""
	tag,sn,hypo = calc['tag'],calc['sn'],calc['hypo']
	dat = datas[tag][sn]
	remember = {
		(sn,'drop_gaussians_points'):dat['drop_gaussians_points'],
		(sn,'sampling'):dat['sampling'],
		(sn,'fields_unity'):dat['fields_unity']}
	ic = InvestigateCurvature(sn=sn,work=work,
		remember=remember,
		protein_abstractor=data['protein_abstractor'][sn]['data'],
		undulations=data['undulations'][sn]['data'],
		**hypo)
	#---package and save by timestamp for now
	attrs,result = copy.deepcopy(ic.finding['attrs']),copy.deepcopy(ic.finding['result'])
	#---final_simplex is not saveable
	try: attrs['bundle'][sn].pop('final_simplex')
	except: pass
	store(obj=result,name='%s.dat'%fn,path=work.paths['post_data_spot'],
		attrs=attrs,verbose=True)

def get_drilldown(fn):
	"""Load data from a drilldown optimization."""
	fl = h5py.File(os.path.join(work.paths['post_data_spot'],'%s.dat'%fn),'r')
	dat = dict([(key,np.array(val)) for key,val in fl.items()])
	fl.close()
	return dat

def plot(name,table):
	"""
	Summary plots for the drilldown.
	"""
	row = table[name]
	item = get_drilldown(name)
	figsize = (10,10)
	axes,fig = square_tiles(1,figsize)
	high_cutoff = table[name]['hypo']['fitting']['high_cutoff']
	ax = axes[0]
	sn = row['sn']
	qs = item['qs']
	#---annotations for blurry_explicit weighting scheme
	if row['hypo']['design'].get('weighting_scheme',None)=='blurry_explicit':
		color = '0.5'
		ratios = item['ratios']
		q_red,ratios_red,q_mapping = blurry_binner(qs,ratios,return_mapping=True)
		ax.scatter(q_red,ratios_red,color='r',marker='o',s=20,zorder=3)
		ax.set_title('%s (weighted residuals)'%sn)
	else: 
		ax.set_title('%s'%sn)
		color = 'k'
	#---plot the oscillator function
	oscillator_function = prepare_oscillator_function()
	ax.scatter(qs,oscillator_function(item['x'][2],qs),marker='x',s=10,color=color,zorder=2)
	#---plot the data
	ax.scatter(qs,item['ratios'],marker='o',s=10,color=color,zorder=2)
	ax.axvline(high_cutoff,color=color)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axhline(1.0,c='k')
	#---inset shows the curvature field
	fields_unity = datas[row['tag']][row['sn']]['fields_unity']
	cf_frame = fields_unity[0]
	fields = np.sum([item['x'][3:][ii]*cf_frame[ii] for ii in range(len(cf_frame))],axis=0)
	axins = inset_axes(ax,width="35%",height="35%",loc=1)
	mag_max = np.abs(fields).max()
	im = axins.imshow(fields.T,origin='lower',interpolation='nearest',
		vmin=mag_max*-1,vmax=mag_max,cmap=mpl.cm.__dict__['RdBu_r'])
	axins.set_xticks([])
	axins.set_yticks([])
	axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	axins.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	cax = inset_axes(axins,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=axins.transAxes,borderpad=0)
	cbar = plt.colorbar(im,cax=cax,orientation="vertical")
	cax.tick_params(axis='y',which='both',left='off',right='off',labelleft='off')
	cax.set_title('%0.3f'%np.abs(fields).max())
	ax.set_ylabel('${H}_{hel}/{H}_{osc}$',fontsize=20)
	ax.set_xlabel('$\| q \|\,({nm}^{-1})$',fontsize=20)
	cax.set_ylabel('$C_0\,({nm}^{-1})$',rotation=-90)
	#---extra details
	label = '\n'.join(['%s = %s'%(key,('%'+'0.%d'%prec+'f')%get_result(prop=key,name=name)) 
		for key,prec in [('kappa',1),('gamma',2),('vibe',5),('error',5)]])
	ax.text(0.05,0.95,label,transform=ax.transAxes,size=16,verticalalignment='top',ha='left')
	#---save the figure
	picturesave('fig.%s'%name,work.plotdir,backup=False,version=False,meta={})

def get_result(prop,name):
	"""Interpret results from the optimizer."""
	if prop=='kappa': return results[name]['x'][0]
	elif prop=='gamma': return results[name]['x'][1]
	elif prop=='vibe': return np.abs(results[name]['x'][2])
	elif prop=='error': 
		ratios = results[name]['ratios']
		qs = results[name]['qs']
		high_cutoff = table[name]['hypo']['fitting']['high_cutoff']
		band = np.all((qs>=0.0,qs<=high_cutoff),axis=0)
		values = ratios[band]
		error = np.sum(np.log10(values.clip(min=machine_eps))**2)/float(len(values))
		return error
	elif prop=='curvature': 
		curvs = results[name]['x'][3:]
		sign = 1.-2*(curvs[np.argmax(np.abs(curvs))]>0.)
		mag = np.abs(curvs).max()
		return sign*mag
	elif prop=='extent':
		return table[name]['hypo']['design']['extents']['extent']
	else: raise Exception('not sure how to get %s'%name)

def summary_plots():
	"""
	Summarize fitted parameters and spectra for several optimizations.
	"""
	table = fetch_hypothesis_table(coordinator)
	sns = work.sns()
	props = ['kappa','gamma','vibe','error','curvature','extent']
	ntiles = len(sns)*len(props)
	#---divide tests by simulation
	keysets = dict([(sn,sorted([key for key in table if table[key]['sn']==sn])) for sn in sns])

	def ncolors(n,cmap='jet'):
		"""Make a series of colors."""
		return [mpl.cm.__dict__[cmap](i) for i in np.linspace(0,1.,n)]

	#---plot properties on a large set of tiles
	for snum,sn in enumerate(sns):
		axes,fig = square_tiles(len(props),figsize=(8,8),hspace=0.4,wspace=0.4,favor_rows=True)
		colors = ncolors(len(keysets[sn]))
		for pnum,prop in enumerate(props):
			ax = axes[pnum]
			for nnum,name in enumerate(keysets[sn]):
				vals = [get_result(prop=prop,name=name)]
				shortname = re.match('^curvature_coupling_drilldown\.(.+)$',name).group(1)
				label = '%s\n%s%s'%(shortname,table[name]['hypo']['design']['curvature_positions']['method'],
					'\nblurry_explicit' if 
					table[name]['hypo']['design'].get('weighting_scheme',None)=='blurry_explicit' else '')
				ax.bar([nnum],vals,width=1.0,color=colors[nnum],label=label)
			ax.set_title(prop)
			ax.set_xticks([])
			if prop=='vibe': ax.set_yscale('log')
			if prop=='curvature': ax.axhline(0.0,c='k',lw=1)
		legend = ax.legend(title='keys',bbox_to_anchor=(1.05,0.,1.,1.),
			loc='lower left',frameon=False,fontsize=12)	
		picturesave('fig.drilldown_summary.%s.bars'%work.namer.short_namer(sn,spot='sims'),
			work.plotdir,backup=False,version=False,meta={},extras=[legend])

	#---plot spectra with the same collors
	for snum,sn in enumerate(sns):
		axes,fig = square_tiles(1,figsize=(8,8))
		ax = axes[0]
		colors = ncolors(len(keysets[sn]))
		for nnum,name in enumerate(keysets[sn]):
			ratios = results[name]['ratios']
			qs = results[name]['qs']
			ax.scatter(qs,ratios,marker='o',s=10,lw=0,color=colors[nnum],alpha=0.65,label=name)
			high_cutoff = table[name]['hypo']['fitting']['high_cutoff']
			ax.axvline(high_cutoff,color=colors[nnum])
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.axhline(1.0,c='k')
		ax.set_title('%s'%(sn))
		picturesave('fig.drilldown_summary.%s.spectra'%work.namer.short_namer(sn,spot='sims'),
			work.plotdir,backup=False,version=False,meta={},extras=[])

def debugger_careful():
	"""
	Note that the development checkpoint held a longer, previous debugger function.
	"""

	"""
	CAREFUL FITTING STRATEGY	
	executive summary:
		generate an objective for zero curvature
		use the optimizer to solve for kappa, gamma, vibe
		plot the undulation spectra the usual way
		plot the spectra the usual way and compare the results?
		??? find some way to merge these?
		"""

	#---imports and get loaders
	from codes.curvature_coupling.curvature_coupling import prepare_objective,formulate_wavevectors
	from codes.curvature_coupling.curvature_coupling import curvature_sum_function,prepare_residual
	from codes.curvature_coupling.curvature_coupling import blurry_binner,gopher,prepare_oscillator_function
	from codes.curvature_coupling.tools import fft_field
	loader_spec = {'module':'codes.curvature_coupling_loader',
		'function':'curvature_coupling_loader_membrane'}
	loader_func = gopher(loader_spec,module_name='module',variable_name='function',work=work)
	oscillator_function = prepare_oscillator_function()

	#---switching
	probe_vibe = False

	#---settings
	sn = 'membrane-v651-enthx8'
	design_name = 'v6_fixed_extent_all_frames_2'
	design_name = 'v2_fixed_extent_2'
	#---!!!!!!!!!!!!!!!!!!! do we need a separate fitting band?
	high_cutoff = high_cutoff_undulate = 1.0
	low_cutoff = 0.0
	midplane_method = 'flat'
	initial_kappa = 30.0
	initial_vibe = 0.0
	initial_curvature = 0.005
	optimize_method = 'Nelder-Mead'
	binner_method = 'explicit'
	weighting_scheme = 'blurry_explicit'
	ndrops = 8
	trial = 0.0
	narrow = True
	midplane_method = 'flat'
	fit_tension = True
	fit_style = ['band,perfect,curvefit',
		'band,perfect,curvefit-crossover', # decent base case
		'band,perfect,fit', # broke
		'band,blurry,curvefit', # dev
		][-1]
	#---get heights
	if 'hqs' not in globals():
		global hqs
		#---! somewhat redundant with the InvestigateCurvature class
		memory = loader_func(midplane_method=midplane_method,
			data=dict(undulations=data['undulations']))
		hqs = memory[(sn,'hqs')]

	#---assemble required data and functions for the objective
	fft_function = fft_field
	residual_function = prepare_residual()
	curvature_fields = datas[design_name][sn]['fields_unity']
	wavevectors_form = formulate_wavevectors(
		vecs=data['undulations'][sn]['data']['vecs'],
		dims=data['undulations'][sn]['data']['mesh'].shape[-2:])
	wavevectors,area = wavevectors_form['wavevectors'],wavevectors_form['area']
	curvature_sum_function = curvature_sum_function
	#---! make this more systematic
	#---! note that band must match binner_method (which defaults to blurry)
	band = np.where(np.all((wavevectors>=low_cutoff,wavevectors<=high_cutoff),axis=0))[0]
	#---subsample hqs for the 100-frame tests
	if design_name in ['v2_fixed_extent_2']:
		frameslice = np.linspace(0,len(hqs)-1,100).astype(int)
	else: frameslice = slice(None,None)

	#---return the objective function decorated for this optimization
	objective = prepare_objective(
		hqs=hqs[frameslice],curvature_fields=curvature_fields,
		wavevectors=wavevectors,area=area,
		curvature_sum_function=curvature_sum_function,fft_function=fft_field,
		band=band,residual_function=residual_function,blurry_binner=blurry_binner,
		binner_method=binner_method,weighting_scheme=weighting_scheme,
		positive_vibe=True,inner_sign=1.0,
		ndrops_uniform=ndrops,fix_curvature=trial)

	def callback(args):
		"""Watch the optimization."""
		global Nfeval
		print('step %s: %s'%(Nfeval,args))
		Nfeval += 1

	global Nfeval
	Nfeval = 0
	initial_conditions = [initial_kappa,0.0,initial_vibe]
	fit = scipy.optimize.minimize(objective,x0=tuple(initial_conditions),
		callback=callback,method=optimize_method)
	print('fitted: %s'%fit)

	#---plot undulations and energy spectra together
	axes,fig = square_tiles(3 if probe_vibe else 2,figsize=(12,8),favor_rows=True)

	#---plot undulations
	ax = axes[0]

	sn = 'membrane-v651-enthx8'
	surf = data['undulations'][sn]['data']['mesh'].mean(axis=0)
	vecs = data['undulations'][sn]['data']['vecs']

	lims = [0.,high_cutoff_undulate]
	colors = {sn:{'binned':'k','fitted':'r','line':'b'}}

	uspec = calculate_undulations(surf,vecs,fit_style=fit_style,lims=lims,
		midplane_method=midplane_method,fit_tension=fit_tension)
	label = r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'+'\n'\
		r'$\mathrm{\gamma='+('{0:1.2E}'.format(uspec.get('sigma',0.0)))+'\:k_BT {nm}^{-2}}$'
	q_binned,energy_binned = uspec['q_binned'][1:],uspec['energy_binned'][1:]
	ax.plot(uspec['q_raw'],uspec['energy_raw'],'.',lw=0,markersize=10,markeredgewidth=0,
		label=None,alpha=0.2,color=colors[sn]['binned'])

	q_fit,energy_fit = np.transpose(uspec['points'])
	ax.plot(q_fit,energy_fit,'.',lw=0,markersize=4,markeredgewidth=0,
		label=label,alpha=1.,zorder=4,color=colors[sn]['fitted'])
	def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
		return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))
	ax.plot(q_binned,hqhq(q_binned,kappa=uspec['kappa'],sigma=uspec['sigma'],
		area=uspec['area']),lw=1,zorder=3,color=colors[sn]['line'])

	art = {'fs':{'legend':12}}
	add_undulation_labels(ax,art=art)
	add_std_legend(ax,loc='upper right',art=art)
	add_axgrid(ax,art=art)

	#---plot energy spectra
	ax = axes[1]

	#---plot the elastic spectrum
	spectrum = objective(fit.x,mode='elastic')
	label = r'$\mathrm{\kappa='+('%.1f'%fit.x[0])+'\:k_BT}$'+'\n'\
		r'$\mathrm{\gamma='+('{0:1.2E}'.format(fit.x[1]))+'\:k_BT {nm}^{-2}}$'
	ax.scatter(wavevectors,spectrum,marker='o',s=20,lw=0,color='k',zorder=2,alpha=0.1,label=label)
	#---blurry binner
	qbin,ebin = blurry_binner(wavevectors,spectrum)
	bandbin = np.where(np.all((qbin>=low_cutoff,qbin<=high_cutoff),axis=0))[0]
	ax.scatter(qbin[bandbin],ebin[bandbin],marker='o',s=20,lw=0,color='r',zorder=3,alpha=1.0)
	#---plot the mode-freezing expectation
	vibe = fit.x[2]
	frozen = oscillator_function(vibe,wavevectors)

	ax.scatter(wavevectors,frozen,marker='o',s=10,lw=0,color='b',zorder=3,alpha=1.0)
	#---blurry binner
	qbin,ebin = blurry_binner(wavevectors,frozen*spectrum)
	bandbin = np.where(np.all((qbin>=low_cutoff,qbin<=high_cutoff),axis=0))[0]
	ax.scatter(qbin[bandbin],ebin[bandbin],marker='o',s=20,lw=0,color='r',zorder=3,alpha=1.0)
	#---all points, corrected in the background
	ax.scatter(wavevectors,frozen*spectrum,marker='o',s=10,lw=0,color='k',zorder=2,alpha=0.1)

	ax.axvline(high_cutoff,color='k')
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axhline(1.0,c='k',zorder=1)
	if narrow: ax.set_ylim(0.1,10)
	add_std_legend(ax,loc='upper right',art=art)

	#---double back and plot the altered
	axes[0].plot(q_binned,hqhq(q_binned,kappa=fit.x[0],sigma=fit.x[1],
		area=uspec['area']),lw=1,zorder=3,color='r')

	#---repeat with a different vibe
	if probe_vibe:
		ax = axes[2]

		#---plot the elastic spectrum
		spectrum = objective(fit.x,mode='elastic')
		label = r'$\mathrm{\kappa='+('%.1f'%fit.x[0])+'\:k_BT}$'+'\n'\
			r'$\mathrm{\gamma='+('{0:1.2E}'.format(fit.x[1]))+'\:k_BT {nm}^{-2}}$'
		ax.scatter(wavevectors,spectrum,marker='o',s=20,lw=0,color='k',zorder=2,alpha=0.1,label=label)
		#---blurry binner
		qbin,ebin = blurry_binner(wavevectors,spectrum)
		bandbin = np.where(np.all((qbin>=low_cutoff,qbin<=high_cutoff),axis=0))[0]
		ax.scatter(qbin[bandbin],ebin[bandbin],marker='o',s=20,lw=0,color='r',zorder=3,alpha=1.0)
		#---plot the mode-freezing expectation
		vibe = 0.5
		frozen = oscillator_function(vibe,wavevectors)

		ax.scatter(qbin,frozen,marker='o',s=10,lw=0,color='b',zorder=3,alpha=1.0)
		ax.scatter(qbin,ebin*frozen,marker='o',s=8,lw=0,color='c',zorder=4,alpha=1.0)
		ax.scatter(qbin,ebin/frozen,marker='o',s=8,lw=0,color='m',zorder=4,alpha=1.0)

		ax.axvline(high_cutoff,color='k')
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.axhline(1.0,c='k',zorder=1)
		add_std_legend(ax,loc='lower right',art=art)

	picturesave('fig.DEBUG5',work.plotdir,backup=False,version=False,meta={},extras=[])

	#---! testing minimization here but moved to the big sweep instead
	if False:

		#---return the objective function decorated for this optimization
		objective_free = prepare_objective(
			hqs=hqs[frameslice],curvature_fields=curvature_fields,
			wavevectors=wavevectors,area=area,
			curvature_sum_function=curvature_sum_function,fft_function=fft_field,
			band=band,residual_function=residual_function,blurry_binner=blurry_binner,
			binner_method=binner_method,weighting_scheme=weighting_scheme,
			positive_vibe=True,inner_sign=1.0,
			ndrops_uniform=0)

		Nfeval = 0
		#---try a free optimization
		initial_conditions = [fit.x[0],fit.x[1],fit.x[2]]+[initial_curvature for i in range(ndrops)]
		fit = scipy.optimize.minimize(objective_free,x0=tuple(initial_conditions),
			callback=callback,method=optimize_method)
		spectrum = objective_free(fit.x,mode='elastic')
		print(fit.x)

@autoplot(plotrun)
def main(switch=0b0110):
	"""
	Main calculation and plot loop.
	"""
	nswitch = 4
	#---decide what to run using a binary argument
	args_switch = [(int(switch)/i%2)==1 for i in [2**i for i in range(nswitch)][::-1]]
	compute,make_plots,make_summary,run_debugger = args_switch
	global table
	#---use all simulations
	sns = work.sns()
	tags = master_tags[round_name]
	if not os.path.isfile(os.path.join(work.postdir,coordinator)): new_coordinator(coordinator)
	hypos = prepare_hypotheses()
	table = fetch_hypothesis_table(coordinator)
	calcs_new = prepare_calcs(tags=tags,hypos=hypos,sns=sns)
	if not compute: status('skipping calculations',tag='STATUS')
	else:
		#---process optimization requests
		while calcs_new:
			calc_this = calcs_new.pop(0)
			if calc_this in table.values(): status('already computed %s'%calc_this,tag='status')
			#---perform a new optimization
			else:
				#---get the base name
				base_fn = calc_namer(**calc_this)
				#---increment the number
				version_nums = sorted([int(re.match('^%s\.v(\d+)\.dat$'%
					base_fn,os.path.basename(i)).group(1)) 
					for i in glob.glob(os.path.join(work.postdir,'%s*'%base_fn))])
				if version_nums!=range(len(version_nums)):
					raise Exception('version numbering error on %s: %s'%(base_fn,version_nums))
				fn = '%s.v%d'%(base_fn,len(version_nums))
				#---check if this already exists
				if os.path.isfile(os.path.join(work.postdir,'%s.dat'%fn)):
					raise Exception('already computed %s'%fn)
				run_optimizer(fn=fn,**calc_this)
				table[fn] = calc_this
				write_hypothesis_table(fn=coordinator,table=table)
				status('computed %s'%fn,tag='status')
	#---load all of the results
	#---! note that the current execution scheme means that I can work on plots and compute simultaneously
	if make_plots or make_summary:
		global results
		results = dict([(name,get_drilldown(name)) for name in table])
	#---main routes to plotters
	if make_plots:
		#---plot when the computations are ready
		for key in table: 
			status(key,tag='plot')
			plot(key,table)
	#---perform summary plots
	if make_summary: 
		individual_reviews_drilldown_height_summary()
		summary_plots()
	#---debugging
	if run_debugger: debugger_careful()

#---iterative development
if __name__=='__replotting__': pass

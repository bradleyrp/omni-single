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
import scipy.interpolate
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patheffects as path_effects
machine_eps = eps = np.finfo(float).eps
import copy,json,glob,time
import joblib
import h5py

def stringer(x,p=5):
	"""A nice way to print numbers."""
	return (' '*(1*(x>0)))+('%.1e'%x 
		if (abs(x)<10**(-1*p) or abs(x)>10**p) else ('{:>%df}'%(p+2+1*(x<0))).format(x))

if 'round_name' not in globals():
	settings = {}
	#---get instructions
	exec(open('calcs/settings_curvature_undulation_coupling_drilldown.py').read(),settings)
	master_tags,prepare_hypotheses,round_name = [settings[k] for k in 
		'master_tags,prepare_hypotheses,round_name'.split(',')]
	#---where to catalog the completed optimizations
	coordinator = 'curvature_coupling_drilldown.json'
	del settings

#---this script has as very ad hoc namespace
global sn,midplane_method,memory

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
		elif self.scipy_optimize_function=='fmin':
			self.fit = scipy.optimize.fmin(self.objective,x0=tuple(init),disp=True)
			self.fitted = self.fit
		else: raise Exception('unclear optimize procedure')
		status(str(self.fitted),tag='result')
	def callback(self,args):
		"""Watch the optimization."""
		args_string = ' '.join([stringer(a) for a in args])
		output = (u'\r' if self.stepno>0 else '\n')+'[OPTIMIZE] step %s: %s'%(self.stepno,args_string)
		if ~self.silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		self.stepno += 1

def prepare_optimization_components():
	"""
	Import individual components for optimization.
	"""
	if 'loader_func' not in globals():
		global loader_func,prepare_objective,formulate_wavevectors,curvature_sum_function
		global prepare_residual,blurry_binner,gopher,prepare_oscillator_function,log_reverse,fft_field
		global loader_spec,oscillator_function
		#---imports and get loaders
		from codes.curvature_coupling.curvature_coupling import prepare_objective,formulate_wavevectors
		from codes.curvature_coupling.curvature_coupling import curvature_sum_function,prepare_residual
		from codes.curvature_coupling.curvature_coupling import blurry_binner,gopher
		from codes.curvature_coupling.curvature_coupling import prepare_oscillator_function,log_reverse
		from codes.curvature_coupling.tools import fft_field
		loader_spec = {'module':'codes.curvature_coupling_loader',
			'function':'curvature_coupling_loader_membrane'}
		loader_func = gopher(loader_spec,module_name='module',variable_name='function',work=work)
		oscillator_function = prepare_oscillator_function()

def get_heights():
	"""Load heights."""
	global sn,memory
	#---get heights
	if 'hqs' not in globals():
		global hqs
		midplane_method = 'flat'
		#---! somewhat redundant with the InvestigateCurvature class
		memory = loader_func(midplane_method=midplane_method,
			data=dict(undulations=data['undulations']))
		hqs = memory[(sn,'hqs')]

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

	prepare_optimization_components()
	#---switching
	probe_vibe = False

	#---settings
	sn = 'membrane-v651-enthx8'
	design_name = 'v6_fixed_extent_all_frames_2'
	design_name = 'v2_fixed_extent_2'
	#---! should we be fitting in two bands
	high_cutoff = high_cutoff_undulate = 1.0
	low_cutoff = 0.0
	midplane_method = 'flat'
	initial_kappa = 30.0
	initial_vibe = 5.0
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

	#---assemble required data and functions for the objective
	get_heights()
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

def debug_compare_alt():
	"""
	Compare standard and alternate master modes.
	"""

	sn = 'membrane-v651-enthx8'
	prepare_optimization_components()
	get_heights()

	#---digression
	axes,fig = square_tiles(1,figsize=(8,8))
	ax = axes[0]

	#---settings
	low_cutoff = 0.0
	midplane_method = 'flat'
	high_cutoff = high_cutoff_undulate = 1.0
	binner_method = 'explicit'
	weighting_scheme = [None,'explicit_blurry'][0]
	#---! needs removed
	initial_kappa = 20.0
	initial_vibe = 0.0
	#---simulation-specific settings
	ndrops = 8

	#curvature = 0.04
	curvature = 0.0
	
	# via: make plot curvature_undulation_coupling_drilldown meta=calcs/specs/curvature_s4_all_drilldown.yaml
	design_name = 'v6_fixed_extent_all_frames_2'
	frameslice = slice(None,None)
	design_name = 'v2_fixed_extent_2'
	frameslice = np.linspace(0,len(hqs)-1,100).astype(int)

	fft_function = fft_field
	residual_function = prepare_residual()
	wavevectors_form = formulate_wavevectors(
		vecs=data['undulations'][sn]['data']['vecs'],
		dims=data['undulations'][sn]['data']['mesh'].shape[-2:])
	wavevectors,area = wavevectors_form['wavevectors'],wavevectors_form['area']
	band = np.where(np.all((wavevectors>=low_cutoff,wavevectors<=high_cutoff),axis=0))[0]

	if 'opt' not in globals():

		curvature_fields = datas[design_name][sn]['fields_unity']

		objective = prepare_objective(
			hqs=hqs[frameslice],curvature_fields=curvature_fields,
			wavevectors=wavevectors,area=area,
			curvature_sum_function=curvature_sum_function,fft_function=fft_field,
			band=band,residual_function=residual_function,blurry_binner=blurry_binner,
			binner_method=binner_method,weighting_scheme=weighting_scheme,
			#---! removed positive vibe
			positive_vibe=True,
			inner_sign=1.0,ndrops_uniform=ndrops,fix_curvature=curvature,
			oscillator_function=prepare_oscillator_function(reverse=True))
		opt = QuickOpt(objective=objective,scipy_optimize_function='minimize',
			silent=True,init=(initial_kappa,0.0,initial_vibe))
		spectrum = objective(opt.fitted,mode='elastic')

		objective_alt = prepare_objective(
			hqs=hqs[frameslice],curvature_fields=curvature_fields,
			wavevectors=wavevectors,area=area,
			curvature_sum_function=curvature_sum_function,fft_function=fft_field,
			band=band,blurry_binner=blurry_binner,
			residual_function=prepare_residual(mode='alt'),
			binner_method=binner_method,weighting_scheme=weighting_scheme,
			#---! removed positive vibe
			positive_vibe=False,
			inner_sign=-1.0,
			oscillator_function=prepare_oscillator_function(reverse=False,positive_vibe=False),
			ndrops_uniform=ndrops,fix_curvature=curvature)
		objective_alt((initial_kappa,0.0,5),debug=True)
		opt_alt = QuickOpt(objective=objective_alt,scipy_optimize_function='fmin',
			silent=True,init=(initial_kappa,0.0,5))
		spectrum_alt = objective_alt(opt_alt.fitted,mode='elastic')

	vibe = opt.fitted[2]
	frozen = oscillator_function(vibe,wavevectors)
	subject = spectrum*frozen
	frozen_alt = prepare_oscillator_function(
		positive_vibe=False)(opt_alt.fitted[2],wavevectors)
	subject_alt = spectrum_alt/frozen_alt
	if 0: ax.plot(wavevectors,subject,
		markersize=5,marker='o',lw=0,color='magenta',zorder=2,alpha=1.0,
		markeredgewidth=1.0,markeredgecolor='k')
	xp,yp = perfect_collapser(wavevectors,subject)
	if 1: ax.plot(xp,yp,
		lw=2,zorder=10,color='magenta',solid_capstyle='round',
		path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])
	xpf,ypf = perfect_collapser(wavevectors,frozen)
	ax.plot(xpf,log_reverse(ypf),'--',c='magenta',lw=1)
	xpf,ypf = perfect_collapser(wavevectors,spectrum)
	ax.plot(xpf,ypf,'-',c='k',lw=1,alpha=0.6,zorder=1)
	if 1: 
		if 0: ax.plot(wavevectors,subject_alt,
			markersize=5,marker='o',lw=0,color='cyan',zorder=2,alpha=1.0,
			markeredgewidth=1.0,markeredgecolor='k')
		xp,yp = perfect_collapser(wavevectors,subject_alt)
		if 1: ax.plot(xp,yp,
			lw=2,zorder=10,color='cyan',solid_capstyle='round',
			path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])
		xpf,ypf = perfect_collapser(wavevectors,frozen_alt)
		ax.plot(xpf,ypf,'--',c='cyan',lw=1)
		xpf,ypf = perfect_collapser(wavevectors,spectrum_alt)
		ax.plot(xpf,ypf,'-',c='k',lw=1,alpha=0.6,zorder=1)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axhline(1.0,c='k')
	ax.set_ylim(0.1,10)
	ax.axvline(high_cutoff,color='k')
	print(opt_alt.fitted)
	objective_alt(opt_alt.fitted,debug=True,mode='elastic')
	picturesave('fig.DEBUG7c',work.plotdir,backup=False,version=False,meta={})

def debug_resurvey(sn,figname='DEBUG11',use_protrusion=True):
	"""
	Run the survey with the alternate master mode.
	"""
	prepare_optimization_components()
	get_heights()

	#---define a landscape
	extent_keys = ['v3_fixed_extent_0.5','v1_fixed_extent_1','v2_fixed_extent_2',
		'v5_fixed_extent_3','v4_fixed_extent_4','v6_fixed_extent_5','v7_fixed_extent_6',
		'v8_fixed_extent_8','v9_fixed_extent_10']
	extents = np.array([datas[k][sn]['spec']['extents']['extent'] for k in extent_keys])
	curvatures = np.arange(-0.1,0.1+0.005,0.01)
	#---moar
	halfsweep = np.concatenate((np.arange(0.01,0.1+0.01,0.01),np.arange(0.2,1.0+0.1,0.1)))
	curvatures = np.concatenate((halfsweep[::-1]*-1,[0.],halfsweep))

	#---settings
	low_cutoff = 0.0
	midplane_method = 'flat'
	high_cutoff = high_cutoff_undulate = 2.0
	binner_method = 'explicit'
	weighting_scheme = [None,'blurry_explicit'][-1]
	#---! needs removed
	initial_kappa = 20.0
	initial_vibe = 0.0
	#---simulation-specific settings
	ndrops = 8

	#---prepare objective functions
	if 'job_toc' not in globals():
		job_toc = {}

		fft_function = fft_field
		residual_function = prepare_residual()
		wavevectors_form = formulate_wavevectors(
			vecs=data['undulations'][sn]['data']['vecs'],
			dims=data['undulations'][sn]['data']['mesh'].shape[-2:])
		wavevectors,area = wavevectors_form['wavevectors'],wavevectors_form['area']
		band = np.where(np.all((wavevectors>=low_cutoff,wavevectors<=high_cutoff),axis=0))[0]

		start = time.time()
		#---prepare objective functions for each extent
		#---! cannot pickle function objects so cannot do this in parallel
		for ii,(design_name,extent) in enumerate(zip(extent_keys,extents)):
			for jj,curvature in enumerate(curvatures):
				print('\n')
				status('preparing objective function',i=ii*len(curvatures)+jj,
					looplen=len(extents)*len(curvatures),start=start,tag='build')
				curvature_fields = datas[design_name][sn]['fields_unity']
				#---! hard coded
				frameslice = np.linspace(0,len(hqs)-1,100).astype(int)
				if use_protrusion: initial = (initial_kappa,0.0,0.0,initial_vibe) 
				else: initial = (initial_kappa,0.0,initial_vibe)
				objective = prepare_objective(
					#---! note that master mode was added after this function was working. check overrides!
					master_mode='standard',
					hqs=hqs[frameslice],curvature_fields=curvature_fields,
					wavevectors=wavevectors,area=area,
					curvature_sum_function=curvature_sum_function,fft_function=fft_field,
					band=band,residual_function=residual_function,blurry_binner=blurry_binner,
					binner_method=binner_method,weighting_scheme=weighting_scheme,
					positive_vibe=True,inner_sign=1.0,
					ndrops_uniform=ndrops,fix_curvature=curvature)
				opt = QuickOpt(objective=objective,silent=True,init=initial)
				job_toc[(extent,curvature)] = opt

	figsize = (12,10)
	contour_interp_pts = 100
	contour_line_skip = 4
	contour_nlevels = 100
	under_color = 'm'
	
	axes,fig = square_tiles(2,figsize,favor_rows=True)
	ax = axes[0]
	raw = np.array([[job_toc[(e,c)].fit.fun for e in extents] for c in curvatures])
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
	error_min = raw.min()*0.1
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
	picturesave('fig.%s'%figname,work.plotdir,backup=False,version=False,meta={})
	i,j = np.unravel_index(raw.argmin(),raw.shape)
	print((curvatures[i],extents[j]))
	print(curvature_extent_error[curvature_extent_error[:,2].argmin()])

def debug_search_alt(sn=None,master_mode=None,subtractor=None):
	"""
	Search using the alternate master mode.
	"""
	
	if sn==None: sn = ['membrane-v651-enthx8','membrane-v1005'][-1]
	if master_mode==None: master_mode = ['standard','alt','comp'][0]
	if master_mode=='subtractor' and type(subtractor)==type(None): raise Exception('need subtractor')

	prepare_optimization_components()
	get_heights()

	#---settings
	midplane_method = 'flat'
	low_cutoff = 0.0
	high_cutoff = high_cutoff_undulate = 1.0
	binner_method = 'explicit'
	weighting_scheme = [None,'explicit_blurry'][0]
	#---! needs removed
	initial_kappa = 20.0
	initial_vibe = 0.0
	#---simulation-specific settings
	ndrops = work.meta[sn]['nprots']

	use_all_frames = False
	if use_all_frames:
		# call with: meta=calcs/specs/curvature_s4_all_drilldown.yaml
		design_name = 'v6_fixed_extent_all_frames_2'
		frameslice = slice(None,None)
	else: 
		design_name = 'v2_fixed_extent_2'
		frameslice = np.linspace(0,len(hqs)-1,100).astype(int)

	fft_function = fft_field
	wavevectors_form = formulate_wavevectors(
		vecs=data['undulations'][sn]['data']['vecs'],
		dims=data['undulations'][sn]['data']['mesh'].shape[-2:])
	wavevectors,area = wavevectors_form['wavevectors'],wavevectors_form['area']
	band = np.where(np.all((wavevectors>=low_cutoff,wavevectors<=high_cutoff),axis=0))[0]
	curvature_fields = datas[design_name][sn]['fields_unity']

	if 'opt' not in globals():
		objective = prepare_objective(
			master_mode=master_mode,
			hqs=hqs[frameslice],curvature_fields=curvature_fields,
			wavevectors=wavevectors,area=area,
			curvature_sum_function=curvature_sum_function,fft_function=fft_field,
			band=band,blurry_binner=blurry_binner,
			binner_method=binner_method,weighting_scheme=weighting_scheme,
			ndrops_uniform=ndrops,fix_curvature=0.,
			**({'subtractor':subtractor} if master_mode=='subtractor' else {}))
		if master_mode in ['standard','alt']: initial = (initial_kappa,0.,5.)
		elif master_mode in ['comp']: initial = (initial_kappa,0.0,0.0,5.,0.)
		elif master_mode=='subtractor': initial = (initial_kappa,0.0,0.)
		opt = QuickOpt(objective=objective,silent=False,init=initial)
		spectrum = objective(opt.fitted,mode='elastic')

	import ipdb;ipdb.set_trace()

	#---digression
	axes,fig = square_tiles(1,figsize=(8,8))
	ax = axes[0]
	if master_mode=='alt':
		frozen = prepare_oscillator_function(reverse=False,positive_vibe=False)(
			opt.fitted[2],wavevectors)
		subject = spectrum/frozen
	elif master_mode in ['standard']:
		frozen = prepare_oscillator_function(reverse=True,positive_vibe=True)(
			opt.fitted[2],wavevectors)
		subject = spectrum/frozen
	elif master_mode=='comp':
		frozen = prepare_oscillator_function(reverse=True,positive_vibe=True)(
			opt.fitted[3],wavevectors)
		subject = spectrum/frozen
	elif master_mode=='subtractor':
		frozen = subtractor + 1.0
		subject = spectrum*(subtractor+1.0)
	if 0: ax.plot(wavevectors,subject,
		markersize=5,marker='o',lw=0,color='magenta',zorder=2,alpha=1.0,
		markeredgewidth=1.0,markeredgecolor='k')
	xp,yp = perfect_collapser(wavevectors,subject)
	if 1: ax.plot(xp,yp,
		lw=2,zorder=10,color='magenta',solid_capstyle='round',
		path_effects=[path_effects.withStroke(linewidth=4,foreground='k')])
	xpf,ypf = perfect_collapser(wavevectors,frozen)
	if master_mode=='subtractor': ax.plot(xpf,ypf,'--',c='magenta',lw=1)
	else: ax.plot(xpf,log_reverse(ypf),'--',c='magenta',lw=1)
	xpf,ypf = perfect_collapser(wavevectors,spectrum)
	ax.plot(xpf,ypf,'-',c='k',lw=1,alpha=0.6,zorder=1)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.axhline(1.0,c='k')
	#ax.set_ylim(0.1,10)
	ax.axvline(high_cutoff,color='k')
	picturesave('fig.DEBUG.search_master_mode_%s.%s%s'%(master_mode,sn,
		'.all_frames' if use_all_frames else ''),
		work.plotdir,backup=False,version=False,meta={})

def debug_protrusion(return_subtractor=False):
	"""
	CRAZY DETAILED FIT METHOD AND AUDIT.
	IT ENDS HERE! THIS IS THE RECKONING!
	"""
	prepare_optimization_components()
	get_heights()
	surf = data['undulations'][sn]['data']['mesh'].mean(axis=0)
	vecs = data['undulations'][sn]['data']['vecs']

	#---settings
	q_min,q_tension_bending_hard,q_bending_tension_hard,q_cut = 0.0,0.3,0.6,2.0
	midplane_method = 'flat'
	#---! settings almost gone
	fit_style = 'band,perfect,curvefit-crossover'
	fit_tension = True

	figsize = (12,8)
	axes,fig = square_tiles(2,figsize,favor_rows=True)
	ax = axes[0]

	#---! original
	if False:
		colors = {sn:{'binned':'k','fitted':'r','line':'b'}}
		uspec = calculate_undulations(surf,vecs,fit_style=fit_style,lims=(q_min,q_cut),
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

	wavevectors_form = formulate_wavevectors(
		vecs=data['undulations'][sn]['data']['vecs'],
		dims=data['undulations'][sn]['data']['mesh'].shape[-2:])
	qs,area = wavevectors,area = wavevectors_form['wavevectors'],wavevectors_form['area']
	hqhq = (np.abs(hqs).reshape((len(hqs),-1)).mean(axis=0)[1:])**2 #---! check this
	band = np.all((qs>=q_min,qs<=q_cut),axis=0)

	def model(q_raw,sigma,kappa,gamma_p):
		kappa,gamma_p = np.abs(kappa),np.abs(gamma_p)
		if sigma>1.0: sigma = 1.0
		elif sigma<-1.0: sigma = -1.0
		######### 
		sigma = 0.0
		return ((1.0)/(area))*(
			(1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps)+(1.)/(gamma_p*q_raw**2+machine_eps))
	def residual(a,b): return (np.log10(a/b)**2).mean()
	def objective((sigma,kappa,gamma_p)): return residual(hqhq[band],model(qs[band],sigma,kappa,gamma_p))
	#---adding harmonic oscillator correction
	def model(q_raw,sigma,kappa,gamma_p,vibe):
		kappa,gamma_p = np.abs(kappa),np.abs(gamma_p)
		if sigma>1.0: sigma = 1.0
		elif sigma<-1.0: sigma = -1.0
		######### 
		sigma = 0.0
		pure = ((1.0)/(area))*(
			(1.)/(kappa*q_raw**4+sigma*q_raw**2+machine_eps)+(1.)/(gamma_p*q_raw**2+machine_eps))
		osc = (vibe*q_raw+machine_eps)*(1./(np.exp(vibe*q_raw)-1)+machine_eps)
		return pure*osc
	def residual(a,b): return (np.log10(a/b)**2).mean()
	def objective((sigma,kappa,gamma_p,vibe)): 
		return residual(hqhq[band],model(qs[band],sigma,kappa,gamma_p,vibe))

	global stepno
	stepno = 1
	def callback(args,silent=False):
		"""Watch the optimization."""
		global stepno
		args_string = ' '.join([stringer(a) for a in args])
		output = (u'\r' if stepno>0 else '\n')+'[OPTIMIZE] step %s: %s'%(stepno,args_string)
		if ~silent:
			sys.stdout.flush()
			sys.stdout.write(output)
		stepno += 1

	initial_conditions = (0.,20.,0.,0.)
	fit = scipy.optimize.minimize(objective,
		x0=initial_conditions,method='Nelder-Mead',callback=callback)

	def hqhqf(q_raw,kappa,exponent=4.0): return 1.0/(area/1.0*(kappa*q_raw**(exponent)))

	#----sort the wavevectors
	qsort = np.argsort(qs)
	fitted = model(qs,*fit.x)
	ax.plot(qs[qsort],fitted[qsort],'-',color='c',zorder=3)
	ax.plot(qs[qsort],hqhq[qsort],'.',color='m',zorder=1)
	ax.plot(qs[qsort],hqhqf(qs[qsort],fit.x[0]),'.',color='g',zorder=1)
	ax.set_xscale('log')
	ax.set_yscale('log')

	ax = axes[1]
	if 0: ax.plot(qs[qsort],hqhq[qsort]*qs**4,'.',color='m',zorder=1)
	#---show the deviation correctly
	ax.plot(qs[qsort],fitted[qsort]*qs[qsort]**4,'.',color='m',zorder=1)
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_ylim(0.1,10.)
	picturesave('fig.RECKON.%s'%(sn),work.plotdir,backup=False,version=False,meta={})
	if return_subtractor: return fitted

@autoplot(plotrun)
def main(switch=0b0000):
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
	if not compute: status('skipping calculations',tag='STATUS')
	else:
		calcs_new = prepare_calcs(tags=tags,hypos=hypos,sns=sns)
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
	#---! bleeding edge development
	compare_alt,resurvey,search_alt,protrusion = 0,0,0,0
	if compare_alt: debug_compare_alt()
	if resurvey: debug_resurvey()
	if search_alt: debug_search_alt()
	if protrusion: debug_protrusion()

#---iterative development
if __name__=='__replotting__':

	sn = ['membrane-v651-enthx8','membrane-v1005'][0]
	#---! history of different things I've worked on 
	if False:
		debug_resurvey(sn=sn,figname='DEBUG665',use_protrusion=False)
		debug_resurvey(sn=sn,figname='DEBUG666',use_protrusion=True)
		debug_search_alt(sn=sn,mode='comp')
		debug_protrusion()
		subtractor = debug_protrusion(return_subtractor=True)
		subtractor = 0.*np.ones(subtractor.shape[0])
		debug_search_alt(sn=sn,master_mode='subtractor',subtractor=subtractor)
		debug_protrusion()
		debug_search_alt(sn=sn,master_mode='standard')
	debug_compare_alt()
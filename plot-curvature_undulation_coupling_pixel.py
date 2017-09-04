#!/usr/bin/env python

"""
Plot all of the upstream loops for the curvature undulation coupling, pixel method.
"""

import copy

#---automatically supervise important globals
variables = ['data','datas','printers','routine']
for key in [v for v in variables if v not in globals()]: globals()[key] = None
#---redundant assignment for clarity (plotname is automagically loaded)
plotname = 'curvature_undulation_coupling_pixel'

###---FUNCTIONS

def register_printer(func):
	"""Add decorated functions to a list of "printers" which are the default routine."""
	global printers
	if printers is None: printers = []
	printers.append(func.__name__)
	return func

def collect_upstream_calculations_over_loop():
	"""
	Some plotting and analysis benefits from checking all calculations in an upstream loop (which is 
	contrary to the original design of )
	!!! This script is a candidate for inclusion in omnicalc.py
	"""
	global plotname,work
	#---! move to a separate function
	plotspecs = work.plots.get(plotname,work.calcs.get(plotname,{})).get('specs',{})
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
	if plotname not in work.plots: work.plots[plotname] = copy.deepcopy(work.calcs[calcname])
	#---load other upstream data
	#---get all upstream curvature sweeps
	upstreams,upstreams_stubs = work.calc_meta.unroll_loops(work.calcs[calcname],return_stubs=True)
	#for u in upstreams_stubs: u['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for unum,upstream in enumerate(upstreams_stubs):
		#---use the whittle option to select a particular calculation
		dat,cal = plotload(calcname,whittle_calc={calcname:upstream['specs']})
		tag = upstreams_stubs[unum]['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
	#---singluar means the typical "focus" of the upstream calculation, plural is everything else
	return dict(datas=datas,calcs=calcs,data=data,calc=calc)

@register_printer
def curvature_field_review(extrema_compare=True,fn_base='fig.curvature_pixel.review'):
	"""
	Review the curvature fields in a single plot.
	"""
	global datas,work
	figsize = 12
	cmap_name = 'RdBu_r'
	ntiles = sum([len(v) for v in datas.values()])
	sns = work.sns()
	cf_key = 'cf_first'
	axes,fig = square_tiles(ntiles=ntiles,figsize=figsize,hspace=0.25,wspace=0.25)
	cmax = max([np.abs(d[sn]['cf']).max() for d in datas.values() for sn in work.sns()])
	tags = sorted(datas.keys())
	for dnum,tag in enumerate(tags):
		dat = datas[tag]
		for snum,sn in enumerate(sns):
			ax = axes[dnum*len(sns)+snum]
			#---!!! ULTRAHACK this needs removed
			if 'datau' not in globals():
				global datau
				datau,_ = plotload('undulations')
			vecs = datau[sn]['data']['vecs'].mean(axis=0)
			if not extrema_compare: cmax_this = np.abs(datas[tag][sn]['cf']).max()
			else: cmax_this = cmax
			ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
				vmax=cmax_this,vmin=-1*cmax_this,cmap=mpl.cm.__dict__[cmap_name],
				extent=[0,vecs[0],0,vecs[0]])
			ax.scatter(*datas[tag][sn]['drop_gaussians_points'].mean(axis=0).T,
				c='k',s=1)
			error = datas[tag][sn]['bundle'][sn]['fun']
			extrema = np.abs(datas[tag][sn]['cf']).max()
			ax.set_title('%s\n%s: error=%.4f, extrema=%.4f'%(
				work.meta.get(sn,{}).get('label',sn),tag,error,extrema))
	picturesave(fn_base,work.plotdir,backup=False,version=True,meta={})

@register_printer
def curvature_field_review_self():
	"""Wrapper for self-comparison."""
	curvature_field_review(fn_base='fig.curvature_pixel.review_self',extrema_compare=False)

###---MAIN

def loader():
	"""Load the data only once."""
	#---you must declare all globals (despite the automation)
	global variables,work,datas,data,routine
	#---reload if not all of the globals in the variables
	if any([globals()[v] is None for v in variables]):
		status('loading upstream data',tag='load')
		#---get the routine for this plot (not updated before replot if you change the metadata)
		#---! could you make this automatically detectable
		#avail = ['curvature_field_review','curvature_field_review']
		plotspecs = work.plots.get(plotname,work.calcs.get(plotname,{})).get('specs',{})
		routine = plotspecs.get('routine',printers)
		combodat = collect_upstream_calculations_over_loop()
		data,datas = combodat['data'],combodat['datas']
	else: status('data are already loaded',tag='load')

def printer():
	"""Print the plots automatically from the routine."""
	global routine
	#---! get routine from the metadata here
	if routine is None: routine = list(printers)	
	#---routine items are function names
	for key in routine: 
		status('running routine %s'%key,tag='printer')
		globals()[key]()

if __name__=='__main__': 
	loader()
	printer()

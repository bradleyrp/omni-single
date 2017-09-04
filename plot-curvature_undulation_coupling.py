#!/usr/bin/env python

"""
Plot all of the upstream loops for the curvature undulation coupling analysis.
!!Ryan ported this script quickly. It needs cleaned, checked, pushed, and tested on ocean. Otherwise suitable.
"""

import copy
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from render.wavevids import plothull
str_types = [str,unicode] if sys.version_info<(3,0) else [str]

#---automatically supervise important globals
variables = ['data','datas','printers','routine','postdat','undulations_name','calcs']
for key in [v for v in variables if v not in globals()]: globals()[key] = None
#---redundant assignment for clarity (plotname is automagically loaded)
plotname = 'curvature_undulation_coupling'

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
	datas,calcs = {},{}
	#---loop over upstream calculations and load each specifically, using plotload with whittle_calc
	for unum,upstream in enumerate(upstreams_stubs):
		#---use the whittle option to select a particular calculation
		dat,cal = plotload(calcname,whittle_calc={calcname:upstream['specs']})
		tag = upstreams_stubs[unum]['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
	#---singluar means the typical "focus" of the upstream calculation, plural is everything else
	return dict(datas=datas,calcs=calcs,data=data,calc=calc)

def plot_hull_and_trial_centers(data,sn,ax,n_instances=None,debug_frame=0,color=None):
	"""Plot the protein hull along with the positions of the trial functions."""
	global postdat
	#---! currently only set for a single protein
	trials = data[sn]['drop_gaussians_points'].transpose(1,0,2)
	nframes = len(trials)
	if n_instances!=None: samples = np.arange(0,nframes,nframes/n_instances)
	else: samples = np.array([debug_frame])
	pts = np.concatenate(trials[samples])
	for ptsnum,pts in enumerate(trials[samples]):
		color_this = mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)) if not color else color
		ax.scatter(*pts.T,s=1,c=color_this)
	griddims = data[sn]['cf'].shape
	vecs,points_protein = [postdat[sn][i] for i in ['vecs','points_protein']]
	for ptsnum,pts in enumerate(points_protein[samples][...,:2]):
		color_this = mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)) if not color else color
		plothull(ax,pts,griddims=griddims,vecs=vecs,lw=0,c=color_this)
	ax.set_xlim((0,vecs[0]))
	ax.set_ylim((0,vecs[1]))
	ax.set_aspect(1.0)

@register_printer
def individual_reviews():
	"""
	Loop over all upstream curvature-undulation coupling calculations and plot a panel of review plots.
	"""
	for tag in datas:
		for sn in work.sns():
			status('reviewing curvature for %s'%sn,tag='plot')
			figsize = (10,10)
			cmap_name = 'RdBu_r'	
			#---this plot is an assortment of the following views
			viewnames = ['average_height','average_height_pbc','neighborhood_static',
				'neighborhood_dynamic','average_field','example_field','example_field_pbc',
				'spectrum','spectrum_zoom']
			axes,fig = square_tiles(len(viewnames),figsize)
			#---several plots use the same data
			vecs = postdat[sn]['vecs']
			griddims = datas[tag][sn]['cf'].shape
			#---PLOT the mean curvature field
			ax = axes[viewnames.index('average_field')]
			cmax = np.abs(datas[tag][sn]['cf']).max()
			cmax_instant = np.abs(datas[tag][sn]['cf_first']).max()
			ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
			mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
			ax.scatter(*mean_trial.T,s=1,c='k')
			ax.set_title(r'$\mathrm{\langle C_0(x,y) \rangle}$'+' (max %.3f)'%cmax)
			#---PLOT a single instance of the neighborhood (good for debugging)
			ax = axes[viewnames.index('neighborhood_static')]
			debug_frame = 0
			plot_hull_and_trial_centers(datas[tag],sn,ax,debug_frame=debug_frame,color='k')
			ax.set_title('neighborhood, static, frame %d'%debug_frame)
			#---PLOT a composite of several frames of the neighborhood to see the dynamics
			ax = axes[viewnames.index('neighborhood_dynamic')]
			plot_hull_and_trial_centers(datas[tag],sn,ax,n_instances=10)
			ax.set_title('neighborhood, dynamic')
			#---PLOT the average height
			mesh = data[undulations_name][sn]['data']['mesh']
			surf = mesh.mean(axis=0).mean(axis=0)
			surf -= surf.mean()
			hmax = np.abs(surf).max()
			ax = axes[viewnames.index('average_height')]
			ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
				extent=[0,vecs[0],0,vecs[1]])
			try: mean_prot_pts = postdat[sn]['points_protein_mean'][:,:2]
			except: mean_prot_pts = data[protein_abstractor_name][sn]['data']['points'].mean(axis=0)[:,:2]
			plothull(ax,[mean_prot_pts],griddims=datas[tag][sn]['cf'].shape,vecs=vecs,c='k',lw=0)
			ax.set_title(r'$\mathrm{\langle z \rangle}$'+'(max %.3f)'%hmax)
			#---PLOT an example curvature field ('first' is hard-coded, instead of saving each frame)
			ax = axes[viewnames.index('example_field')]
			example_frame = 0
			cf_first = datas[tag][sn]['cf_first']
			ax.imshow(cf_first.T,origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
			#---there might be some points that are off the chart here? worth checking
			ax.set_xlim((0,vecs[0]))
			ax.set_ylim((0,vecs[1]))
			mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2)[example_frame]
			ax.scatter(*mean_trial.T,s=1,c='k')
			ax.set_title('example field (max %.3f)'%cmax_instant,fontsize=10)
			#---PLOT periodic view of the example field
			ax = axes[viewnames.index('example_field_pbc')]
			ax.imshow(np.tile(cf_first.T,(3,3)),origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],
				extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
			ax.set_title('example field (max %.3f)'%cmax_instant)
			#---PLOT periodic view of the average height
			ax = axes[viewnames.index('average_height_pbc')]
			ax.imshow(np.tile(surf.T,(3,3)),origin='lower',interpolation='nearest',
				vmax=hmax,vmin=-1*hmax,cmap=mpl.cm.__dict__[cmap_name],
				extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
			ax.set_title(r'$\mathrm{\langle z \rangle}$'+'(max %.3f)'%hmax)
			#---PLOT spectrum
			ax = axes[viewnames.index('spectrum')]
			ax.scatter(datas[tag][sn]['qs'],datas[tag][sn]['ratios'],s=4,c='k',alpha=0.25)
			#---! high cutoff is hard-coded here but needs to be removed to the yaml. we need to get default
			hicut = calcs[tag][sn]['calcs']['specs'].get('fitting',{}).get('high_cutoff',1.0)
			qs = datas[tag][sn]['qs']
			band = qs<=hicut
			ax.scatter(datas[tag][sn]['qs'][band],datas[tag][sn]['ratios'][band],s=10,c='k',alpha=1.0)
			ax.axhline(1.0,c='k')
			ax.set_xscale('log')
			ax.set_yscale('log')
			error = datas[tag][sn]['bundle'][sn]['fun']
			ax.set_title('full spectrum')
			ax.grid(True)
			ax.axvline(hicut,ymin=0.0,ymax=1.0,c='k',lw=2)
			#---PLOT spectrum, relevant section
			ax = axes[viewnames.index('spectrum_zoom')]
			qs = datas[tag][sn]['qs']
			band = qs<=hicut
			ys = datas[tag][sn]['ratios'][band]
			ax.scatter(qs[band],ys,s=20,c='k',alpha=1.0,clip_on=False)
			ax.axhline(1.0,c='k')
			ax.set_xscale('log')
			ax.set_yscale('log',subsy=[])
			ax.set_yticks([min(ys),1.0,max(ys)])
			#---intransigent ticks!
			ax.set_yticklabels([('%.2f'%i if type(i)!=bool else '') for i in [min(ys),False,max(ys)]])
			ax.set_xlim(min(qs),hicut)
			error = datas[tag][sn]['bundle'][sn]['fun']
			ax.set_title('spectrum (error %.5f)'%error)
			#---no tick marks on anything
			for ax in axes:
				ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
				ax.tick_params(axis='y',which='both',left='off',right='off',labelbottom='off')
			#---the metadata for this plot comes from the design section
			try: meta = calcs[tag][sn]['calcs']['specs']['specs']
			#---! custom upstream calculations for e.g. dextran project put the specs one level down
			except: meta = calcs[tag][sn]['calcs']['specs']
			#---add high cutoff (from fitting parameters if defined) to the meta
			meta['high_cutoff'] = hicut
			picturesave('fig.coupling_review.%s'%sn,work.plotdir,backup=False,version=True,meta=meta)

###---MAIN

def loader():
	"""Load the data only once."""
	#---you must declare all globals here and above in the "automatic supervise" section
	global variables,work,datas,data,routine,postdat,undulations_name,calcs
	#---reload if not all of the globals in the variables
	if any([globals()[v] is None for v in variables]):
		status('loading upstream data',tag='load')
		#---load sequence from the previous version of this plot script
		plotspecs = work.plots[plotname].get('specs',{})
		routine = plotspecs.get('routine',printers)
		calcname = plotspecs.get('calcname',plotname)
		protein_abstractor_name = plotspecs.get('protein_abstractor_name','protein_abstractor')
		undulations_name = plotspecs.get('undulations_name','undulations')
		#---new method for getting all upstream calculations in the loop
		combodat = collect_upstream_calculations_over_loop()
		data,datas,calcs = combodat['data'],combodat['datas'],combodat['calcs']
		#---extra loading compared to the pixel method from which this was derived, in order to use the new style
		#---...plotting scheme with register_printer
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

#---load and print
if __name__=='__main__': 
	loader()
	printer()

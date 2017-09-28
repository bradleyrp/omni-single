#!/usr/bin/env python

"""
Plot all of the upstream loops for the curvature undulation coupling analysis.
"""

import copy
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from render.wavevids import plothull
from codes.looptools import basic_compute_loop

#---function names to plot or None for all
plot_super.routine = None

@autoload(plot_super)
def loader():
	"""Load data."""
	#---only load once
	if 'data' not in globals():
		global data,datas,routine,postdat,undulations_name,calcs,protein_abstractor_name
		plotspecs = work.plots[plotname].get('specs',{})
		calcname = plotspecs.get('calcname',plotname)
		#---new method for getting all upstream calculations in the loop
		combodat = work.collect_upstream_calculations_over_loop(plotname)
		data,datas,calcs = combodat['data'],combodat['datas'],combodat['calcs']
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

def individual_reviews_plotter(viewnames,out_fn,figsize=(10,10)):
	"""
	Loop over all upstream curvature-undulation coupling calculations and plot a panel of review plots.
	"""
	for tag in datas:
		for sn in work.sns():
			status('reviewing curvature for %s'%sn,tag='plot')
			cmap_name = 'RdBu_r'	
			axes,fig = square_tiles(len(viewnames),figsize,hspace=0.4,wspace=0.4)
			#---several plots use the same data
			vecs = postdat[sn]['vecs']
			griddims = datas[tag][sn]['cf'].shape
			#---shared variables for several plots
			cmax = np.abs(datas[tag][sn]['cf']).max()
			cmax_instant = np.abs(datas[tag][sn]['cf_first']).max()
			hicut = calcs[tag][sn]['calcs']['specs'].get('fitting',{}).get('high_cutoff',1.0)
			#---PLOT the mean curvature field
			if 'average_field' in viewnames:
				ax = axes[viewnames.index('average_field')]
				ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
				ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title(r'$\mathrm{\langle C_0(x,y) \rangle}$'+' (max %.3f)'%cmax)
			#---PLOT a single instance of the neighborhood (good for debugging)
			if 'neighborhood_static' in viewnames:
				ax = axes[viewnames.index('neighborhood_static')]
				debug_frame = 0
				plot_hull_and_trial_centers(datas[tag],sn,ax,debug_frame=debug_frame,color='k')
				ax.set_title('neighborhood, static, frame %d'%debug_frame)
			#---PLOT a composite of several frames of the neighborhood to see the dynamics
			if 'neighborhood_dynamic' in viewnames:
				ax = axes[viewnames.index('neighborhood_dynamic')]
				plot_hull_and_trial_centers(datas[tag],sn,ax,n_instances=10)
				ax.set_title('neighborhood, dynamic')
			#---PLOT the average height
			if 'average_height' in viewnames:
				mesh = data[undulations_name][sn]['data']['mesh']
				surf = mesh.mean(axis=0).mean(axis=0)
				surf -= surf.mean()
				hmax = np.abs(surf).max()
				ax = axes[viewnames.index('average_height')]
				im = ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
					extent=[0,vecs[0],0,vecs[1]],vmax=hmax,vmin=-1*hmax)
				try: mean_prot_pts = postdat[sn]['points_protein_mean'][:,:2]
				except: mean_prot_pts = data[protein_abstractor_name][sn]['data']['points'].mean(axis=0)[:,:2]
				plothull(ax,[mean_prot_pts],griddims=datas[tag][sn]['cf'].shape,vecs=vecs,c='k',lw=0)
				ax.set_title('height profile')
				ax.set_xlabel('x (nm)')
				ax.set_ylabel('y (nm)')
				#---colorbar
				axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.05,0.,1.,1.),
					bbox_transform=ax.transAxes,borderpad=0)
				cbar = plt.colorbar(im,cax=axins,orientation="vertical")
				axins.set_ylabel(r'$\mathrm{\langle z \rangle \, (nm)}$',rotation=270,labelpad=20)
				axins.tick_params(axis='y',left='off',right='off',labelright='on')
			#---PLOT an example curvature field ('first' is hard-coded, instead of saving each frame)
			if 'example_field' in viewnames:
				ax = axes[viewnames.index('example_field')]
				example_frame = 0
				cf_first = datas[tag][sn]['cf_first']
				im = ax.imshow(cf_first.T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				#---colorbar
				axins = inset_axes(ax,width="5%",height="100%",loc=3,bbox_to_anchor=(1.05,0.,1.,1.),
					bbox_transform=ax.transAxes,borderpad=0)
				cbar = plt.colorbar(im,cax=axins,orientation="vertical")
				axins.set_ylabel(r'$\mathrm{C_0\,({nm}^{-1})}$',rotation=270,labelpad=20)
				axins.tick_params(axis='y',left='off',right='off',labelright='on')
				#---there might be some points that are off the chart here? worth checking
				ax.set_xlim((0,vecs[0]))
				ax.set_ylim((0,vecs[1]))
				ax.set_xlabel('x (nm)')
				ax.set_ylabel('y (nm)')
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2)[example_frame]
				ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title('curvature field',fontsize=10)
			#---PLOT periodic view of the example field
			if 'example_field_pbc' in viewnames:
				ax = axes[viewnames.index('example_field_pbc')]
				ax.imshow(np.tile(cf_first.T,(3,3)),origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],
					extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
				ax.set_title('example field (max %.3f)'%cmax_instant)
			#---PLOT periodic view of the average height
			if 'average_height_pbc' in viewnames:
				ax = axes[viewnames.index('average_height_pbc')]
				ax.imshow(np.tile(surf.T,(3,3)),origin='lower',interpolation='nearest',
					vmax=hmax,vmin=-1*hmax,cmap=mpl.cm.__dict__[cmap_name],
					extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
				ax.set_title(r'$\mathrm{\langle z \rangle}$'+'(max %.3f)'%hmax)
			#---PLOT spectrum
			if 'spectrum' in viewnames:
				ax = axes[viewnames.index('spectrum')]
				ax.scatter(datas[tag][sn]['qs'],datas[tag][sn]['ratios'],s=4,c='k',alpha=0.25)
				#---! high cutoff is hard-coded here but needs to be removed to the yaml. we need to get default
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
			if 'spectrum_zoom' in viewnames:
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
			picturesave('fig.%s.%s'%(out_fn,sn),work.plotdir,backup=False,version=True,meta=meta)

@autoplot(plot_super)
def individual_reviews():
	"""
	Plot a few versions of the main figure.
	"""
	plotspec = {
		'coupling_review':{
			'viewnames':['average_height','average_height_pbc','neighborhood_static',
				'neighborhood_dynamic','average_field','example_field','example_field_pbc',
				'spectrum','spectrum_zoom']},
		'coupling_review.simple':{
			'viewnames':['average_height','example_field'],'figsize':(6,6)}}
	for out_fn,details in plotspec.items(): 
		individual_reviews_plotter(out_fn=out_fn,**details)

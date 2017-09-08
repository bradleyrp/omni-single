#!/usr/bin/env python

from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from render.wavevids import plothull
str_types = [str,unicode] if sys.version_info<(3,0) else [str]

if 'data' not in globals():
	avail = ['curvature_field_review','individual_reviews']
	#---alternate method for plotting controls
	plotspecs = work.plots[plotname].get('specs',{})
	routine = plotspecs.get('routine',avail)
	calcname = plotspecs.get('calcname',plotname)
	protein_abstractor_name = plotspecs.get('protein_abstractor_name','protein_abstractor')
	undulations_name = plotspecs.get('undulations_name','undulations')
	#---the curvature undulation coupling data are notably absent from the upstream calculations
	#---...because we pull all upstream sweeps here for comparison
	data,calc = plotload(plotname)
	#---get all upstream curvature sweeps
	ups = work.calc_meta.unroll_loops(work.calcs[calcname],return_stubs=True)[1]
	for up in ups: up['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for unum,up in enumerate(ups):
		#---temporarily set the plots
		#---if the specs value is a string it is PROBABLY in a loop?
		#---! THIS CODE WAS ABANDONED. FOR THE DEXTRAN PROJECT JUST HIDE UPSTREAM PARAMETER SWEEPS
		#---! make this systematic later
		if type(up['specs']['design']) in str_types:
			new_plot_specs = work.calcs[calcname]['specs']['design']['loop'][up['specs']['design']]
		#---no loop means we just pass along the upstream specs to the plots temporarily for plotload
		else: new_plot_specs = up['specs']
		work.plots[plotname]['calculation'] = {plotname:new_plot_specs}
		#---look up the right upstream data
		dat,cal = plotload(calcname)
		tag = up['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
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

if 'curvature_field_review' in routine:
	figsize = (10,10)
	cmap_name = 'RdBu_r'	
	layout = {'out':{'grid':[len(datas),1]},'ins':[{'grid':[1,len(work.sns())]} for d in datas]}
	axes,fig = panelplot(layout,figsize=figsize)
	cmax = max([np.abs(d[sn]['cf']).max() for d in datas.values() for sn in work.sns()])
	for dnum,(tag,dat) in enumerate(datas.items()):
		for snum,sn in enumerate(work.sns()):
			ax = axes[dnum][snum]
			ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name])
			error = datas[tag][sn]['bundle'][sn]['fun']
			ax.set_title('%s\n%s: %.4f'%(work.meta[sn].get('label',sn),tag,error))
	picturesave('fig.curvature_optimized_meta',work.plotdir,backup=False,version=True,
		meta={'maximum C_0':cmax,'sns':work.sns(),'designs':datas.keys()})

if 'individual_reviews' in routine:

	"""
	Reviewing the various features of the curvature coupling fits.
	"""

	def plot_hull_and_trial_centers(ax,n_instances=None,debug_frame=0,color=None):
		"""Plot the protein hull along with the positions of the trial functions."""
		global postdat
		#---! currently only set for a single protein
		trials = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2)
		nframes = len(trials)
		if n_instances!=None: samples = np.arange(0,nframes,nframes/n_instances)
		else: samples = np.array([debug_frame])
		pts = np.concatenate(trials[samples])
		for ptsnum,pts in enumerate(trials[samples]):
			color_this = mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)) if not color else color
			ax.scatter(*pts.T,s=1,c=color_this)
		griddims = datas[tag][sn]['cf'].shape
		#if data[protein_abstractor_name][sn]['data']['vecs']=='readymade_meso_v1':
		#	vecs = data[undulations_name][sn]['data']['vecs'].mean(axis=0)
		#else: 
		#vecs = data[protein_abstractor_name][sn]['data']['vecs'].mean(axis=0)
		#points_protein = data[protein_abstractor_name][sn]['data']['points']
		#---alternate handling for mesoscale
		#if points_protein=='readymade_meso_v1':
		#points_protein = data[protein_abstractor_name][sn]['data']['points_all']
		vecs,points_protein = [postdat[sn][i] for i in ['vecs','points_protein']]
		for ptsnum,pts in enumerate(points_protein[samples][...,:2]):
			color_this = mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)) if not color else color
			plothull(ax,pts,griddims=griddims,vecs=vecs,lw=0,c=color_this)
		ax.set_xlim((0,vecs[0]))
		ax.set_ylim((0,vecs[1]))
		ax.set_aspect(1.0)

	for tag in datas:
		for sn in work.sns():
			status('reviewing curvature for %s'%sn,tag='plot')
			figsize = (10,10)
			cmap_name = 'RdBu_r'	
			#---this plot is an assortment of views
			viewnames = ['average_height','average_height_pbc','neighborhood_static',
				'neighborhood_dynamic','average_field','example_field','example_field_pbc',
				'spectrum','spectrum_zoom']
			axes,fig = square_tiles(len(viewnames),figsize)
			#---several plots use the same data
			#---custom handling for the mesoscale data which lacks vectors from the nanogel
			#if data[protein_abstractor_name][sn]['data']['vecs']=='readymade_meso_v1':
			#	vecs = data[undulations_name][sn]['data']['vecs'].mean(axis=0)
			#else: 
			#vecs = data[protein_abstractor_name][sn]['data']['vecs'].mean(axis=0)
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
			plot_hull_and_trial_centers(ax,debug_frame=debug_frame,color='k')
			ax.set_title('neighborhood, static, frame %d'%debug_frame)
			#---PLOT a composite of several frames of the neighborhood to see the dynamics
			ax = axes[viewnames.index('neighborhood_dynamic')]
			plot_hull_and_trial_centers(ax,n_instances=10)
			ax.set_title('neighborhood, dynamic')
			#---PLOT the average height
			mesh = data[undulations_name][sn]['data']['mesh']
			surf = mesh.mean(axis=0).mean(axis=0)
			surf -= surf.mean()
			hmax = np.abs(surf).max()
			ax = axes[viewnames.index('average_height')]
			ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
				extent=[0,vecs[0],0,vecs[1]])
			#if data[protein_abstractor_name][sn]['data']['points']=='readymade_meso_v1':
			#	mean_prot_pts = data[protein_abstractor_name][sn]['data']['points_all'].mean(axis=0)[:,:2]
			#else: 
			#mean_prot_pts = data[protein_abstractor_name][sn]['data']['points'].mean(axis=0)[:,:2]
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
			####### hicut = work.plots[plotname].get('fitting',{}).get('high_cutoff',1.0)
			hicut = work.calcs[calcname]['specs'].get('fitting',{}).get('high_cutoff',1.0)
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
			#---! added specs???
			meta = calcs[tag][sn]['calcs']['specs']['specs']['design']
			#---add high cutoff (from fitting parameters if defined) to the meta
			meta['high_cutoff'] = hicut
			picturesave('fig.coupling_review.%s'%sn,work.plotdir,backup=False,version=True,meta=meta)

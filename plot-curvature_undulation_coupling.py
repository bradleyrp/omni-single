#!/usr/bin/env python

from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from render.wavevids import plothull

avail = ['curvature_field_review','individual_reviews']
routine = work.plots[plotname].get('specs',{}).get('routine',avail)

if 'data' not in globals():
	#---the curvature undulation coupling data are notably absent from the upstream calculations
	#---...because we pull all upstream sweeps here for comparison
	data,calc = plotload(plotname)
	#---get all upstream curvature sweeps
	ups = work.calc_meta.unroll_loops(work.calcs[plotname],return_stubs=True)[1]
	for up in ups: up['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for unum,up in enumerate(ups):
		#---temporarily set the plots
		work.plots[plotname]['calculation'] = {plotname:up['specs']}
		#---look up the right upstream data
		dat,cal = plotload(plotname)
		tag = up['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])

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
		vecs = data['protein_abstractor'][sn]['data']['vecs'].mean(axis=0)
		points_protein = data['protein_abstractor'][sn]['data']['points']
		for ptsnum,pts in enumerate(points_protein[samples][...,:2]):
			color_this = mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)) if not color else color
			plothull(ax,[pts],griddims=griddims,vecs=vecs,lw=0,c=color_this)
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
				'neighborhood_dynamic','average_field','example_field','example_field_pbc','spectrum']
			axes,fig = square_tiles(len(viewnames),figsize)
			#---several plots use the same data
			vecs = data['protein_abstractor'][sn]['data']['vecs'].mean(axis=0)
			griddims = datas[tag][sn]['cf'].shape
			#---PLOT the mean curvature field
			ax = axes[viewnames.index('average_field')]
			cmax = np.abs(datas[tag][sn]['cf']).max()
			ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
			mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
			ax.scatter(*mean_trial.T,s=1,c='k')
			ax.set_title(r'$\mathrm{\langle C_0(x,y) \rangle}$')
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
			mesh = data['undulations'][sn]['data']['mesh']
			surf = mesh.mean(axis=0).mean(axis=0)
			surf -= surf.mean()
			hmax = np.abs(surf).max()
			ax = axes[viewnames.index('average_height')]
			ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
				extent=[0,vecs[0],0,vecs[1]])
			mean_prot_pts = data['protein_abstractor'][sn]['data']['points'].mean(axis=0)[:,:2]
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
			ax.set_title('example field',fontsize=10)
			#---PLOT periodic view of the example field
			ax = axes[viewnames.index('example_field_pbc')]
			ax.imshow(np.tile(cf_first.T,(3,3)),origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],
				extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
			ax.set_title('example field (max %.3f)'%cmax)
			#---PLOT periodic view of the average height
			ax = axes[viewnames.index('average_height_pbc')]
			ax.imshow(np.tile(surf.T,(3,3)),origin='lower',interpolation='nearest',
				vmax=hmax,vmin=-1*hmax,cmap=mpl.cm.__dict__[cmap_name],
				extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
			ax.set_title(r'$\mathrm{\langle z \rangle}$'+'(max %.3f)'%hmax)
			#---PLOT spectrum
			ax = axes[viewnames.index('spectrum')]
			ax.scatter(datas[tag][sn]['qs'],datas[tag][sn]['ratios'],s=4,c='k',alpha=0.5)
			ax.axhline(1.0,c='k')
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.set_title(r'energy spectrum (relative)')
			ax.grid(True)
			#---no tick marks on anything
			for ax in axes:
				ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
				ax.tick_params(axis='y',which='both',left='off',right='off',labelbottom='off')
			#---the metadata for this plot comes from the design section
			meta = calcs[tag][sn]['calcs']['specs']['design']
			picturesave('fig.coupling_review.%s'%sn,work.plotdir,backup=False,version=True,meta=meta)

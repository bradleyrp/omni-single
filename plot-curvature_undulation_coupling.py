#!/usr/bin/env python

from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from render.wavevids import plothull

avail = ['curvature_field_review','individual_reviews']
routine = work.plots[plotname].get('specs',{}).get('routine',avail)

if 'data' not in globals():
	data,calc = plotload(plotname)
	#except:	data,calc = None,None
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

	def plot_hull_and_trial_centers(ax,n_instances=None,debug_frame=0):
		"""Plot the protein hull along with the positions of the trial functions."""
		#---! currently only set for a single protein
		trials = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2)
		nframes = len(trials)
		if n_instances!=None: samples = np.arange(0,nframes,nframes/n_instances)
		else: samples = np.array([debug_frame])
		pts = np.concatenate(trials[samples])
		for ptsnum,pts in enumerate(trials[samples]):
			ax.scatter(*pts.T,s=1,c=mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)))
		griddims = datas[tag][sn]['cf'].shape
		vecs = data['protein_abstractor'][sn]['data']['vecs'].mean(axis=0)
		points_protein = data['protein_abstractor'][sn]['data']['points']
		for ptsnum,pts in enumerate(points_protein[samples][...,:2]):
			plothull(ax,[pts],
				griddims=griddims,vecs=vecs,lw=0,
				c=mpl.cm.__dict__['jet'](float(ptsnum)/len(samples)))
		ax.set_xlim((0,vecs[0]))
		ax.set_ylim((0,vecs[1]))
		ax.set_aspect(1.0)

	tag = datas.keys()[0]
	tag = 'v9'
	for sn in work.sns():
		status('plotting for %s'%sn)
		figsize = (10,10)
		cmap_name = 'RdBu_r'	
		ntiles = 4
		axes,fig = square_tiles(ntiles,figsize)
		#---global values for consistency between reviews
		#---plot the mean curvature field
		ax = axes[3]
		cmax = np.abs(datas[tag][sn]['cf']).max()
		ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
			vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name])
		mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
		ax.scatter(*mean_trial.T,s=1,c='k')
		#---plot a single instance of the neighborhood (good for debugging)
		ax = axes[1]
		plot_hull_and_trial_centers(ax,debug_frame=0)
		#---plot a composite of several frames of the neighborhood to see the dynamics
		ax = axes[2]
		plot_hull_and_trial_centers(ax,n_instances=10)
		#---plot the average height
		mesh = data['undulations'][sn]['data']['mesh']
		surf = mesh.mean(axis=0).mean(axis=0)
		surf -= surf.mean()
		hmax = np.abs(surf).max()
		ax = axes[0]
		ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'])
		mean_prot_pts = data['protein_abstractor'][sn]['data']['points'].mean(axis=0)[:,:2]
		plothull(ax,[mean_prot_pts],griddims=datas[tag][sn]['cf'].shape,
			vecs=data['protein_abstractor'][sn]['data']['vecs'].mean(axis=0))
		#---no tick marks on anything
		for ax in axes:
			ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
			ax.tick_params(axis='y',which='both',left='off',right='off',labelbottom='off')
		picturesave('fig.TESTING.%s'%sn,work.plotdir,backup=False,version=True,meta={})

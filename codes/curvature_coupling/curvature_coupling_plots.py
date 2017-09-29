#!/usr/bin/env python

"""
Curvature undulation coupling plots.
"""

import re
from codes.curvature_coupling.curvature_coupling import InvestigateCurvature
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from render.wavevids import plothull
from codes.looptools import basic_compute_loop

#---plotting environment
import builtins
for key in builtins._plotrun_specials: 
	globals()[key] = builtins.__dict__[key]

def center_by_protein(sn,surfs,static=False):
	"""
	Align any field to the position of a single protein.
	Assumes that  the protein points are available.
	The static option moves a single field to the center by the average protein position.
	"""
	dat_protein = data[protein_abstractor_name][sn]['data']
	protein_noh = np.array([ii for ii,i in enumerate(dat_protein['names']) if not re.match('^H',i)])
	protein_pts = dat_protein['points_all']
	if protein_pts.shape[1]!=1: raise Exception('cannot center multiple proteins!')
	protein_pts = protein_pts.mean(axis=1)
	protein_cogs = protein_pts[:,protein_noh].mean(axis=1)
	ivecs = dat_protein['vecs']
	protein_cogs += (protein_cogs<0.0)*ivecs - (protein_cogs>=ivecs)*ivecs
	nframes,nx,ny = surfs.shape
	shifts = (protein_cogs[:,:2]/(ivecs.mean(axis=0)[:2]/np.array([nx,ny]))).round().astype(int)
	surfs_shifted = []
	reindexer = np.concatenate(np.transpose(np.meshgrid(*[np.arange(i) for i in [nx,ny]])))
	if not static:
		for fr,surf in enumerate(surfs):
			shift_this = reindexer+(shifts[fr]-(np.array([nx,ny])/2.0).astype(int))
			surfs_shifted.append(surf[shift_this.T[0]%nx,shift_this.T[1]%ny].reshape((nx,ny)))
		return np.array(surfs_shifted)
	else:
		if len(surfs)!=1: raise Exception('the length of surfaces must be one for static')
		else: surf = surfs[0]
		#---single frame version based on the mean protein cog
		shift = (protein_cogs[:,:2].mean(axis=0)/
			(ivecs.mean(axis=0)[:2]/np.array([nx,ny]))).round().astype(int)
		shift_this = reindexer+(shift-(np.array([nx,ny])/2.0).astype(int))
		return surf[shift_this.T[0]%nx,shift_this.T[1]%ny].reshape((nx,ny))

def add_colorbar(ax,im,title,title_top=True):
	#---colorbar
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	if not title_top: axins.set_ylabel(title,rotation=270,labelpad=20)
	else: cbar.set_label(title,labelpad=-20,y=1.1,rotation=0)
	axins.tick_params(axis='y',left='off',right='off',labelright='on')

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

def plothull_protein_center_average(ax,sn,vecs,tag):
	"""Plot the protein hull in the center of the box."""
	dat_protein = data[protein_abstractor_name][sn]['data']
	protein_noh = np.array([ii for ii,i in enumerate(dat_protein['names']) 
		if not re.match('^H',i)])
	protein_pts = dat_protein['points_all']
	if protein_pts.shape[1]!=1: raise Exception('cannot center multiple proteins!')
	mean_prot_pts = protein_pts[:,0].mean(axis=0)[:,:2]
	prot_center = mean_prot_pts - mean_prot_pts.mean(axis=0) + vecs[:2]/2.
	plothull(ax,[prot_center],griddims=datas[tag][sn]['cf'].shape,
		vecs=vecs,c='k',lw=0)

def individual_reviews_plotter(viewnames,out_fn,seep=None,figsize=(10,10),
	horizontal=False,wspace=0.4,**kwargs):
	"""
	Loop over all upstream curvature-undulation coupling calculations and plot a panel of review plots.
	"""
	#---seeping namespace
	if seep!=None:
		for key,val in seep.items(): globals()[key] = val
	for tag in datas:
		for sn in work.sns():
			status('reviewing curvature for %s'%sn,tag='plot')
			cmap_name = 'RdBu_r'
			square_tiles_args = {}
			if wspace!=None: square_tiles_args.update(wspace=wspace)
			axes,fig = square_tiles(len(viewnames),figsize,hspace=0.4,
				favor_rows=horizontal,**square_tiles_args)
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
				im = ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
				ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title('average field')
				add_colorbar(ax,im,title=r'$\mathrm{\langle C_0(x,y) \rangle}$')
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
				im = ax.imshow(surf.T,origin='lower',
					interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
					extent=[0,vecs[0],0,vecs[1]],vmax=hmax,vmin=-1*hmax)
				try: mean_prot_pts = postdat[sn]['points_protein_mean'][:,:2]
				except: mean_prot_pts = data[protein_abstractor_name][sn][
					'data']['points'].mean(axis=0)[:,:2]
				plothull(ax,[mean_prot_pts],griddims=datas[tag][sn]['cf'].shape,vecs=vecs,c='k',lw=0)
				ax.set_title('height profile')
				ax.set_xlabel('x (nm)')
				ax.set_ylabel('y (nm)')
				add_colorbar(ax,im,title=r'$\mathrm{\langle z \rangle \, (nm)}$')
			#---PLOT an example curvature field ('first' is hard-coded, instead of saving each frame)
			if 'example_field' in viewnames:
				ax = axes[viewnames.index('example_field')]
				example_frame = 0
				cf_first = datas[tag][sn]['cf_first']
				im = ax.imshow(cf_first.T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				add_colorbar(ax,im,title=r'$\mathrm{C_0\,({nm}^{-1})}$')
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
				im = ax.imshow(np.tile(cf_first.T,(3,3)),origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],
					extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
				ax.set_title('curvature field (max %.3f)'%cmax_instant)
				add_colorbar(ax,im,title=r'$\mathrm{C_0\,({nm}^{-1})}$')
			#---PLOT periodic view of the average height
			if 'average_height_pbc' in viewnames:
				ax = axes[viewnames.index('average_height_pbc')]
				im = ax.imshow(np.tile(surf.T,(3,3)),origin='lower',interpolation='nearest',
					vmax=hmax,vmin=-1*hmax,cmap=mpl.cm.__dict__[cmap_name],
					extent=[-vecs[0],2*vecs[0],-vecs[1],2*vecs[1]])
				ax.set_title('average height')
				add_colorbar(ax,im,title=r'$\mathrm{\langle z \rangle}$')
				#---track the protein
				protein_pts = data[protein_abstractor_name][sn]['data']['points_all']
				if protein_pts.shape[1]!=1: raise Exception('only one protein')
				prot_traj = protein_pts[:,0].mean(axis=1)[:,:2]
				for shift in np.concatenate(np.transpose(np.meshgrid([-1,0,1],[-1,0,1])))*vecs[:2]:
					ax.scatter(prot_traj[:,0]+shift[0],prot_traj[:,1]+shift[1],c='k',s=1)
				ax.set_xlim((-vecs[0],2*vecs[0]))
				ax.set_ylim((-vecs[1],2*vecs[1]))
			#---PLOT spectrum
			if 'spectrum' in viewnames:
				ax = axes[viewnames.index('spectrum')]
				ax.scatter(datas[tag][sn]['qs'],datas[tag][sn]['ratios'],s=4,c='k',alpha=0.25)
				#---! high cutoff is hard-coded here but needs to be removed to the yaml
				#---! ...we need to get the default
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
			#---PLOT average height profile RELATIVE TO THE PROTEIN POSITION
			if 'average_height_center' in viewnames:
				mesh = data[undulations_name][sn]['data']['mesh']
				surfs = mesh.mean(axis=0)#.mean(axis=0)
				surfs_shifted = center_by_protein(sn,surfs)
				surf = np.mean(surfs_shifted,axis=0)
				surf -= surf.mean()
				hmax = np.abs(surf).max()
				ax = axes[viewnames.index('average_height_center')]
				im = ax.imshow(surf.T,origin='lower',interpolation='nearest',cmap=mpl.cm.__dict__['RdBu_r'],
					extent=[0,vecs[0],0,vecs[1]],vmax=hmax,vmin=-1*hmax)
				try: mean_prot_pts = postdat[sn]['points_protein_mean'][:,:2]
				except: mean_prot_pts = data[protein_abstractor_name][sn]['data']['points'].mean(axis=0)[:,:2]
				#---! still need to move the protein
				if False: plothull(ax,[mean_prot_pts],griddims=datas[tag][sn]['cf'].shape,
					vecs=vecs,c='k',lw=0)
				plothull_protein_center_average(ax,sn,vecs,tag)
				ax.set_title('height profile (centered)')
				ax.set_xlabel('x (nm)')
				ax.set_ylabel('y (nm)')
				add_colorbar(ax,im,title=r'$\mathrm{\langle z \rangle \, (nm)}$')
			#---PLOT the average field shifted as if the protein were in the center
			#---we call this "naive" because it takes a single, average curvature field and moves it 
			#---...to the average position of the protein. if the protein crosses the periodic bounds, 
			#---...then moving *anything* to the average position of the protein makes no sense
			#---...so see the explicit version below called "curvature_field_center"
			if 'average_field_center' in viewnames:
				ax = axes[viewnames.index('average_field_center')]
				nframes = data[undulations_name][sn]['data']['mesh'].shape[1]
				cf_centered = center_by_protein(sn,np.array([datas[tag][sn]['cf']]),static=True)
				im = ax.imshow(cf_centered.T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
				if False: ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title('curvature field (centered,naive)')
				add_colorbar(ax,im,title=r'$\mathrm{C_0\,({nm}^{-1})}$')
			#---PLOT the average curvature field, shifted to the instantaneous position of the protein
			#---...and then averaged. this makes slightly more sense than "average_field_center" but we are 
			#---...still assuming an average and then just smearing it around. see "curvature_field_center"
			if 'average_field_center_smear' in viewnames:
				ax = axes[viewnames.index('average_field_center_smear')]
				nframes = data[undulations_name][sn]['data']['mesh'].shape[1]
				cf_centered = center_by_protein(sn, np.tile(datas[tag][sn]['cf'],(nframes,1,1))).mean(axis=0)
				im = ax.imshow(cf_centered.T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
				if False: ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title('curvature field (centered,smear)')
				add_colorbar(ax,im,title=r'$\mathrm{C_0\,({nm}^{-1})}$')
			#---PLOT the average of the protein-centered instantaneous curvature fields
			#---note that you have to add "store_instantaneous_fields: True" to the design of the calculation
			if 'curvature_field_center' in viewnames:
				ax = axes[viewnames.index('curvature_field_center')]
				nframes = data[undulations_name][sn]['data']['mesh'].shape[1]
				cf_centered = center_by_protein(sn,datas[tag][sn]['cfs']).mean(axis=0)
				im = ax.imshow(cf_centered.T,origin='lower',interpolation='nearest',
					vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name],extent=[0,vecs[0],0,vecs[1]])
				mean_trial = datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2).mean(axis=0)
				if False: ax.scatter(*mean_trial.T,s=1,c='k')
				ax.set_title('curvature field (centered)')
				ax.set_xlabel('x (nm)')
				ax.set_ylabel('y (nm)')
				add_colorbar(ax,im,title=r'$\mathrm{\langle C_0 \rangle \,({nm}^{-1})$')
				plothull_protein_center_average(ax,sn,vecs,tag)
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

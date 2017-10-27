#!/usr/bin/env python

"""
Plot routines that use the lipid_mesh calculation.
"""

@autoload(plotrun)
def load():
	"""Load the data."""
	#---condition for loading the data
	if 'data' not in globals(): 
		#---load once
		global data,calc,data_prot
		data,calc = plotload(plotname,work)
		data_prot,_ = plotload('protein_abstractor',work)

@autoplot(plotrun)
def plot_curvature_maps():
	"""Loop over simulations and make curvature plots."""
	spacing = work.plots['lipid_mesh'].get('specs',{}).get('curvature_map_spacing',2.0)
	for sn in work.sns(): plot_curvature(sn,spacing=spacing)

def plot_curvature(sn,**kwargs):
	"""
	Plot a map of the curvature.
	"""
	global curvs_map_all
	dat = data[sn]['data']
	nframes = int(dat['nframes'])
	def get(mn,fr,name): return dat['%d.%d.%s'%(mn,fr,name)]
	spacing = kwargs.pop('spacing',1.0)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	mvecs = np.mean([get(0,fr,'vec') for fr in range(nframes)],axis=0)
	ngrid = np.round(mvecs/spacing).astype(int)
	#---curvature of each leaflet
	curvs_map_all = [np.zeros((ngrid[0],ngrid[1])) for mn in range(2)]
	for mn in range(2):
		curvs = np.zeros((ngrid[0],ngrid[1]))
		curvs_counts = np.zeros((ngrid[0],ngrid[1]))
		nmol = int(dat['%d.1.nmol'%mn])
		for fr in range(nframes):
			simps = get(mn,fr,'simplices')
			pts = get(mn,fr,'points')
			vec = get(mn,fr,'vec')
			curv = get(mn,fr,'mean')
			#---repack the points in the box
			pts_repack = (pts-np.floor(pts/vec)*vec)
			pts_rounded = (pts_repack/(vec/ngrid)).astype(int)[:,:2]
			curvs[pts_rounded[:nmol,0],pts_rounded[:nmol,1]] += curv[:nmol]
			curvs_counts[pts_rounded[:nmol,0],pts_rounded[:nmol,1]] += 1
		obs = np.where(curvs_counts>0)
		means = curvs[obs]/curvs_counts[obs]
		curvs_map_all[mn][obs] = means

	from mpl_toolkits.axes_grid1.inset_locator import inset_axes
	import scipy
	import scipy.ndimage

	#---average the leaflets
	curvs_map = np.mean(curvs_map_all,axis=0)
	#---! removed the following line which sets nan (white on imshow). not sure what it is for
	#---ironically back to nan for the imshow
	if False: curvs_map[np.where(curvs_counts==0)] = np.nan
	#---! it might be nice to smooth the plots but beware!
	if False: curvs_map = scipy.ndimage.gaussian_filter(curvs_map,sigma=(4,4),order=0)
	vmax = max([abs(j) for j in [curvs_map.max(),curvs_map.min()]])
	ax = plt.gca()
	im = ax.imshow(curvs_map,cmap=mpl.cm.RdBu_r,vmax=vmax,vmin=-1*vmax,origin='lower',
		extent=[0,mvecs[0],0,mvecs[1]])
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	cbar.set_label('$\mathrm{C_0\,({nm}^{-1})$',labelpad=-40,y=1.05,rotation=0)
	ax.set_xlabel('x (nm)')
	ax.set_ylabel('y (nm)')
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	axins.tick_params(axis='y',which='both',left='off',right='off',labelright='on')
	ax.set_title('curvature map')
	if sn in data_prot:
		from render.wavevids import plothull
		points_all = data_prot[sn]['data']['points_all'] 
		points_all_mean_time = points_all.mean(axis=0)
		plothull(ax,points_all_mean_time[...,:2],griddims=curvs_map.shape,vecs=mvecs,c='k',lw=0)	
	picturesave('fig.%s.%s'%('curvature',sn),work.plotdir,backup=False,version=True,meta={})
	plt.close()

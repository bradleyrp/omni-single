#!/usr/bin/env python

"""
Plot routines that use the lipid_mesh calculation.
"""

#! major problem that it loads more upstream data than necessary

import time
import scipy
import scipy.interpolate
import scipy.ndimage
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from base.compute_loop import basic_compute_loop

@autoload(plotrun)
def load():
	"""Load the data."""
	#---condition for loading the data
	if 'data' not in globals(): 
		#---load once
		global data,calc,data_prot
		data,calc = plotload(plotname,work)
		data_prot,_ = plotload('protein_abstractor',work)
		#---compute once
		global survey
		survey = compute_curvature_distributions()

def interpolate_curvatures(fr):
	"""Interpolate curvatures in parallel for the distributions."""
	global raw,fine
	interp = scipy.interpolate.griddata(raw[fr][:,:2],raw[fr][:,2],
		(fine[fr].T[0],fine[fr].T[1]),method='cubic')
	#---! eliminate nan from the edges
	interp[np.isnan(interp)] = 0.0
	return interp

def compute_curvature_distributions():
	"""Compute curvature distributions."""
	spacing = work.plots['lipid_mesh'].get('specs',{}).get('curvature_map_spacing',2.0)
	global data
	survey = {}
	for snum,sn in enumerate(work.sns()):
		dat = data[sn]['data']
		nframes = int(dat['nframes'])
		def get(mn,fr,name): return dat['%d.%d.%s'%(mn,fr,name)]
		ngrid = np.round(data[sn]['data']['vecs'].mean(axis=0)[:2]/spacing).astype(int)
		nmols = [int(dat['%d.1.nmol'%mn]) for mn in range(2)]
		start = time.time()
		survey[sn] = np.zeros(ngrid)
		for mn in range(2):
			#---formulate X,Y,curvature points
			xy = np.array([get(mn,fr,'points')[:nmols[mn],:2] for fr in range(nframes)])
			curvatures = np.array([get(mn,fr,'mean') for fr in range(nframes)])
			global raw,fine
			raw = np.concatenate((np.transpose(xy,(2,0,1)),
				np.reshape(curvatures,(1,nframes,-1)))).transpose(1,2,0)
			#---location of grid points on the unit square
			prop_pts = np.transpose(np.meshgrid(
				range(0,ngrid[0]),range(0,ngrid[1])))/np.array(ngrid).astype(float)
			#---using the new trick for getting points over frames via vectors and points from the unit square
			fine = (np.tile(np.reshape(prop_pts,
				(ngrid[0],ngrid[1],1,2)),(nframes,1))*dat['vecs'][:,:2]).transpose((2,0,1,3))
			interp_accumulate = basic_compute_loop(interpolate_curvatures,
				[dict(fr=fr) for fr in range(nframes)])
			survey[sn] += np.sum(interp_accumulate,axis=0)
		survey[sn] = survey[sn]/(2.*nframes)
	return survey

@autoplot(plotrun)
def plot_curvature_maps():
	"""Loop over simulations and make curvature plots."""
	spacing = work.plots['lipid_mesh'].get('specs',{}).get('curvature_map_spacing',2.0)
	blur = work.plots['lipid_mesh'].get('specs',{}).get('curvature_map_blur',0.0)
	for sn in work.sns(): plot_curvature_basic(sn,spacing=spacing,blur=blur)

def plot_curvature_basic(sn,**kwargs):
	"""
	Plot a map of the curvature.
	Note: you must use a bin spacing of 2.0 or 5.0 and not 0.5 otherwise the curvature ends up biased 
	very strongly in one direction (this is a methodology error of some kind, related to binnning.)
	A Gaussian blur is necessary to get usable information from the map, however this washes out the 
	curvature magnitudes leaving only a relative pattern.
	"""
	dat = data[sn]['data']
	nframes = int(dat['nframes'])
	def get(mn,fr,name): return dat['%d.%d.%s'%(mn,fr,name)]
	spacing = kwargs.pop('spacing',2.0)
	sigma = kwargs.pop('blur',0.0)
	if kwargs: raise Exception('unprocessed kwargs %s'%kwargs)
	mvecs = np.mean([get(0,fr,'vec') for fr in range(nframes)],axis=0)
	ngrid = np.round(mvecs/spacing).astype(int)
	#---curvature of each leaflet
	curvs_mean_all = []
	curvs_map_all = [np.zeros((ngrid[0],ngrid[1])) for mn in range(2)]
	for mn in range(2)[::-1]:
		curvs = np.zeros((ngrid[0],ngrid[1]))
		curvs_counts = np.zeros((ngrid[0],ngrid[1]))
		nmol = int(dat['%d.1.nmol'%mn])
		curvs_mean = []
		for fr in range(nframes):
			simps = get(mn,fr,'simplices')
			pts = get(mn,fr,'points')
			vec = get(mn,fr,'vec')
			curv = get(mn,fr,'mean')
			#---previously implemented a filter
			if False:
				filt = np.where(np.abs(curv)<0.5)
				pts,curv = pts[filt],curv[filt]
			#---repack the points in the box
			pts_repack = (pts-np.floor(pts/vec)*vec)
			pts_rounded = (pts_repack/(vec/ngrid)).astype(int)[:,:2]
			curvs[pts_rounded[:nmol,0],pts_rounded[:nmol,1]] += curv[:nmol]
			curvs_mean.append(curv[:nmol].mean())
			curvs_counts[pts_rounded[:nmol,0],pts_rounded[:nmol,1]] += 1
		curvs_mean_all.append(curvs_mean)
		obs = np.where(curvs_counts>0)
		means = curvs[obs]/curvs_counts[obs]
		curvs_map_all[mn][obs] = means
	#---average the leaflets
	curvs_map = np.mean(curvs_map_all,axis=0)
	#---! removed the following line which sets nan (white on imshow). not sure what it is for
	#---ironically back to nan for the imshow
	if False: curvs_map[np.where(curvs_counts==0)] = np.nan
	#---! it might be nice to smooth the plots but beware!
	if True: curvs_map = scipy.ndimage.gaussian_filter(curvs_map,sigma=(sigma,sigma),order=0)
	vmax = max([abs(j) for j in [curvs_map.max(),curvs_map.min()]])
	ax = plt.gca()
	im = ax.imshow(curvs_map,cmap=mpl.cm.RdBu_r,vmax=vmax,vmin=-1*vmax,origin='lower',
		extent=[0,mvecs[0],0,mvecs[1]])
	axins = inset_axes(ax,width="5%",height="100%",loc=3,
		bbox_to_anchor=(1.05,0.,1.,1.),bbox_transform=ax.transAxes,borderpad=0)
	cbar = plt.colorbar(im,cax=axins,orientation="vertical")
	cbar.set_label(r'$C_0\,({nm}^{-1})$',labelpad=-40,y=1.05,rotation=0)
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

@autoplot(plotrun)
def plot_curvature_distributions(*kwargs):
	"""Plot curvature distributions for all simulations."""
	spacing = work.plots['lipid_mesh'].get('specs',{}).get('curvature_map_spacing',2.0)
	def renorm_unit(vals): return 0.5+(vals/np.abs(vals).max())/2.
	ax = plt.gca()
	for sn in work.sns():
		#---ignore ridiculous curvatures which indicate a mesh error by requesting specific bins
		bins = np.linspace(-0.5,0.5,100)
		counts,vals = np.histogram(survey[sn].reshape(-1),bins=bins)
		ax.plot((vals[1:]+vals[:-1])/2.,counts,label=work.meta.get(sn,{}).get('label',sn))
	ax.axvline(0.0,lw=1,c='k',zorder=1)
	ax.legend()
	picturesave('fig.%s'%('curvature.distribution'),work.plotdir,backup=False,version=True,meta={})
	plt.close()

if __name__=='__replotting__':

	import scipy
	import scipy.interpolate

	#! incoming
	spacing = 2.0
	sn = work.sns()[0]

	dat = data[sn]['data']
	nframes = int(dat['nframes'])
	def get(mn,fr,name): return dat['%d.%d.%s'%(mn,fr,name)]
	ngrid = np.round(data[sn]['data']['vecs'].mean(axis=0)[:2]/spacing).astype(int)
	nmols = [int(dat['%d.1.nmol'%mn]) for mn in range(2)]
	def renorm_unit(vals): return 0.5+(vals/np.abs(vals).max())/2.

	start = time.time()
	if 'survey_debug' not in globals():
		survey_debug = np.zeros(ngrid)
		for mn in range(2):
			#---formulate X,Y,curvature points
			xy = np.array([get(mn,fr,'points')[:nmols[mn],:2] for fr in range(nframes)])
			curvatures = np.array([get(mn,fr,'mean') for fr in range(nframes)])
			raw = np.concatenate((np.transpose(xy,(2,0,1)),
				np.reshape(curvatures,(1,nframes,-1)))).transpose(1,2,0)
			#---location of grid points on the unit square
			prop_pts = np.transpose(np.meshgrid(
				range(0,ngrid[0]),range(0,ngrid[1])))/np.array(ngrid).astype(float)
			#---using the new trick for getting points over frames via vectors and points from the unit square
			fine = (np.tile(np.reshape(prop_pts,
				(ngrid[0],ngrid[1],1,2)),(nframes,1))*dat['vecs'][:,:2]).transpose((2,0,1,3))
			for fr in range(nframes):
				status('interpolating',i=mn*nframes+fr,looplen=2*nframes,start=start)
				interp = scipy.interpolate.griddata(raw[fr][:,:2],raw[fr][:,2],
					(fine[fr].T[0],fine[fr].T[1]),method='cubic')
				#---! eliminate nan from the edges
				interp[np.isnan(interp)] = 0.0
				survey_debug += interp
		survey_debug = survey_debug/(2.*nframes)
	bins = np.linspace(-0.5,0.5,100)
	#! got good lucking distributions which were very slightly off-center if you plot interp instead, below
	#! need to debug the single interp with a spike slightly off-center
	#! need to debug the way it gets totally washed out when you sum it up
	#! and even then, we would expect that according to the washed out plots, a whole survey of curvature in one distribution will be really insufficient
	counts,vals = np.histogram(survey_debug.reshape(-1),bins=bins)
	plt.plot((vals[1:]+vals[:-1])/2.,counts);
	plt.axvline(0.0)
	picturesave('debug8b',work.plotdir);

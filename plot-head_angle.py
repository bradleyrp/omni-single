#!/usr/bin/env python

from scipy.stats import kde
import scipy.spatial
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import subprocess

#---block: import the post-processed data	
if 'data' not in globals(): 
	sns,(data,calc) = work.sns(),plotload(plotname,work)
	#---! hacking right now to get times
	_,(data_contacts,_) = work.sns(),plotload('contacts',work)

#---block: settings
fsbase = 20
nbins = 50
nlevels = 25
xmin,xmax,ymin,ymax = -180,180,-180,180
color_map_name = 'gnuplot2 pink jet'.split()[-1]
plot_simple = True
angle_ticks = np.arange(-150,200,50)
cmap = mpl.cm.get_cmap(color_map_name)
routine = ['simple','head_angle_vs_distance','target_vs_distance'][:1]

#---block: utility functions
def angle_image_format(ax,sn=None):
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	plt.setp(ax.xaxis.get_majorticklabels(),rotation=90,fontsize=fsbase-2)
	plt.setp(ax.yaxis.get_majorticklabels(),rotation=0,fontsize=fsbase-2)
	ax.axhline(0,lw=1.5,c='w',zorder=4)
	ax.axvline(0,lw=1.5,c='w',zorder=4)
	#---! turning off twinx for now because it messes up the proportions
	if sn and False:
		ax_side = ax.twinx()
		ax_side.set_yticks([])
		ax_side.set_ylabel(work.meta[sn]['composition_name'],fontsize=fsbase,rotation=270,labelpad=25)

#---block: add the kernel density estimates
if 'postdat' not in globals():

	#---distances in Angstroms where -1 means no cutoff (sort this for the gif later)
	distances_sweep = list(np.arange(10.0,21.0,1.0))+[-1]
	lenscale = 10.0

	#---load extra data for distances calculation
	#---! please make a note of how this works in the meta file
	#---! this is obviously contintent on having a protein
	data_lipids,_ = plotload('lipid_abstractor',work)
	data_proteins,_ = plotload('protein_abstractor',work)

	postdat,postdat_dists = {},{}
	for pnum,sn in enumerate(sns):

		#---compute nearest distances
		nframes,nlipids = data_lipids[sn]['data']['points'].shape[:2]
		inds_pip2 = np.where(data_lipids[sn]['data']['resnames']==work.meta[sn]['ptdins_resname'])[0]
		distances = np.zeros((nframes,len(inds_pip2)))
		valid_frames = np.ones(nframes).astype(bool)
		#---get nearest distances for each frame
		for fr in range(nframes):
			status('finding closest distance between lipids and proteins',i=fr,looplen=nframes,tag='compute')
			#---get protein and lipid points for this frame
			pts_fore = lipid_pts = data_lipids[sn]['data']['points'][fr][inds_pip2]
			#---works for one protein only right now
			nprot = 0
			pts_back = prot_pts = data_proteins[sn]['data']['points_all'][fr][nprot]
			vec = data_lipids[sn]['data']['vecs'][fr]
			try: tree = scipy.spatial.ckdtree.cKDTree(pts_back,boxsize=np.concatenate((vec,vec)))
			except: 
				#---very rarely we have frames with PBC problems
				valid_frames[fr] = False
				continue
			close,nns = tree.query(pts_fore,k=1)
			distances[fr] = close
		postdat_dists[sn] = {'distances':distances,'valid_frames':valid_frames}

		#---loop over distance cutoffs
		postdat[sn] = {}
		for distance_cutoff in distances_sweep:

			status('kernel density estimates %s'%sn,looplen=len(sns),i=pnum,tag='plot')
			theta,phi = [data[sn]['data'][i] for i in ['theta','phi']]
			#---HACK for backwards PHI in symmetric systems
			#---! deprecated: factor = -1.0 if work.meta[sn]['composition_name']=='symmetric' else 1.0
			#---! ...however this needs checked!
			factor = 1.0
			phi = phi*factor
			raw = np.array([theta,phi]).T
			if distance_cutoff!=-1:
				subsel = np.where(postdat_dists[sn]['distances']<=distance_cutoff/lenscale)
				catted = np.transpose((theta[subsel],phi[subsel]))
			else:
				#---no averaging: everything goes into one bin
				catted = np.concatenate(raw)
			if len(catted)==0: raise Exception('no targets within %.1f'%distance_cutoff)
			#---KDE
			x,y = catted.T
			k = kde.gaussian_kde(catted.T)
			xi, yi = np.mgrid[xmin:xmax:nbins*1j,ymin:ymax:nbins*1j]
			zi = k(np.vstack([xi.flatten(), yi.flatten()]))
			postdat[sn][distance_cutoff] = dict(xi=xi,yi=yi,zi=zi,x=x,y=y)

#---block: make a thumbnail plot
if 'simple' in routine:

	plotspec = {'fs_suptitle':20}
	gif_fns = []

	#---sweep over distance cutoffs
	for distance_cutoff in distances_sweep:

		zmax = dict([(sn,postdat[sn][distance_cutoff]['zi'].max()) for sn in sns])
		axes,fig = square_tiles(len(sns),figsize=(10,10),wspace=0.5,hspace=0.5)
		for pnum,sn in enumerate(sns):
			status('plotting %s'%sn,looplen=len(sns),i=pnum,tag='plot')
			#---! previously used a custom function for getting the right axis
			ax = axes[pnum]
			xi,yi,zi = [postdat[sn][distance_cutoff][key] for key in ['xi','yi','zi']]
			ax.pcolormesh(xi,yi,zi.reshape(xi.shape),vmin=0,vmax=zmax[sn],cmap=cmap)
			cs = ax.contourf(xi,yi,zi.reshape(xi.shape),vmax=zmax[sn],vmin=0,
				levels=np.linspace(0,zmax[sn],nlevels),extend='both',origin='lower',
				lw=2,zorder=3,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
			#---! whitespace problem
			ax.set_aspect('equal')
			#---got latex errors? remove underscores from the name or use dollar signs
			title = work.meta[sn].get('label',sn)+(
				'\n(cholesterol)' if work.meta[sn].get('cholesterol',False) else '')
			ax.set_title(title)
			angle_image_format(ax,sn)
			ax.set_xlim(xmin,xmax)
			ax.set_ylim(ymin,ymax)
			ax.set_xlabel(r'$\mathrm{\theta}$',fontsize=fsbase-2,labelpad=-15)
			ax.set_ylabel(r'$\mathrm{\phi}$',fontsize=fsbase-2,labelpad=-15)
		suptitle = fig.suptitle(r'$\mathrm{{PIP}_{2}}$ inositol tilt $\mathrm{(r\leq%.f\AA)}$'%
			distance_cutoff,fontsize=plotspec['fs_suptitle'])
		fn = 'fig.head_angle_thumbs%s'%('' if distance_cutoff==-1 else '.cutoff_%.1f'%distance_cutoff)
		gif_fns.append(fn)
		picturesave(fn,work.plotdir,backup=False,version=False,meta={},extras=[suptitle])
		plt.close()
	gif_fns_resort = ['%s.png'%i for i in gif_fns][::-1]+['%s.png'%gif_fns[0] for i in range(len(gif_fns))]
	subprocess.call("convert -delay 20 -loop 0 %s fig.head_angle_thumbs_by_distance.gif"%(' '.join(gif_fns_resort)),cwd=work.plotdir,shell=True)

#---block: head_angle vs distance correlation
if 'head_angle_vs_distance' in routine:

	#---works for one protein only right now
	nprot = 0
	xmin,xmax = 0.,5.0
	ymin,ymax = -180.0,180.0
	cutoff = 1.1
	angle_distance_style = ['means','all_each'][-1]

	#---create plots for both angles
	for angle in ['theta','phi']:

		axes,fig = square_tiles(len(sns),figsize=(10,10))
		for pnum,sn in enumerate(sns):
			status('plotting %s'%sn,looplen=len(sns),i=pnum,tag='plot')
			distances,valid_frames = [postdat_dists[sn][i] for i in ['distances','valid_frames']]
			#---correlate angle and distance for each lipid for each frame
			ax = axes[pnum]
			if angle_distance_style=='means':
				reducer = lambda x: x.mean(axis=0)
				ys = reducer(data[sn]['data'][angle][valid_frames])
				xs = reducer(distances[valid_frames])
				ax.scatter(xs,ys)
			elif angle_distance_style=='all_each':
				nlipids = distances[valid_frames].shape[1]
				colors = [mpl.cm.get_cmap('jet')(i) for i in np.linspace(0,1.0,nlipids)]
				np.random.shuffle(colors)
				for nl in range(nlipids):
					ys = data[sn]['data'][angle][valid_frames,nl]
					xs = distances[valid_frames,nl]
					ax.scatter(xs,ys,color=colors[nl])
				#---save the entire set of points for later
				reducer = lambda x: x.reshape(-1)
				xs,ys = [reducer(i) for i in [distances[valid_frames],data[sn]['data'][angle]]]
			else: raise Exception('unclear style')
			ax.set_xlim((xmin,xmax))
			ax.set_ylim((ymin,ymax))
			ax.set_aspect(abs((xmax-xmin)/(ymax-ymin))/1.0)
			title = work.meta[sn].get('label',sn)+(
				'\n(cholesterol)' if work.meta[sn].get('cholesterol',False) else '')
			ax.set_title(title)
			ax.set_xlabel('distance (nm)')
			ax.set_ylabel({'theta':r'$\mathrm{\theta}$','phi':r'$\mathrm{\phi}$'}[angle])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.axvline(cutoff,ymin=0,ymax=1.0,color='k',lw=1)
			ax.axhline(0.0,xmin=0,xmax=1,c='k',lw=1.5)
			axins = inset_axes(ax,width="30%",height="30%",loc=1)
			nears_means = [ys[xs<cutoff].mean(),ys[xs>=cutoff].mean()]
			nears_stds = [ys[xs<cutoff].std()/np.sqrt(np.sum(xs<cutoff)),
				ys[xs>=cutoff].std()/np.sqrt(np.sum(xs>=cutoff))]
			axins.bar([0.5,1.5],nears_means,yerr=nears_stds)
			axins.set_xticks([0.5,1.5])
			axins.set_xticklabels([r'$\mathrm{%s %.1f}$'%(i,cutoff) for i in ['<','\geq']])
			axins.axhline(0.0,xmin=0,xmax=1,c='k',lw=1.5)
		picturesave('fig.angle_by_prot_r_%s'%angle,work.plotdir,backup=False,version=False,meta={})
		plt.close()

#---block: make a target-vs-distance plot
if 'target_vs_distance' in routine:

	plotspec = {'figsize':(8,8),'fs_legend':14}

	axes,fig = square_tiles(len(sns),figsize=plotspec['figsize'],wspace=0.4,hspace=0.4)

	for snum,sn in enumerate(sns):
		ax = axes[snum]
		title = work.meta[sn].get('label',sn)+(
			'\n(cholesterol)' if work.meta[sn].get('cholesterol',False) else '')
		ax.set_title(title)
		for distance_cutoff in [j for j in distances_sweep if j>0]:
			valid_frames = postdat_dists[sn]['valid_frames']
			distances = postdat_dists[sn]['distances'][valid_frames]
			series = (distances<(distance_cutoff/lenscale)).sum(axis=1)
			#---! hacking with contacts for now until we add times back to the other calculations
			this_times = data_contacts[sn]['data']['times'][valid_frames]/1000.0
			ax.plot(this_times,series,color='black',
				alpha=float(distance_cutoff)/max(distances_sweep),
				label=r'$\mathrm{r\leq%.1f\AA}$'%distance_cutoff)
		ax.set_xlabel('time (ns)')
		ax.set_ylabel('nearby $\mathrm{{PIP}_{2}}$')

	legend = axes[-1].legend(loc='upper left',fontsize=plotspec['fs_legend'],
		bbox_to_anchor=(1.05,0.0,1.,1.),shadow=True,fancybox=True,title='centers of mass')	

	picturesave('fig.head_angle_thumbs_distance.cutoff_%.1f'%distance_cutoff,
		work.plotdir,backup=False,version=False,meta={},extras=[legend])

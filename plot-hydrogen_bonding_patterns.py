#!/usr/bin/env python

"""
Maps of hydrogen bonding paterns in lipid bilayers.
"""

from config import bash

routine_all = ['test']
routine = work.plots.get('hydrogen_bonding_patterns',{}).get('routine',routine_all)

if 'test' in routine:
	#---! dev
	if 'data' not in globals():
		data,calc = plotload('hydrogen_bonding_patterns',work)

	#---set smooth_window to zero for no smoothing
	#---! smoothing is broken because PBCs. DO NOT USE SMOOTHING.
	smooth_window = 0

	bond_spec_keys = 'resname_1 resid_1 name_1 resname_2 resid_2 name_2 name_h'
	sns = work.sns()

	#---loop over simulations
	for sn in sns:
		lipid_mesh = work.plotload('ptdins_snapshots_mesh',sns=[sn])[0][sn]['data']
		out_dn = os.path.join(work.plotdir,'mov.hydrogen_bonding_pattern.%s'%sn)
		if os.path.isdir(out_dn): status('refusing to overwrite %s'%out_dn,tag='warning')
		else:
			os.mkdir(out_dn)
			dat_lipids = data['lipid_abstractor'][sn]['data']
			dat_hbonds = data['hydrogen_bonding'][sn]['data']
			#---get the top monolayer
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			#---construct several mappings between relevant indices
			def mapback(seq):
				"""Hash a list of numbers back to their indices."""
				return dict([(v,k) for k,v in zip(np.arange(len(seq)),seq)])
			#---we use 'i' to refer to the master index in lipid abstractor resname, resid, mono lists
			#---the m2i mapping converts top (0) monolayer index into master index
			#---! NOTE THAT SAPI has a different top monolayer index in the lipid abstractor data !!!
			#---! ...this might explain the flipped head-tail angle plots, which paul said were intuitive
			#---! ...but which always flummoxed ryan and possibly david
			top_mono = work.meta[sn].get('index_top_monolayer',0)
			m2i = np.where(data['lipid_abstractor'][sn]['data']['monolayer_indices']==top_mono)[0]
			resids,resnames = [data['lipid_abstractor'][sn]['data'][k] for k in ['resids','resnames']]
			#---the r2i mapping converts resid into master index
			r2i = mapback(resids)
			#---the i2m mapping converts master index into mesh index
			i2m = mapback(m2i)
			#---the m2r mapping converts mesh index into resid and r2m reverses it
			m2r = np.array([resids[j] for j in m2i])
			r2m = mapback(m2r)
			resnames_all = np.array([resnames[j] for j in m2i])
			#---get the hydrogen bonds
			bonds,obs,valid_frames = [dat_hbonds[k] for k in ['bonds','observations','valid_frames']]
			nframes = len(valid_frames)
			#---pre-cache the points for smoothing
			if smooth_window: 
				raise Exception('smoothing has been removed because PBCs')
				smooth_window = 20
				nmol = len(m2i)
				#---note that depending on the PBC links we get a variable number of 
				points_inside = np.array([lipid_mesh['%d.%d.points'%(top_mono,fr)][:nmol] 
					for fr in range(nframes)])
				windows = np.array([np.arange(j,j+smooth_window) 
					for j in np.arange(0,nframes-smooth_window)])
				points_inside_smooth = np.array([points_inside[w].mean(axis=0) for w in windows])
			for frameno,fr in enumerate(valid_frames):
				status('drawing for %s'%sn,i=fr,looplen=nframes,tag='render')
				dot_size_small,dot_size_med,dot_size_large = 8,16,25
				#---copy points so that we can smooth the points in the box
				points = np.array(lipid_mesh['%d.%d.points'%(top_mono,fr)])
				if smooth_window: points[:nmol] = points_inside_smooth[fr]
				ax = plt.subplot(111)
				vec = dat_lipids['vecs'][fr]
				#---get the hydrogen bonds for this plot
				resname_inds = np.array([bond_spec_keys.split().index(j) 
					for j in ['resid_%d'%i for i in [1,2]]])
				#---use frameno here because obs is over valid_frames only
				bonds_this_all = bonds[np.where(obs[frameno])[0]][:,resname_inds].astype(int)
				#---filter out intramolecular hyrogen bonds
				bonds_this = bonds_this_all[np.where(bonds_this_all[:,0]!=bonds_this_all[:,1])[0]]
				#---! note no check for resname type
				#---get bonds that are in the mesh
				bonds_this_valid = np.where([np.all([np.in1d(i,r2m.keys()) 
					for i in b]) for b in bonds_this])[0]
				ghosts = lipid_mesh['%d.%d.ghost_ids'%(top_mono,fr)]
				#---draw a line for each bond
				for lnum,link_resids in enumerate(bonds_this[bonds_this_valid]):
					verts = np.array([r2m[j] for j in link_resids])
					points_this = points[verts,:2]
					if not np.any(points_this.ptp(axis=0)>=vec[:2]/2.):
						ax.plot(*points_this.T,c='k',lw=1,zorder=2)
					else:
						targs = np.where(np.in1d(ghosts,verts))[0]
						for target in targs:
							for urp in link_resids:
								possible_long_bond = np.array([target,r2m[urp]])
								try:
									if not np.any((points[possible_long_bond].ptp(axis=0)>=vec/2.)[:2]):
										ax.plot(*points[possible_long_bond][:,:2].T,c='k',lw=1,zorder=2)		
								except: pass
				#---formulate the color list
				colors = np.array([colorize(work.meta[sn],resname=r) for r in resnames_all[ghosts]])
				#---stray dots are transparent
				ax.scatter(*points[:,:2].T,s=dot_size_small,color=colors,alpha=0.35,zorder=3,lw=0)
				#---active dots have color
				highlit = np.array([r2m[j] for j in np.unique(bonds_this[bonds_this_valid])])
				#---note that only the non-ghost points can be highlit so the color list below is complete
				if highlit.max()>=m2r.shape[0]:
					raise Exception('some highlit points are not in the mesh keys!')
				ax.scatter(*points[highlit,:2].T,lw=0.25,s=dot_size_small,zorder=4,edgecolor='k',
					c=[colorize(work.meta[sn],resname=r) for r in resnames_all[highlit]])
				#---highlight the PIP2
				highlit_pip2 = np.where(np.in1d(resnames_all[highlit],
					work.vars['selectors']['resnames_PIP2']+['PtdIns']))[0]
				ax.scatter(*points[highlit[highlit_pip2],:2].T,lw=0.25,
					s=dot_size_large,zorder=5,edgecolor='k',
					c=[colorize(work.meta[sn],resname=r) for r in resnames_all[highlit[highlit_pip2]]])
				#---plot a box
				lims_prop = 0.1
				ax.set_xlim((vec[0]*-lims_prop,vec[0]*(1.+lims_prop)))
				ax.set_ylim((vec[1]*-lims_prop,vec[1]*(1.+lims_prop)))
				box_prop = (lims_prop*2+1.0)/lims_prop
				ax.axhline(0,xmin=1./box_prop,xmax=(box_prop-1.0)/box_prop,c='k',lw=1)
				ax.axhline(vec[1],xmin=1./box_prop,xmax=(box_prop-1.0)/(box_prop),c='k',lw=1)
				ax.axvline(0,ymin=1./box_prop,ymax=(box_prop-1.0)/(box_prop),c='k',lw=1)
				ax.axvline(vec[0],ymin=1./box_prop,ymax=(box_prop-1.0)/(box_prop),c='k',lw=1)
				ax.set_aspect('equal')
				ax.set_xticks([])
				ax.set_yticks([])
				ax.tick_params(axis=u'both',which=u'both',length=0)
				ax.axis('off')
				picturesave('snap.%05d'%frameno,out_dn,
					backup=False,version=True,meta={},extras=[],loud=False)
			#---render when complete
			try:
				# https://superuser.com/questions/1005315/interpolation-with-ffmpeg
				cmd = 'ffmpeg -i "snap.%05d.v1.png" '+'mov.hydrogen_bonding_pattern.%s'%sn+'.mp4'
				bash(cmd,cwd=out_dn)
			except: status('failed to render the video. try "%s" in %s'%(cmd,out_dn))
		del lipid_mesh

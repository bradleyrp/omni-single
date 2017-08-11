#!/usr/bin/env python

"""
Consolidated snapshot codes.
"""

import shutil
import scipy
import scipy.stats
try: from codes import vmdmake
except: raise Exception(
	'cannot import vmdmake. try `git clone http://github.com/bradleyrp/amx-vmd %s`'%
	os.path.join(os.getcwd(),'calcs','codes','vmdmake'))
#---common tools are housed in the snapshotter module
from codes.snapshotter import *
import matplotlib.patheffects as path_effects
#---! for some reason we are importing post twice on this plot. because two spots?
#---one-off settings for multiple sections: head_angle_snaps and head_angle_snaps_detail
kwargs_vmdmake = {'CGBONDSPATH':'~/libs/cg_bonds.tcl','GMXDUMP':'~/libs/gmxdump'}

routine = [
	#---head_snaps_detail depends on head_angle_snaps
	'head_angle_snaps','head_angle_snaps_detail',
	'common_hydrogen_bonding',
	'common_hydrogen_bonding_specific',
	'consolidated_hydrogen_bonding'][-1:]

if 'data' not in globals():
	data,calc = plotload(plotname)
	sns = work.sns()
	if False:
		sns,(data_tilt,calc) = work.sns(),plotload('head_angle',work)
		#---! not even used!
		data_salt,calc_salt = plotload('salt_bridges',work)
		#---! was thinking of checking proximity to box
		#---! ...with the lipid abstractor: data_points,calc_points = plotload('lipid_abstractor',work)

if 'head_angle_snaps' in routine:

	sns_special = ['membrane-v531','membrane-v532']
	sns_top = [s for s in sns if s not in sns_special]
	#---plot settings
	fsbase = 20
	nbins = 50
	nlevels = 25
	xmin,xmax,ymin,ymax = -180,180,-180,180
	color_map_name = 'gnuplot2 pink jet'.split()[-1]
	plot_simple = False
	angle_ticks = np.arange(-150,200,50)
	cmap = mpl.cm.get_cmap(color_map_name)
	#---snapshotting positions
	positions = ['modal','modal_low']
	#---snapshot naming tag
	tag_head_angle = 'v5'
	#---modal ranking number
	exemplar_rank = 1

	if 'postdat_head_angle' not in globals():
		postdat_head_angle = {}
		for pnum,sn in enumerate(sns):
			status('kernel density estimates %s'%sn,looplen=len(sns),i=pnum,tag='plot')
			theta,phi = [data_tilt[sn]['data'][i] for i in ['theta','phi']]
			#---HACK for backwards PHI in symmetric systems
			factor = -1.0 if work.meta[sn]['composition_name']=='symmetric' else 1.0
			phi = phi*factor
			raw = np.array([theta,phi]).T
			#---no averaging: everything goes into one bin
			catted = np.concatenate(raw)
			#---KDE
			x,y = catted.T
			k = scipy.stats.kde.gaussian_kde(catted.T)
			xi, yi = np.mgrid[xmin:xmax:nbins*1j,ymin:ymax:nbins*1j]
			zi = k(np.vstack([xi.flatten(), yi.flatten()]))
			postdat_head_angle[sn] = dict(xi=xi,yi=yi,zi=zi,x=x,y=y)

	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],'fig.head_angle.%s'%tag_head_angle)
	if not os.path.isdir(tempdir): os.mkdir(tempdir)
	status('snapshots dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')

if 'head_angle_snaps' in routine:

	#---pick select postions for plotting
	reps = {}
	for sn in sns_special:
		xi,yi,zi,x,y = [postdat_head_angle[sn][key] for key in ['xi','yi','zi','x','y']]
		nframes,nlipids = data_tilt[sn]['data']['theta'].shape
		for pos in positions:
			if pos=='modal':
				ind = np.unravel_index(np.argmax(zi.reshape(xi.shape)),xi.shape)
				theta,phi = xi[ind],yi[ind]
			#---somewhat hackish
			elif pos=='modal_low':
				#---actually ultra hackish
				xyz = np.transpose((xi.reshape(-1),yi.reshape(-1),zi));
				theta,phi,_ = next(i for i in sorted(xyz,key=lambda a_entry:a_entry[2])[::-1]  
					if i[0]<0 and i[1]<0)
			frame,resid_rel = exemplar_lipid(theta,phi,x,y,nframes,nlipids,rank=exemplar_rank)
			reps[(sn,pos)] = {'frame':frame,'resid_rel':resid_rel,'theta':theta,'phi':phi}
			
	rendered_fns = {}
	for sn in sns_special[:]:

		#---modified the calculation output to include the slices for use-cases like this one
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(calc['extras'][sn]['slice_path'],suf))
			for suf in ['gro','xtc']]
		#---get the tpr from the raw data
		tpr = work.raw.get_last(sn,subtype='tpr')

		#---precompute snapshot filenames and plot if any are missing
		snapshot_fns = [snapshot_namer(sn=work.raw.prefixer(sn_this),tag=tag,pos=pos) 
			for (sn_this,pos),val in reps.items() if sn_this==sn]
		if any([not os.path.isfile(os.path.join(tempdir,i+'.png')) for i in snapshot_fns]):

			#---each render needs to start from the beginning because we move things around for the good side
			#---loop over desired lipid/frame pairs for rendering
			for key in [r for r in reps if r[0]==sn]: 

				import ipdb;ipdb.set_trace()

				#---run VMD to make the snapshots via vmdmake
				view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
					frames='',res=(2000,2000),**kwargs_vmdmake)
				view.do('load_dynamic','standard','bonder')

				val = reps[key]
				frame,resid_rel,angle = [val[i] for i in ['frame','resid_rel','theta']]
				sn,pos = key
				view['snapshot_filename'] = snapshot_namer(sn=work.raw.prefixer(sn),pos=pos,tag=tag)
				view.command("set lipids_select [atomselect top \"resname %s\"]"%
					work.meta[sn]['ptdins_resname'])
				view.command("set my_index [lindex [lsort -unique [$lipids_select get resid]] %d]"%
					resid_rel)
				headsel = '(%s)'%' or '.join(['name %s'%i for i in work.vars['head_atoms'].split()])
				view.select(**{'me':"resname %s and resid $my_index and not hydrogen and %s"%(
					work.meta[sn]['ptdins_resname'],headsel),
					'smooth':False,'style':'licorice','goodsell':True})
				view.command("animate goto %d"%frame)
				view.do('reset','xview')
				view.command('set reverse_view %d'%(0 if angle>0 else 1))
				#---! EXTREMELY IMPORTANT NOTE: DON'T MESS WITH SELECTIONS WITHOUT REMEMBERING
				#---! ...that the GOODSIDE code rotates the molecules to get the right view
				#---! ...since we cannot rotate the camera
				view.command(re.sub('PIP2RESNAME',work.meta[sn]['ptdins_resname'],vmd_render_good_side))
				view.set_color_cursor({'MG':'pink','Cal':'blue'}[work.meta[sn]['cation']])
				view.select(**{'near_ions':
					"name %s and within 4.6 of (resname %s and resid $my_index)"%(
					work.meta[sn]['cation'],
					work.meta[sn]['ptdins_resname']),
					'smooth':False,'style':'Beads 0.6 12.0','goodsell':True,'color_specific':True})
				#---skipping water because it wasn't drawing bonds after transformations, other problems
				if False: view.select(**{'hydration_shell':
					"water and (same residue as (within 3.01 of (name %s "%work.meta[sn]['cation']+
					"and within 4.6 of (resname %s and resid $my_index))))"%(
					work.meta[sn]['ptdins_resname']),
					'smooth':False,'style':'cpk_water','goodsell':True})
				view.command('scale by 0.9')
				view.do('snapshot')
				view.command("mol delrep 1 top")
				view.command("mol delrep 0 top")
				#---! YOU MUST TERMINATE THE REPS IN THIS LOOP FOR ALL SELECTIONS!
				view.show(quit=True)

if 'head_angle_snaps_detail' in routine:

	if 'head_angle_snaps' not in routine: 
		raise Exception('you must include head_angle_snaps in routine before head_angle_snaps_detail')

	#---rows and columns on the left side
	nrl,ncl = 4,2
	axes,fig = panelplot(figsize=(16,12),
		layout={'out':{'grid':[2,2],'wspace':0.0,'hspace':0.3,'wratios':[1,2]},
		'ins':[
			{'grid':[2,2],'hratios':[4,1][::-1],'wratios':[1,4][::-1],'wspace':0.1,'hspace':0.1},
			{'grid':[1,2],'wspace':0},
			{'grid':[2,2],'hratios':[4,1][::-1],'wratios':[1,4][::-1],'wspace':0.1,'hspace':0.1},
			{'grid':[1,2],'wspace':0},
			]})
	axn = {'contour':1,'y':0,'x':3,'del':2}
	axn = {'contour':2,'y':3,'x':0,'del':1}

	zmax = max([postdat_head_angle[sn]['zi'].max() for sn in sns_special])
	for ss,sn in enumerate(sns_special):
		pnum = [0,2][ss]
		ax = axes[pnum][axn['contour']]
		xi,yi,zi,x,y = [postdat_head_angle[sn][key] for key in ['xi','yi','zi','x','y']]
		ax.pcolormesh(xi,yi,zi.reshape(xi.shape),vmin=0,vmax=zmax,cmap=cmap)
		ax.set_xlim(xmin,xmax)
		ax.set_ylim(ymin,ymax)

		labels = {'theta':r'tilt angle $\mathrm{\theta}$',
			'phi':r'rotation angle $\mathrm{\phi}$'}
		ax.set_xlabel(labels['theta'],fontsize=fsbase)
		ax.set_ylabel(labels['phi'],fontsize=fsbase)

		cs = ax.contourf(xi,yi,zi.reshape(xi.shape),vmax=zmax,vmin=0,levels=np.linspace(0,zmax,nlevels),
			extend='both',origin='lower',lw=2,zorder=3,cmap=cmap,extent=[xmin,xmax,ymin,ymax])
		xbins = np.arange(xmin,xmax,(xmax-xmin)/nbins)
		ybins = np.arange(ymin,ymax,(ymax-ymin)/nbins)

		angle_image_format(ax,fsbase=fsbase)

		ax = axes[pnum][axn['x']]
		ax.hist(x,bins=xbins,color='gray',lw=1,edgecolor='w')

		ax.spines['left'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])

		ax = axes[pnum][axn['y']]
		ax.hist(y,bins=ybins,orientation='horizontal',color='gray',lw=1,edgecolor='w')

		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.yaxis.set_ticks_position('none')
		ax.xaxis.set_ticks_position('none')
		ax.set_xticks([])
		ax.set_yticks([])

		fig.delaxes(axes[pnum][axn['del']])
		for ax in [axes[pnum][axn[i]] for i in ['x','y','contour'][-1:]]:
			ax.set_xticks(angle_ticks)
			ax.set_yticks(angle_ticks)

	patches,counter = [],0
	tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
	for ss,sn in enumerate(sns_special):
		pnum = [1,3][ss]
		for posnum,pos in enumerate(positions):
			fn = os.path.join(tempdir,snapshot_namer(sn=work.raw.prefixer(sn),tag=tag,pos=pos)+'.png')
			im = mpl.image.imread(fn)
			ax = axes[pnum][posnum]
			ax.imshow(im)
			ax.axis('off')
			tb = ax.text(0.5,0,chr(ord('A')+counter),
				bbox=tagbox,rotation=0,ha="center",va="top",color='k',
				fontsize=fsbase,transform=ax.transAxes)
			patches.append(tb)
			if posnum == 0:
				tb = ax.text(0.1,0.9,work.meta[sn]['ion_label'],
					bbox=tagbox,rotation=0,ha="center",va="center",color='k',
					fontsize=fsbase+4,transform=ax.transAxes,zorder=10)
				patches.append(tb)
			#---add the marker to the histogram
			ax = axes[[0,2][ss]][axn['contour']]
			theta,phi = [reps[(sn,pos)][i] for i in ['theta','phi']]
			tb = ax.text(theta,phi,chr(ord('A')+counter),
				rotation=0,ha="center",va="center",color='k',
				fontsize=fsbase-4)
			tb.set_path_effects([path_effects.Stroke(linewidth=4,foreground='w'),
				path_effects.Normal()])
			#bb = tb.get_bbox_patch()
			#bb.set_boxstyle("round,pad=0.1")
			counter += 1

	#---saving the snapshot tag here so we can keep track of the fix to exemplar_lipid above, fixed on v4
	picturesave('fig.head_angle_detail',work.plotdir,backup=False,version=True,
		meta={'tag':tag_head_angle,'exemplar_rank':exemplar_rank},extras=patches)
	plt.close()

def render_lipid_pair(**kwargs):
	"""
	Abstracted the vmdmake routine for looping.
	"""
	#---default settings
	hbond_style = ['dashed','cylinder'][-1]
	licorice_thick = 0.4
	hbond_cylinder_radius = 0.15
	goodsell_ions = False
	goodsell_lipids = True
	#---unpack the kwargs
	sn = kwargs.pop('sn')
	frame = kwargs.pop('frame')
	tempdir = kwargs.pop('tempdir')
	hydrogen_bonds = kwargs.pop('hydrogen_bonds',{})
	pair = kwargs.pop('pair',{})
	snapshot_fn = kwargs.pop('fn')
	associates = kwargs.pop('associates',None)
	if not pair: raise Exception('render_lipid_pair requires a pair')
	gro,xtc,tpr = [kwargs.pop(i) for i in 'gro xtc tpr'.split()]
	if kwargs: raise Exception('unprocessed arguments: %s'%kwargs)
	if False:
		index = kwargs['index']
		#---! backwards compatibility
		rn1,ri1,an1,rn2,ri2,an2 = [kwargs[key] for key in 
			'resname_1 resid_1 name_1 resname_2 resid_1 name_2'.split()]
		code = rn1,ri1,an1,hname2,rn2,ri2,an2
		snapshot_fn = snapshot_namer(tag=tag,sn=sn_this,index=index,
			code='_'.join(code),rank=rank)
		gro = work.slice(sn_this)['current']['all']['gro']
		xtc = work.slice(sn_this)['current']['all']['xtc']
		tpr = work.get_last_start_structure(sn_this,part_name='tpr')
	#---run VMD to make the snapshots via vmdmake
	view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
		frames='',res=(4000,4000),**kwargs_vmdmake)
	view.do('load_dynamic','standard','bonder')
	view.command("animate goto %d"%(frame))
	selections = []
	#---show both lipids
	selections.append("(resname %s and resid %s) and not hydrogen"%(pair['resname_1'],pair['resid_1']))
	view.select(**{'partner_1':selections[-1],'style':
		'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids})
	selections.append("(resname %s and resid %s) and not hydrogen"%(pair['resname_2'],pair['resid_2']))
	view.select(**{'partner_2':selections[-1],
		'style':'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids})
	#---show associates
	if associates:
		for assoc_num,assoc in enumerate(associates):
			selections.append("(resname %s and resid %s) and not hydrogen"%(assoc['resname'],assoc['resid']))
			view.select(**{'assoc_%d'%assoc_num:selections[-1],'style':
				'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids})
	#---for each hydrogen bond, show the lipids with the hydrogen (both are plotted to include the bond)
	this_style = {'style':'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids}
	for hnum,hbond in enumerate(hydrogen_bonds):
		#---the first residue in the pair is the donor and we draw it here
		selstring = "(resname %s and resid %s)"%(hbond['resname_1'],hbond['resid_1'])
		selections.append('(%s and not hydrogen) or (%s and name %s)'%(selstring,selstring,hbond['name_h']))
		select_args = dict(**{'partner_1_hbond_%d'%hnum:selections[-1]})
		select_args.update(**this_style)
		view.select(**select_args)
		#---also set the endpoints for the bond itself
		#---! currently set for heavy-heavy dashed line below
		#view.select(**{'partner_1_hbond_%d_ref'%hnum:'(%s and name %s)'%(selstring,hbond['name_h'])})
		view.command('set partner_1_hbond_%d_ref [atomselect top "%s"]'%(
			hnum,'(%s and name %s)'%(selstring,hbond['name_1'])))
		#---draw the other lipid
		selstring = "(resname %s and resid %s)"%(hbond['resname_2'],hbond['resid_2'])
		selections.append('(%s and not hydrogen)'%selstring)
		select_args = dict(**{'partner_2_hbond_%d'%hnum:selections[-1]})
		select_args.update(**this_style)
		view.select(**select_args)
		#view.select(**{'partner_2_hbond_%d_ref'%hnum:'(%s and name %s)'%(selstring,hbond['name_2'])})
		view.command('set partner_2_hbond_%d_ref [atomselect top "%s"]'%(
			hnum,'(%s and name %s)'%(selstring,hbond['name_2'])))
	#---beads for the acceptor/donor
	if False:
		selections.append("((resname %s and resid %s and name %s) and not hydrogen)"%(rn1,ri1,an1))
		view.select(**{
			'partner_1_bead':selections[-1],
			'style':'VDW 0.400000 12.000000','goodsell':True})
		#---beads for the acceptor/donor
		selections.append("((resname %s and resid %s and name %s) and not hydrogen)"%(rn2,ri2,an2))
		view.select(**{
			'partner_2_bead':selections[-1],
			'style':'VDW 0.400000 12.000000','goodsell':True})
	if False:
		#---bead for the hydrogen (I think residue 2 is always the donor but this makes it sure)
		selections.append("(resid %s or resid %s) and name %s"%(ri1,ri2,hname2))
	#---set the selections
	view.command('set me [atomselect top "%s"]'%' or '.join(['(%s)'%i for i in selections]))
	if False:
		view.select(**{
			'hydrogen_bead':selections[-1],
			'style':'VDW 0.400000 12.000000','goodsell':True})
	#---move everything
	#---draw a line for the hydrogen bond because the hydrogen bond rep is too permissive 
	#---...and dashed is hard to see
	#---below is a method for hydrogen bonding but VMD is too permissive 
	#---...and allows carbon to be the donor for whatever reason. perhaps it can be modified.
	if False:
		#---show both in case you want to do a hydrogen bond representation
		view.select(**{'both_partners':
			"((resname %s and resid %s)) or ((resname %s and resid %s))"%(rn1,ri1,rn2,ri2),
			'style':'HBonds 4.000000 30.000000 3.000000'})
		view.command("mol modselect 2 0 ((resname %s and resid %s)) or ((resname %s and resid %s))"%
			(rn1,ri1,rn2,ri2))
		view.command("mol modstyle 2 0 HBonds 3.000000 20.000000 4.000000")
		view.command("mol modcolor 2 top ColorID 16")
	#---set selections for the vmd_render_good_side_bond method over a bond (or at least two points)
	#---we could use representative points like "C14" or "P" but the conformations are too dynamic
	for resid_ind in [1,2]:
		view.command("set %s [atomselect top \"(resname %s and resid %s and name %s)\"]"%(
			'partner_%d_bead'%resid_ind,pair['resname_%d'%resid_ind],
			pair['resid_%d'%resid_ind],pair['name_%d'%resid_ind]))
	#---show the good side
	view.command(vmd_render_good_side_bond)
	#---see vmdmake.py for vmd colors
	cation = work.meta[sn].get('cation_relevant',work.meta[sn]['cation'])
	view.set_color_cursor({'MG':'pink','Cal':'blue','NA':'green','K':'tan'}.get(cation,'purple'))
	#---show ions near both residues
	for resid_ind in [1,2]:
		view.select(**{'near_ions_%d'%resid_ind:
			"name %s and within 4.6 of (resname %s and resid %s)"%(
			cation,work.meta[sn]['ptdins_resname'],pair['resid_%d'%resid_ind]),
			'smooth':False,'style':'VDW 0.4 12.0','goodsell':goodsell_ions,'color_specific':True})
	#---emphazize the hydrogen bonds
	for hnum,hbond in enumerate(hydrogen_bonds):
		if hbond_style=='cylinder':
			view.command('draw color black')
			view.command(("draw cylinder [expr [$partner_1_hbond_%d_ref get {x y z}]] "+
				"[expr [$partner_2_hbond_%d_ref get {x y z}]] radius %.3f")%(
				hnum,hnum,hbond_cylinder_radius))
		elif hbond_style=='dashed':
			view.command(("draw line [expr [$partner_1_hbond_%d_ref get {x y z}]] "+
				"[expr [$partner_2_hbond_%d_ref get {x y z}]] style dashed")%(hnum,hnum))
		else: raise Exception('unclear hydrogen bond style %s'%hbond_style)
	#---render to disk
	view['snapshot_filename'] = snapshot_fn
	view.do('snapshot')
	view.show(quit=True)

if 'common_hydrogen_bonding' in routine:

	#---settings
	tag_hbonds = 'v5'
	overwrite_snaps = True
	sns_group = list(sns)
	#---! override because tired of dealing with upside-down lipids!
	sns_group = [sn for sn in work.sns() if work.meta[sn]['composition_name']=='asymmetric']
	n_most_popped_bonds = 3
	#---number of snapshots to do. we walk down the list of ranked commonest bonds, which is 
	#---...partly arbitrary (see notes below) hence  this just gives us more snapshots more or less
	#---you can specify explicit nranked list to debug things very carefully
	nranked = 10 # or try a list e.g. [0]

	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],'fig.hydrogen_bonding.%s'%tag_hbonds)
	#---only make figures if the folder is empty. it only takes a few minutes -- no more coding logic
	if (os.path.isdir(tempdir) and not overwrite_snaps): status('found %s and refusing to overwrite'%tempdir)
	else:
		try: os.mkdir(tempdir)
		except: status('directory exists so we are overwriting')
		status('snapshots dropping to %s'%tempdir,tag='note')
		if not overwrite_snaps: status('delete the snapshot directory to rebuild them')
		#---we no longer precache identities here
		#---note that this section has been ported in from plot-hydrogen_bonding.py from legacy
		#---...omnicalc however the underlying hydrogen bond data structure changed a lot since it was
		#---...written so we collect the most populated hydrogen bonds using some new code below.
		#---the original code was taking the most common specific bond including resid, resname, atom
		#---...however in this version we are taking the most common resname pairs. this reproduces 
		#---...something similar to the original snapshotter, however it is also not composition-
		#---...weighted, and in later sections we will use better selection criteria
		#---the original code was from plot-hydrogen_bonds.py and plot-head_angle_contours.py in legacy
		snapshot_catalog = {}
		#---get the top three hydrogen bonds from each simulation
		for sn in sns_group:
			bonds,obs,valid_frames = [data['hydrogen_bonding'][sn]['data'][i] 
				for i in ['bonds','observations','valid_frames']]
			#---discard atom names
			bonds_red = bonds[:,np.array([0,1,3,4])]
			#---this step induces a sort of sorts, so that actually the line where we define the ranking
			#---...ends up being in order over the mean of the sum of observations
			bonds_inds,bonds_counts = uniquify(bonds_red)
			#---subselect lipids and filter out intramolecular bonds
			lipid_resnames = np.array(work.vars['selectors']['resnames_lipid_chol'])
			#---further subselection
			lipid_resnames = np.array([work.meta[sn]['ptdins_resname']])
			subsel = np.where(np.all((
				np.in1d(bonds_red[bonds_inds][:,0],lipid_resnames),
				np.in1d(bonds_red[bonds_inds][:,2],lipid_resnames),
				bonds_red[bonds_inds][:,1]!=bonds_red[bonds_inds][:,3]
				),axis=0))[0]
			#---having filtered for desired bonds and already uniquified, we find the most common
			#---...after we map each unique reduced bond back to all of the actual bonds (with atoms)
			#---note that this step is slow
			status('slow where loop to map reduced bonds to observations for %s'%sn,tag='compute')
			map_red_to_obs = [np.where(np.all(bonds_red==i,axis=1))[0] 
				for i in bonds_red[bonds_inds][subsel]]
			#---loop over rankings. these are somewhat arbitrary because we have subselected so hard
			ranking = np.argsort([obs.T[i].sum(axis=0).mean() for i in map_red_to_obs])[::-1]
			nranked = range(nranked) if type(nranked)==int else nranked
			for rank_num in nranked:
				#---now we have the commonest bonds so we can find the most-representative frame
				#---the ranking is over the map_red_to_obs which indexes the whole list
				#---...so we consult the whole obs list to get the particular bonds that have the reduced
				#---...bond and then identify the frame with the max, noting however that many frames
				#---...will probably have the max since our unique residue-residue bond probably only has
				#---...a handful of different unique bonds over the atoms
				fr = np.argmax(obs.T[map_red_to_obs[ranking[rank_num]]].sum(axis=0))
				#---! wow. where to begin. this frame is wrong and it was really tough to figure out why
				frame_actual = valid_frames[fr]
				#---for each of the observed hydrogen bonds at this frame we save instructions to render it
				#---first we get the handful of bonds that for our resid-resid pair
				bonds_by_pair = bonds[map_red_to_obs[ranking[rank_num]]][
					np.where(obs.T[map_red_to_obs[ranking[rank_num]]].sum(axis=1))]
				#---then at the frame we selected (admittedly probably many frames have a couple bonds)
				#---...we figure out which of the handful are actually on that frame
				bonds_in_frame = np.where(obs.T[map_red_to_obs[ranking[rank_num]]].T[fr])
				bond_spec_keys = 'resname_1 resid_1 name_1 resname_2 resid_2 name_2 name_h'
				bond_spec = dict(hydrogen_bonds=[
					dict([(i,bond[ii]) for ii,i in enumerate(bond_spec_keys.split())]) 
					for bond in bonds_by_pair[bonds_in_frame]])
				#---construct a description of the lipid pair, first checking that there is only one pair
				inds,counts = uniquify(bonds_by_pair[bonds_in_frame][:,np.array([0,1,3,4])])
				if not len(inds)==1: raise Exception('non-unique lipid pairing: %s'%
					bonds_by_pair[bonds_in_frame])
				lipid_pair_spec = dict([(i,bonds_by_pair[bonds_in_frame][inds[0]][ii]) 
					for ii,i in enumerate(bond_spec_keys.split())])
				#---now that we have all of the instructions we render it		
				filetag = 'fig.snapshot.%s.fr%d.%s_%s_o%d'%(sn,fr,
					lipid_pair_spec['resid_1'],lipid_pair_spec['resid_2'],rank_num)
				render_spec = dict(sn=sn,pair=lipid_pair_spec,frame=frame_actual,tempdir=tempdir,fn=filetag)
				render_spec.update(**work.get_gmx_sources(sn=sn,calc=calc))
				render_spec.update(**bond_spec)
				#---check the lipid distance to see if we are broken across PBCs
				resids = [np.where(data['lipid_abstractor'][sn]['data']['resids']==int(
					bond_spec['hydrogen_bonds'][0]['resid_%d'%i]))[0][0] for i in [1,2]]
				com_distance = data['lipid_abstractor'][sn]['data']['points'][
					frame_actual][resids].ptp(axis=0)
				if np.any(com_distance>=data['lipid_abstractor'][sn]['data']['vecs'][frame_actual]/2.):
					status('lipids are broken over PBCs and this rendering algo cannot work so skipping',
						tag='warning')
					continue
				import ipdb;ipdb.set_trace()
				render_lipid_pair(**render_spec)
				status('done rendering simulation %s'%sn,i=rank_num,looplen=len(nranked),tag='render')
				snapshot_catalog[(sn,rank_num)] = render_spec

if 'common_hydrogen_bonding_specific' in routine:

	overwrite_snaps = True
	sns_group = list(sns)
	#---! override because tired of dealing with upside-down lipids!
	sns_group = [sn for sn in work.sns() if work.meta[sn]['composition_name']=='asymmetric']
	#---find pairs of a particular type that are the most common
	trawler = ('ptdins','ptdins')
	trawler = ('ptdins','CHL1')
	nranked = 10 # or try a list e.g. [0]
	tag_hbonds = 'v5_%s_%s'%(trawler)
	#---! eventually we need a loop over the trawler, tag, and nranked sections
	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],'fig.hydrogen_bonding.%s'%tag_hbonds)
	#---only make figures if the folder is empty. it only takes a few minutes -- no more coding logic
	if (os.path.isdir(tempdir) and not overwrite_snaps): status('found %s and refusing to overwrite'%tempdir)
	else:
		try: os.mkdir(tempdir)
		except: status('directory exists so we are overwriting')
		status('snapshots dropping to %s'%tempdir,tag='note')
		if not overwrite_snaps: status('delete the snapshot directory to rebuild them')
		bond_spec_keys = 'resname_1 resid_1 name_1 resname_2 resid_2 name_2 name_h'
		for sn in sns_group:
			dat = data['hydrogen_bonding'][sn]['data']
			resnames = dict(zip(dat['resnames'],dat['resnames']))
			resnames['ptdins'] = work.meta[sn]['ptdins_resname']
			resname_pair = [resnames[j] for j in trawler]
			#---start with the complete set of bonds
			bonds = dat['bonds']
			obs = dat['observations']
			#---first filter for intermolecular bonds
			bonds_inter = bonds[:,1]!=bonds[:,4]
			sel_inter = np.where(bonds_inter)
			bondlist = bonds[sel_inter]
			#---get unique intermolecular bonds
			#---! THEY ARE ALREADY UNIQUE. SWITCHING TO A NEW, FIFTH METHOD WITH A NEW TAG
			bondlist_inds,bondlist_counts = uniquify(bondlist)
			#---filter for bonds where the two residue names are set above by trawler
			sel_resname_match = np.where(np.any([np.all((bondlist[bondlist_inds][:,
				np.array([0,3])]==resname_pair[::i]),axis=1) for i in [1,-1]],axis=0))[0]
			#---sort bonds by average occupancy
			resort_by_mean = obs.mean(axis=0)[bondlist_inds[sel_resname_match]].argsort()[::-1]
			modal_bond_inds = bondlist_inds[sel_resname_match][resort_by_mean]
			#---loop over rankings
			nranked = range(nranked) if type(nranked)==int else nranked
			for rank_num in nranked[:len(modal_bond_inds)]:
				modal_bond_ind = modal_bond_inds[rank_num]
				modal_bond = bondlist[modal_bond_ind]
				#---now that we have the modal bond we can find the modal instance of the bond in the list
				modal_bond_pair = modal_bond[np.array([0,2])],modal_bond[np.array([3,5])]
				matched_bonds = np.where(np.any([np.all([np.all(
					bondlist[bondlist_inds][:,np.array(i)]==modal_bond_pair[::k][j],axis=1) 
					for i,j in [([0,2],0),([3,5],1)]],axis=0) for k in [1,-1]],axis=0))[0]
				modal_matched_bonds_ind = np.argsort(obs.sum(axis=0)[
					sel_inter][bondlist_inds][matched_bonds])[-1]
				#---now we choose the actual observed bond that is the commonest
				modal_bond = bonds[sel_inter][bondlist_inds][matched_bonds][modal_matched_bonds_ind]
				#---this bond is actually kind of sparse so now we check for total bonds between 
				#---...residues on this list
				bonds_sub = bonds[sel_inter][bondlist_inds][sel_resname_match]
				#---we could use uniquify but we want both directions so we sort first
				modal_resname_pair_ind = uniquify(np.array([sorted(i) 
					for i in bonds_sub[:,np.array([1,4])]]))[0][0]
				#---probably not a coincidence that this is usually zero, since we are sorting this upstream
				#---get the resid pair for the commonst bond between all residues that match trawler 
				#---...regardless of atom
				modal_resname_pair = bonds_sub[modal_resname_pair_ind,np.array([1,4])]
				sel_resid_pair = np.where(np.any([np.all(
					[bondlist[bondlist_inds][sel_resname_match][:,i]==modal_resname_pair[::k][j] 
					for i,j in [(1,0),(4,1)]],axis=0) for k in [1,-1]],axis=0))[0]
				#---get the frame where these two are bound the most
				fr = np.argmax(obs.T[sel_inter][bondlist_inds][sel_resname_match][
					sel_resid_pair].T.sum(axis=1))
				import ipdb;ipdb.set_trace()
				#---which two bonds we can expect at this frame (note there may be other frames with equal 
				#---...multiplicities)
				which_these_bonds_inds = obs.T[sel_inter][bondlist_inds][
					sel_resname_match][sel_resid_pair].T[fr]
				which_these_bonds = bondlist[bondlist_inds][sel_resname_match][
					sel_resid_pair][np.where(which_these_bonds_inds)]
				#---! continuing development after a while
				#---prepare the pair specification
				valid_frames = data['hydrogen_bonding'][sn]['data']['valid_frames']
				frame_actual = valid_frames[fr]
				#---only one of the bonds is necessary to specify the lipids
				lipid_pair_spec = dict([(i,which_these_bonds[0][ii]) 
					for ii,i in enumerate(bond_spec_keys.split())])
				filetag = 'fig.snapshot.%s.fr%d.%s_%s_o%d'%(sn,fr,
					lipid_pair_spec['resid_1'],lipid_pair_spec['resid_2'],rank_num)
				render_spec = dict(sn=sn,pair=lipid_pair_spec,frame=frame_actual,tempdir=tempdir,fn=filetag)
				render_spec.update(**work.get_gmx_sources(sn=sn,calc=calc))
				bond_spec = dict(hydrogen_bonds=[dict([(i,bond[ii]) 
					for ii,i in enumerate(bond_spec_keys.split())]) 
					for bond in which_these_bonds])
				render_spec.update(**bond_spec)
				#---check the lipid distance to see if we are broken across PBCs
				resids = [np.where(data['lipid_abstractor'][sn]['data']['resids']==int(
					bond_spec['hydrogen_bonds'][0]['resid_%d'%i]))[0][0] for i in [1,2]]
				com_distance = data['lipid_abstractor'][sn]['data']['points'][
					frame_actual][resids].ptp(axis=0)
				if np.any(com_distance>=data['lipid_abstractor'][sn]['data']['vecs'][frame_actual]/2.):
					status('lipids are broken over PBCs and this rendering algo cannot work so skipping',
						tag='warning')
					#continue
				render_lipid_pair(**render_spec)

if 'consolidated_hydrogen_bonding' in routine:

	#---extra settings
	overwrite_snaps = False
	#---! override because tired of dealing with upside-down lipids!
	sns_group = [sn for sn in work.sns() if work.meta[sn]['composition_name']=='asymmetric']

	#---specify a tag for each method
	methods = {
		#'v6_common':{'nranked':10,'sns':['membrane-v532']},
		#'v6_ptdins_chl1':{'nranked':10,'resname_filter':['PtdIns','CHL1'],'sns':['membrane-v532']},
		'v6_ptdins_chl1_with_chl1':{'nranked':10,'resname_filter':['PtdIns','CHL1'],
			'sns':['membrane-v532'],'barnacles':['CHL1'],'associates':['CHL1']},
		'v6_ptdins_with_chl1':{'nranked':10,'resname_filter':['PtdIns'],
			'sns':['membrane-v532'],'barnacles':['CHL1'],'associates':['CHL1']},
		}

	#---one rendering run per method
	for method_tag,method in methods.items():

		#---store the snapshots in the post_plot_spot according to the tag
		tempdir = os.path.join(work.paths['post_plot_spot'],'fig.hydrogen_bonding.%s'%method_tag)
		#---only make figures if the folder is empty. it only takes a few minutes -- no more coding logic
		if (os.path.isdir(tempdir) and not overwrite_snaps): 
			status('found %s and refusing to overwrite'%tempdir)
			continue
		try: os.mkdir(tempdir)
		except: status('directory exists so we are overwriting',tag='warning')
		status('snapshots dropping to %s'%tempdir,tag='note')
		if not overwrite_snaps: status('delete the snapshot directory to rebuild them')

		snapshot_catalog = {}
		#---get the top three hydrogen bonds from each simulation
		for sn in method.get('sns',work.sns()):
			bonds,obs,valid_frames = [data['hydrogen_bonding'][sn]['data'][i] 
				for i in ['bonds','observations','valid_frames']]
			#---discard atom names
			bonds_red = bonds[:,np.array([0,1,3,4])]
			#---this step induces a sort of sorts, so that actually the line where we define the ranking
			#---...ends up being in order over the mean of the sum of observations
			bonds_inds,bonds_counts = uniquify(bonds_red)
			#---subselect lipids and filter out intramolecular bonds
			lipid_resnames = np.array(work.vars['selectors']['resnames_lipid_chol'])
			#---further subselection
			resname_filter = method.get('resname_filter',None)
			if resname_filter: 
				lipid_resnames = [l for l in lipid_resnames if l in resname_filter]
				#---alias "PtdIns" back to the right residue name
				lipid_resnames = [work.meta[sn]['ptdins_resname'] if l=='PtdIns' 
					else l for l in resname_filter]
			nranked = method['nranked']
			#---perform the selection
			subsel = np.where(np.all((
				np.in1d(bonds_red[bonds_inds][:,0],lipid_resnames),
				np.in1d(bonds_red[bonds_inds][:,2],lipid_resnames),
				bonds_red[bonds_inds][:,1]!=bonds_red[bonds_inds][:,3]
				),axis=0))[0]
			#---having filtered for desired bonds and already uniquified, we find the most common
			#---...after we map each unique reduced bond back to all of the actual bonds (with atoms)
			#---note that this step is slow
			status('slow where loop to map reduced bonds to observations for %s'%sn,tag='compute')
			map_red_to_obs = [np.where(np.all(bonds_red==i,axis=1))[0] 
				for i in bonds_red[bonds_inds][subsel]]
			#---loop over rankings. these are somewhat arbitrary because we have subselected so hard
			ranking = np.argsort([obs.T[i].sum(axis=0).mean() for i in map_red_to_obs])[::-1]
			nranked = range(nranked) if type(nranked)==int else nranked
			for rank_num in nranked:
				#---now we have the commonest bonds so we can find the most-representative frame
				#---the ranking is over the map_red_to_obs which indexes the whole list
				#---...so we consult the whole obs list to get the particular bonds that have the reduced
				#---...bond and then identify the frame with the max, noting however that many frames
				#---...will probably have the max since our unique residue-residue bond probably only has
				#---...a handful of different unique bonds over the atoms
				fr = np.argmax(obs.T[map_red_to_obs[ranking[rank_num]]].sum(axis=0))
				#---! wow. where to begin. this frame is wrong and it was really tough to figure out why
				frame_actual = valid_frames[fr]
				#---for each of the observed hydrogen bonds at this frame we save instructions to render it
				#---first we get the handful of bonds that for our resid-resid pair
				bonds_by_pair = bonds[map_red_to_obs[ranking[rank_num]]][
					np.where(obs.T[map_red_to_obs[ranking[rank_num]]].sum(axis=1))]
				#---then at the frame we selected (admittedly probably many frames have a couple bonds)
				#---...we figure out which of the handful are actually on that frame
				bonds_in_frame = np.where(obs.T[map_red_to_obs[ranking[rank_num]]].T[fr])
				bond_spec_keys = 'resname_1 resid_1 name_1 resname_2 resid_2 name_2 name_h'
				bond_spec = dict(hydrogen_bonds=[
					dict([(i,bond[ii]) for ii,i in enumerate(bond_spec_keys.split())]) 
					for bond in bonds_by_pair[bonds_in_frame]])
				#---construct a description of the lipid pair, first checking that there is only one pair
				inds,counts = uniquify(bonds_by_pair[bonds_in_frame][:,np.array([0,1,3,4])])
				if not len(inds)==1: raise Exception('non-unique lipid pairing: %s'%
					bonds_by_pair[bonds_in_frame])
				lipid_pair_spec = dict([(i,bonds_by_pair[bonds_in_frame][inds[0]][ii]) 
					for ii,i in enumerate(bond_spec_keys.split())])
				#---now that we have all of the instructions we render it		
				filetag = 'fig.snapshot.%s.fr%d.%s_%s_o%d'%(sn,frame_actual,
					lipid_pair_spec['resid_1'],lipid_pair_spec['resid_2'],rank_num)
				render_spec = dict(sn=sn,pair=lipid_pair_spec,frame=frame_actual,tempdir=tempdir,fn=filetag)
				render_spec.update(**work.get_gmx_sources(sn=sn,calc=calc))
				render_spec.update(**bond_spec)
				#---check the lipid distance to see if we are broken across PBCs
				resids = [np.where(data['lipid_abstractor'][sn]['data']['resids']==int(
					bond_spec['hydrogen_bonds'][0]['resid_%d'%i]))[0][0] for i in [1,2]]
				com_distance = data['lipid_abstractor'][sn]['data']['points'][
					frame_actual][resids].ptp(axis=0)
				if np.any(com_distance>=data['lipid_abstractor'][sn]['data']['vecs'][frame_actual]/2.):
					status('lipids are broken over PBCs and this rendering algo cannot work so skipping',
						tag='warning')
					continue
				#---identify the barnacles: lipids nearby that we care about and which also appear in the 
				#---...listing of hydrogen bonds
				barnacles = method.get('barnacles',None)
				if barnacles:
					this_bond = bonds_by_pair[0]
					resname_1,resname_2,resid_1,resid_2 = [this_bond[bond_spec_keys.split().index(i)] for i in ['resname_1','resname_2','resid_1','resid_2']]
					subsel_barnacles = np.where(
						np.any((
							np.all((np.in1d(bonds_red[bonds_inds][:,0],[resname_1,resname_2]),np.in1d(bonds_red[bonds_inds][:,1],[resid_1,resid_2]),np.in1d(bonds_red[bonds_inds][:,2],method['barnacles']),),axis=0),
							np.all((np.in1d(bonds_red[bonds_inds][:,2],[resname_1,resname_2]),np.in1d(bonds_red[bonds_inds][:,3],[resid_1,resid_2]),np.in1d(bonds_red[bonds_inds][:,0],method['barnacles']),),axis=0),
							),axis=0)
						)[0]
					#---only render the bonds if the barnacle has a bond in this frame
					for subind in subsel_barnacles:
						if obs.T[bonds_inds][subind][fr]==1.0:
							bond = bonds[bonds_inds][subind]
							render_spec['hydrogen_bonds'].append(dict([(i,bond[ii]) 
								for ii,i in enumerate(bond_spec_keys.split())]))
				#---check for associates
				associates = method.get('associates',None)
				if associates:
					simplices = data['lipid_mesh'][sn]['data']['0.%d.simplices'%frame_actual]
					resname_resids = np.transpose((data['lipid_abstractor'][sn]['data']['resnames'],data['lipid_abstractor'][sn]['data']['resids']))
					#---loop over partners
					for xxx in [1,4]:
						#---get resid for the first residue in the pair
						ind = list(resname_resids[:,1]).index(bonds_by_pair[0][xxx])
						#---reverse mapping to monolayers
						#---! note this was extremely hasty
						indm = np.where(np.where(data['lipid_abstractor'][sn]['data']['monolayer_indices']==0)[0]==ind)[0][0]
						neighborsm = np.unique(simplices[np.where(np.any(simplices==indm,axis=1))[0]])
						ghosts = data['lipid_mesh'][sn]['data']['0.%d.ghost_ids'%frame_actual][neighborsm]
						neighbors = np.where(data['lipid_abstractor'][sn]['data']['monolayer_indices']==0)[0][ghosts]
						near_neighbors = np.where(np.in1d(data['lipid_abstractor'][sn]['data']['resnames'][neighbors],associates))[0]
						data['lipid_mesh'][sn]['data'].keys()
						if len(near_neighbors)>0: render_spec['associates'] = []
						#---include all nearest neighbors in the visualization
						for resid in neighbors[near_neighbors]:
							render_spec['associates'].append(dict(
								resname=data['lipid_abstractor'][sn]['data']['resnames'][resid],
								resid=str(data['lipid_abstractor'][sn]['data']['resids'][resid])))
				#---also identify 
				render_lipid_pair(**render_spec)
				status('done rendering simulation %s'%sn,i=rank_num,looplen=len(nranked),tag='render')
				snapshot_catalog[(sn,rank_num)] = render_spec
				#if rank_num==3: 
				#	import ipdb;ipdb.set_trace()
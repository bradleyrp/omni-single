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
from base.tools import dictsum
import matplotlib.patheffects as path_effects
#---! for some reason we are importing post twice on this plot. because two spots?
#---one-off settings for multiple sections: head_angle_snaps and head_angle_snaps_detail
#---note that RPB has modified cg_bonds.tcl from the GMX 5 version to include SETTLE bonds in water (on dark)
kwargs_vmdmake = {'CGBONDSPATH':'~/libs/cg_bonds.tcl','GMXDUMP':'/usr/local/gromacs-5.0.4/bin/gmx'}
bond_spec_keys = 'resname_1 resid_1 name_1 resname_2 resid_2 name_2 name_h'

routine = [
	#---head_snaps_detail depends on head_angle_snaps
	'head_angle_snaps','head_angle_snaps_detail',
	#---this hydrogen bonding snapshot code works great and generate many images
	'consolidated_hydrogen_bonding',
	#---attempt to separate the identification of relevant pairs from the plot aesthetics
	'systematic_snapshotter','arrange_snapshots'][-2:-1]

if 'data' not in globals():
	data,calc = plotload(plotname)
	sns = work.sns()
	#---removed head-angle data here. also removed salt_bridges
	#---considered using the lipid_abstractor or mesh objects to check proximity to the box edge
	#---use a global variable for large objects
	memory = {}

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
			counter += 1

	#---saving the snapshot tag here so we can keep track of the fix to exemplar_lipid above, fixed on v4
	picturesave('fig.head_angle_detail',work.plotdir,backup=False,version=True,
		meta={'tag':tag_head_angle,'exemplar_rank':exemplar_rank},extras=patches)
	plt.close()

###---CODE REORGY

import makeface
#---get an automacs landscape with a little help from the user
try: mod = makeface.import_remote('amx/amx')
except: raise Exception('please clone a copy of automacs next to omni in `amx`. '
	'you must also run `make setup all` from that directory to get force field files.')
mod['state'].force_field = 'charmm'
GMXStructure = mod['GMXStructure']

def mapback(seq):
	"""Hash a list of numbers back to their indices."""
	#---this function is a candidate for omni/base/tools.py
	return dict([(v,k) for k,v in zip(np.arange(len(seq)),seq)])

#---custom water coloring
custom_water_coloring = """
set vdw_thick 0.4
set inner_cutoff 2.0
set outer_cutoff 5.0
set counter [molinfo top get numreps]
foreach pivot_ind [$pivots list] {
set pivot [atomselect top "([$pivots text]) and (index $pivot_ind)"]
mol selection [$pivot text]
mol addrep top
mol modstyle $counter 0 VDW $vdw_thick 12.000000
incr counter
set com_pivot [measure center $pivot]
set selector [atomselect top "same residue as (name OW and within 4.6 of ([$pivot text]))"]
foreach ind1 [$selector list] {
set distal [atomselect top "water and same residue as index $ind1"]
set com_distal [measure center $distal]
set val [expr ([veclength [vecsub $com_distal $com_pivot]] - $inner_cutoff ) / ( $outer_cutoff - $inner_cutoff )]
set newsel [atomselect top "water and same residue as index $ind1"]
mol selection [$newsel text]
mol addrep top
mol modstyle $counter 0 Licorice 0.300000 10.000000 10.000000
$newsel set user $val
mol modcolor $counter 0 User
incr counter
}}
"""

def render_lipid_pair(**kwargs):
	"""
	Abstracted the vmdmake routine for looping.
	"""
	global special_ion_color
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
	lipid_color_specific = kwargs.pop('lipid_color_specific',None)
	if not pair: raise Exception('render_lipid_pair requires a pair')
	gro,xtc,tpr = [kwargs.pop(i) for i in 'gro xtc tpr'.split()]
	lipid_atoms = kwargs.pop('lipid_atoms',None)
	lipid_style = kwargs.pop('lipid_style',{})
	show_water = kwargs.pop('show_water',False)
	lipid_colors = kwargs.pop('lipid_colors',False)
	if kwargs: raise Exception('unprocessed arguments: %s'%kwargs)
	#---run VMD to make the snapshots via vmdmake
	view = vmdmake.VMDWrap(site=tempdir,gro=gro,xtc=xtc,tpr=tpr,
		frames='',res=(4000,4000),**kwargs_vmdmake)
	view.do('load_dynamic','standard','bonder')
	view.command("animate goto %d"%(frame))
	selections = []
	#---change water oxygen color for show_water (starting with v15)
	view.command('color Name O pink')
	#---show both lipids
	sel1 = "(resname %s and resid %s) and not hydrogen"%(pair['resname_1'],pair['resid_1'])
	sel2 = "(resname %s and resid %s) and not hydrogen"%(pair['resname_2'],pair['resid_2'])
	if lipid_atoms is not None:
		sel1 = '(%s) and (%s)'%(sel1,' or '.join(['name %s'%i for i in lipid_atoms]))
		sel2 = '(%s) and (%s)'%(sel2,' or '.join(['name %s'%i for i in lipid_atoms]))
	selections.append(sel1)
	#---disallow the two lipid coloring schemes to work at the same time
	if lipid_color_specific and lipid_colors:
		raise Exception('incompatible arguments: lipid_color_specific and lipid_colors')
	if lipid_color_specific: view.set_color_cursor(lipid_color_specific)
	if lipid_colors: view.set_color_cursor(lipid_colors[pair['resname_1']])
	view.select(**dictsum({'partner_1':selections[-1],
		'style':'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids},lipid_style,
		{'color_specific':True} if lipid_color_specific or lipid_colors else {}))
	selections.append(sel2)
	if lipid_colors: view.set_color_cursor(lipid_colors[pair['resname_2']])
	view.select(**dictsum({'partner_2':selections[-1],
		'style':'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids},lipid_style,
		{'color_specific':True} if lipid_color_specific or lipid_colors else {}))
	#---show associates
	if associates:
		for assoc_num,assoc in enumerate(associates):
			if lipid_colors: view.set_color_cursor(lipid_colors[assoc['resname']])
			selections.append("(resname %s and resid %s) and not hydrogen"%(assoc['resname'],assoc['resid']))
			view.select(**dictsum({'assoc_%d'%assoc_num:selections[-1],'style':
				'Licorice %.3f 12.000000 12.000000'%licorice_thick,'goodsell':goodsell_lipids},
				{'color_specific':True} if lipid_color_specific or lipid_colors else {}))
	#---for each hydrogen bond, show the lipids with the hydrogen (both are plotted to include the bond)
	if lipid_color_specific:
		view.set_color_cursor(lipid_color_specific)
	this_style = dictsum({'style':'Licorice %.3f 12.000000 12.000000'%licorice_thick,
		'goodsell':goodsell_lipids},lipid_style,{'color_specific':True} if lipid_color_specific else {})
	for hnum,hbond in enumerate(hydrogen_bonds):
		#---the first residue in the pair is the donor and we draw it here
		selstring = "(resname %s and resid %s)"%(hbond['resname_1'],hbond['resid_1'])
		if lipid_atoms is not None: 
			selstring = '(%s) and (%s)'%(selstring,' or '.join(['name %s'%i for i in lipid_atoms]))
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
		if lipid_atoms is not None: 
			selstring = '(%s) and (%s)'%(selstring,' or '.join(['name %s'%i for i in lipid_atoms]))
		selections.append('(%s and not hydrogen)'%selstring)
		select_args = dict(**{'partner_2_hbond_%d'%hnum:selections[-1]})
		select_args.update(**this_style)
		view.select(**select_args)
		view.command('set partner_2_hbond_%d_ref [atomselect top "%s"]'%(
			hnum,'(%s and name %s)'%(selstring,hbond['name_2'])))
	#---set the selections
	#---! debugging this now
	view.command('set me [atomselect top "%s"]'%' or '.join(['(%s)'%i for i in selections[:2]]))
	#---removed the VMD HBonds call because it is too permissive (allows too many atom types)
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
	#---trying to fix color
	#---previously used the special_ion_color when show_water but reverting to original colors for v15
	if lipid_colors or (show_water and False): view.set_color_cursor(special_ion_color)
	elif show_water: view.set_color_cursor({'MG':'red','Cal':'blue'}.get(cation,'purple'))
	else: view.set_color_cursor({'MG':'pink','Cal':'blue','NA':'green','K':'tan'}.get(cation,'purple'))
	#---show ions near both residues
	for resid_ind in [1,2]:
		view.select(**{'near_ions_%d'%resid_ind:
			"name %s and within 4.6 of (resname %s and resid %s)"%(
			cation,work.meta[sn]['ptdins_resname'],pair['resid_%d'%resid_ind]),
			'smooth':False,'style':'VDW 0.4 12.0','goodsell':goodsell_ions,'color_specific':True})
	#---option to show waters near the ions
	if show_water:
		if show_water.get('distance_color',False):
			for resid_ind in [1,2]:
				view.command('set pivots [atomselect top "[$near_ions_%d text]"]'%resid_ind)
				view.command(custom_water_coloring)
		#---deprecated method which passed the show water dictionary through to the view select
		elif show_water.get('deprecated_method',False):
			#---! this representation HAS MAJOR PROBLEMS. cg_bonds.tcl can be used to fix the extralong bonds
			#---! ...in PIP2 but it erases the bonds with water. we always just use the VDW representation
			#---! ...with large spheres. enormous amount of wasted time trying to fix this annoying bug.
			#---! ...need to get some kind of way of loading the TPR natively in VMD
			for resid_ind in [1,2]:
				view.select(**dictsum({'near_water_near_ions_%d'%resid_ind:
					"same residue as (water and within 3.0 of "+
						"(name %s and within 4.6 of (resname %s and resid %s)))"%(
					cation,work.meta[sn]['ptdins_resname'],pair['resid_%d'%resid_ind]),
					'smooth':False,'style':'licorice','goodsell':True},show_water))
		#---note that version 8 used normal water coloring but the other colors were not quite right
		#---...while version 14 colored waters by distance but probably went out too far
		#---deprecated method 
		elif show_water.get('hydration_shell',False):
			for resid_ind in [1,2]:
				view.select(**{'near_water_near_ions_%d'%resid_ind:
					"same residue as (water and within 3.05 of "+
						"(name %s and within 4.6 of (resname %s and resid %s)))"%(
					cation,work.meta[sn]['ptdins_resname'],pair['resid_%d'%resid_ind]),
					'smooth':False,'style':'licorice','goodsell':True})
		else: raise Exception('no valid show_water methods')
	#---emphasize the hydrogen bonds
	for hnum,hbond in enumerate(hydrogen_bonds):
		if hbond_style=='cylinder':
			if not show_water: view.command('draw color black')
			#---changed cylinder from silver to cyan for v15
			else: view.command('draw color iceblue')
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

def get_pairs(sn,resname_filter=None,nranked=10,any_combos=False,**kwargs):
	"""
	Find the most common pairs of a particular type.
	"""
	global data,memory
	bonds,obs,valid_frames = [data['hydrogen_bonding'][sn]['data'][i] 
		for i in ['bonds','observations','valid_frames']]
	#---discard atom names
	bonds_red = bonds[:,np.array([0,1,3,4])]
	#---this step induces a sort of sorts, so that actually the line where we define the ranking
	#---...ends up being in order over the mean of the sum of observations
	bonds_inds,bonds_counts = uniquify(bonds_red)
	#---subselect lipids and filter out intramolecular bonds
	lipid_resnames = np.array(work.vars['selectors']['resnames_lipid_chol'])
	#---further subselection by resname_filter
	if resname_filter: 
		lipid_resnames = [l for l in lipid_resnames if l in resname_filter]
		#---alias "PtdIns" back to the right residue name
		lipid_resnames = [work.meta[sn]['ptdins_resname'] if l=='PtdIns' 
			else l for l in resname_filter]
	#---two kinds of selections any combinations or specific pairings
	if len(lipid_resnames)==1 or any_combos:
		#---perform the selection
		subsel = np.where(np.all((
			np.in1d(bonds_red[bonds_inds][:,0],lipid_resnames),
			np.in1d(bonds_red[bonds_inds][:,2],lipid_resnames),
			bonds_red[bonds_inds][:,1]!=bonds_red[bonds_inds][:,3]
			),axis=0))[0]
	#---choosing two lipids implies that we want to do a pair otherwise we must set any_combos
	elif len(lipid_resnames)==2:
		subsel = np.where(
			np.any((
				np.all((
					bonds_red[bonds_inds][:,0]==lipid_resnames[0],
					bonds_red[bonds_inds][:,2]==lipid_resnames[1],
					bonds_red[bonds_inds][:,1]!=bonds_red[bonds_inds][:,3]
				),axis=0),
				np.all((
					bonds_red[bonds_inds][:,0]==lipid_resnames[1],
					bonds_red[bonds_inds][:,2]==lipid_resnames[0],
					bonds_red[bonds_inds][:,1]!=bonds_red[bonds_inds][:,3]
				),axis=0),
			),axis=0))[0]
	else: raise Exception('unclear major selection criterion')
	#---having filtered for desired bonds and already uniquified, we find the most common
	#---...after we map each unique reduced bond back to all of the actual bonds (with atoms)
	#---note that this step is slow
	status('slow where loop to map reduced bonds to observations for %s'%sn,tag='compute')
	map_red_to_obs = [np.where(np.all(bonds_red==i,axis=1))[0] 
		for i in bonds_red[bonds_inds][subsel]]
	#---load the mesh object once perform simulation to save memory
	if 'lipid_mesh_%s'%sn not in memory:
		del memory
		memory = {}
		status('fetching mesh object for %s'%sn,tag='load')
		memory['lipid_mesh_%s'%sn] = work.plotload('load_mesh',sns=[sn])[0][sn]['data']
		lipid_mesh = memory['lipid_mesh_%s'%sn]
	#---loop over rankings. these are somewhat arbitrary because we have subselected so hard
	ranking = np.argsort([obs.T[i].sum(axis=0).mean() for i in map_red_to_obs])[::-1]
	nranked = range(nranked) if type(nranked)==int else nranked
	pairspec = {}
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
		pairspec[rank_num] = dict([(i,bonds_by_pair[bonds_in_frame][inds[0]][ii]) 
			for ii,i in enumerate(bond_spec_keys.split())])
		#---keep the hydrogen bonds list for later
		pairspec[rank_num]['hydrogen_bonds'] = bond_spec['hydrogen_bonds']
		pairspec[rank_num]['bonds_by_pair'] = bonds_by_pair
		#---other useful information besides just the bond
		pairspec[rank_num]['frame'] = fr
		pairspec[rank_num]['frame_actual'] = frame_actual
	return pairspec

def pair_to_filetag(sn,pair,rank_num):
	"""Standard snapshot namer."""
	#---! needs more information
	filetag = 'fig.snapshot.%s.fr%d.%s_%s_o%d'%(sn,pair['frame_actual'],
		pair['resid_1'],pair['resid_2'],rank_num)
	return filetag

def get_snapshots_folder(supername,overwrite_snaps=False,**kwargs):
	"""Generate a snapshot folder from a method name."""
	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],'fig.hydrogen_bonding.%s'%supername)
	#---only make figures if the folder is empty. it only takes a few minutes -- no more coding logic
	if (os.path.isdir(tempdir) and not overwrite_snaps): 
		status('found %s and refusing to overwrite'%tempdir,tag='warning')
		return None
	try: os.mkdir(tempdir)
	except: status('directory exists so we are overwriting',tag='warning')
	status('snapshots dropping to %s'%tempdir,tag='note')
	return tempdir

def get_lipid_atoms(method,gro,resname):
	"""Atom names for a subset of a lipid."""
	struct = GMXStructure(gro)
	atom_names = struct.atom_names[np.where(np.all((
		struct.residue_indices==np.unique(
			struct.residue_indices[np.where(struct.residue_names==resname)])[0],
		struct.residue_names==resname),axis=0))]
	if method=='no_tail':
		#---include C31 and C21 otherwise oxygens dangle. otherwise exclude C21 and C31 onward
		#---! this method should be refined for residues other than PIP2
		return np.array([a for a in atom_names if not re.match(r'^(C[23][^1]|C[23]1\d)',a)])
	else: raise Exception('invalid method %s'%method)

def render_machine_mesh_review(render_spec,r2m,zoom=False):
	"""
	Render an image of the mesh for a particular pair of snapshot subjects. 
	This helps us understand some of the more unusual configurations.
	"""
	global memory
	#---intervention to draw meshes
	sn = render_spec['sn']
	this_frame = render_spec['frame']
	top_mono = work.meta[sn].get('index_top_monolayer',0)
	highlight_s,lowlight_s = 60,30
	ax = plt.subplot(111)
	resnames_all = data['lipid_abstractor'][sn]['data']['resnames'][
		np.where(data['lipid_abstractor'][sn]['data']['monolayer_indices']==top_mono)[0]]
	#---! the following load is repetitive with some code elsewhere
	#---load the mesh object once perform simulation to save memory
	if 'lipid_mesh_%s'%sn not in memory:
		del memory
		memory = {}
		status('fetching mesh object for %s'%sn,tag='load')
		memory['lipid_mesh_%s'%sn] = work.plotload('load_mesh',sns=[sn])[0][sn]['data']
	lipid_mesh = memory['lipid_mesh_%s'%sn]
	#---we ignore all other lipids if we want to make the zoomed version
	if not zoom: ax.scatter(
		lipid_mesh['%d.%d.points'%(top_mono,this_frame)][:,0],
		lipid_mesh['%d.%d.points'%(top_mono,this_frame)][:,1],s=lowlight_s,alpha=0.4,lw=0,
		color=[colorize(work.meta[sn],resname=r) for r in resnames_all],zorder=2,clip_on=False)
	#---! removed: ax.scatter(*mesh['0.%d.points'%fr][:,:2][ghosts].T,s=4,c='b')
	ax.set_xlim((0,lipid_mesh['vecs'][this_frame][0]))
	ax.set_ylim((0,lipid_mesh['vecs'][this_frame][1]))
	ax.set_aspect('equal')
	ax.set_xticks([])
	ax.set_yticks([])
	ax.tick_params(axis=u'both',which=u'both',length=0)
	#---render the hydrogen bonds
	for hbond in render_spec.get('hydrogen_bonds',[]):
		indices = np.array([-1,-1])
		for pnum,partner in enumerate([1,2]):
			resid = int(hbond['resid_%d'%partner])
			resname = hbond['resname_%d'%partner]
			index = r2m[resid]
			indices[pnum] = index
			ax.scatter(*lipid_mesh['%d.%d.points'%(top_mono,this_frame)][:,:2][index].T,s=highlight_s,
				c=colorize(work.meta[sn],resname=resname),zorder=5,edgecolor='k',lw=2,clip_on=False)
		ax.plot(*lipid_mesh['%d.%d.points'%(top_mono,this_frame)][:,:2][indices].T,zorder=4,lw=1,c='k')
	for assoc in render_spec.get('associates',[]):
		resid = int(assoc['resid'])
		resname = assoc['resname']
		index = r2m[resid]
		ax.scatter(*lipid_mesh['%d.%d.points'%(top_mono,this_frame)][:,:2][index].T,s=highlight_s,
			c=colorize(work.meta[sn],resname=resname),zorder=3,edgecolor='k',lw=2,clip_on=False)
	ax.axis('off')
	#---save the mesh
	tag_mesh_fn = re.sub('snapshot','snapshot_mesh_review%s'%('_zoom' if zoom else ''),render_spec['fn'])
	picturesave(tag_mesh_fn,render_spec['tempdir'],backup=False,version=True,meta={},extras=[])

def render_machine(**method):
	"""
	Handle the request to render a snapshot by figuring out which part of the simulation to plot, how 
	it should look, and where to plot it.
	"""
	global memory
	sn = method['sn']
	overwrite_snaps = method.get('overwrite_snaps',True)
	nranked = method.get('nranked',1)
	nranked = range(nranked) if type(nranked)==int else nranked
	#---prepare filename and directory
	folder = get_snapshots_folder(overwrite_snaps=overwrite_snaps,**method)
	if not folder:
		status('folder exists so we are skipping this method! (see overwrite_snaps flag)',tag='warning')
		return
	#---collect most-populous/relevant pairs
	pairspec = get_pairs(**method)
	for rank_num in nranked:
		status('rendering snapshot %d from %s'%(rank_num,nranked),tag='render')
		pair = pairspec[rank_num]
		#---prepare the render
		filetag = pair_to_filetag(sn=method['sn'],pair=pair,rank_num=rank_num)
		render_spec = dict(sn=method['sn'],pair=pair,
			frame=pair['frame_actual'],tempdir=folder,fn=filetag,
			hydrogen_bonds=pair['hydrogen_bonds'])
		render_spec.update(**work.get_gmx_sources(sn=method['sn'],calc=calc))
		#---apply modifiers
		render_modifiers = method.get('render_modifiers',{})
		if type(render_modifiers) in [list,tuple]: 
			render_modifiers = dict([(k,True) for k in render_modifiers])
		associates,barnacles,lipid_color_specific = None,None,None
		for mod in render_modifiers:			
			if mod=='no_tail':
				render_spec['lipid_atoms'] = get_lipid_atoms('no_tail',
					gro=render_spec['gro'],resname=work.meta[method['sn']]['ptdins_resname'])
			elif mod=='associates': associates = render_modifiers['associates']
			elif mod=='barnacles': barnacles = render_modifiers['barnacles']
			elif mod=='lipid_color_specific': 
				render_spec['lipid_color_specific'] = render_modifiers['lipid_color_specific']
			elif mod=='lipid_colors' and render_modifiers['lipid_colors']:
				render_spec['lipid_colors'] = render_modifiers['lipid_colors']
			else: raise Exception('invalid modifier %s'%mod)
		#---pass along important flags
		for key in ['lipid_style','show_water']:
			if key in method: render_spec[key] = method[key]
		#---check the lipid distance to see if we are broken across PBCs
		resids = [np.where(data['lipid_abstractor'][sn]['data']['resids']==int(
			pair['hydrogen_bonds'][0]['resid_%d'%i]))[0][0] for i in [1,2]]
		com_distance = data['lipid_abstractor'][sn]['data']['points'][
			pair['frame_actual']][resids].ptp(axis=0)
		if np.any(com_distance>=data['lipid_abstractor'][sn]['data']['vecs'][pair['frame_actual']]/2.):
			status('lipids are broken over PBCs and this rendering algo cannot work so skipping',
				tag='warning')
			continue
		#---identify the barnacles: lipids nearby that we care about and which also appear in the 
		#---...listing of hydrogen bonds
		if barnacles:
			this_bond = pair['bonds_by_pair'][0]
			resname_1,resname_2,resid_1,resid_2 = [this_bond[bond_spec_keys.split().index(i)] 
				for i in ['resname_1','resname_2','resid_1','resid_2']]
			#---recapitulate some of the original parsing from get_pairs
			bonds,obs = [data['hydrogen_bonding'][sn]['data'][i] for i in ['bonds','observations']]
			#---discard atom names
			bonds_red = bonds[:,np.array([0,1,3,4])]
			#---this step induces a sort of sorts, so that actually the line where we define the ranking
			#---...ends up being in order over the mean of the sum of observations
			bonds_inds,bonds_counts = uniquify(bonds_red)
			subsel_barnacles = np.where(
				np.any((
					np.all((
						np.in1d(bonds_red[bonds_inds][:,0],[resname_1,resname_2]),
						np.in1d(bonds_red[bonds_inds][:,1],[resid_1,resid_2]),
						np.in1d(bonds_red[bonds_inds][:,2],barnacles),),axis=0),
					np.all((np.in1d(bonds_red[bonds_inds][:,2],[resname_1,resname_2]),
						np.in1d(bonds_red[bonds_inds][:,3],[resid_1,resid_2]),
						np.in1d(bonds_red[bonds_inds][:,0],barnacles),),axis=0),
					),axis=0)
				)[0]
			#---only render the bonds if the barnacle has a bond in this frame
			for subind in subsel_barnacles:
				if obs.T[bonds_inds][subind][pair['frame']]==1.0:
					bond = bonds[bonds_inds][subind]
					render_spec['hydrogen_bonds'].append(dict([(i,bond[ii]) 
						for ii,i in enumerate(bond_spec_keys.split())]))
		#---INDEXING SEQUENCE
		#---construct several mappings between relevant indices
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
		def get_mesh_neighbors(index):
			"""Get the neighbors for a vertex."""
			neighbors = np.unique(simplices[np.where(np.any(simplices==index,axis=1))[0]])
			ghosts = lipid_mesh['%d.%d.ghost_ids'%(top_mono,pair['frame'])][neighbors]
			return ghosts
		#---check for associates, lipids nearby on the mesh
		if associates:
			#---! the following load is repetitive with some code elsewhere
			#---load the mesh object once perform simulation to save memory
			if 'lipid_mesh_%s'%sn not in memory:
				del memory
				memory = {}
				status('fetching mesh object for %s'%sn,tag='load')
				memory['lipid_mesh_%s'%sn] = work.plotload('load_mesh',sns=[sn])[0][sn]['data']
			lipid_mesh = memory['lipid_mesh_%s'%sn]
			simplices = lipid_mesh['%d.%d.simplices'%(top_mono,pair['frame'])]
			resname_resids = np.transpose((
				data['lipid_abstractor'][sn]['data']['resnames'],
				data['lipid_abstractor'][sn]['data']['resids']))
			#---loop over which of the partners in the pair
			for wp in [1,4]:
				#---get resid for this member of the pair
				ind_mesh  = r2m[int(pair['bonds_by_pair'][0][wp])]
				neighbors_mesh = get_mesh_neighbors(ind_mesh)
				neighbors_resids = np.array([m2r[i] for i in get_mesh_neighbors(ind_mesh)])
				neighbors_resnames = [data['lipid_abstractor'][sn]['data']['resnames'][r2i[i]] 
					for i in neighbors_resids]
				#---reverse mapping to monolayers
				near_neighbors = np.where(np.in1d(neighbors_resnames,associates))[0]
				if len(near_neighbors)>0 and 'associates' not in render_spec:
					render_spec['associates'] = []
				#---include all nearest neighbors in the visualization
				for resid in neighbors_resids[near_neighbors]:
					render_spec['associates'].append(dict(
						resname=resnames[r2i[resid]],
						resid=str(resids[r2i[resid]])))
		#---render the mesh for reference
		render_machine_mesh_review(render_spec=render_spec,r2m=r2m)
		render_machine_mesh_review(render_spec=render_spec,r2m=r2m,zoom=True)
		#---if associates are broken over PBCs we render the mesh for reference but not the snapshot
		resids_this = np.array([r2i[k] for k in np.concatenate((
			[render_spec['pair']['resid_%d'%j] for j in [1,2]],
			[i['resid_%d'%j] for i in render_spec['hydrogen_bonds'] for j in [1,2]],
			[i['resid'] for i in render_spec.get('associates',[])],
			)).astype(int)])
		com_distance = data['lipid_abstractor'][method['sn']]['data']['points'][pair['frame_actual']][
			resids_this].ptp(axis=0)
		if np.any(com_distance>=data['lipid_abstractor'][method['sn']]['data']['vecs'][
			pair['frame_actual']]/2.):
			status('associates are broken over PBCs so we rendered the mesh but not the snapshot',
				tag='warning')
			continue
		#---render the snapshot
		render_lipid_pair(**render_spec)

###---MAIN

if 'arrange_snapshots' in routine:
	"""
	Arrange some snapshots into a panel plot.
	"""
	import matplotlib.image as mpimg

	def get_blank_border(img):
		"""Return the border limits for removing whitespace."""
		nonblank = np.any(img!=1.0,axis=2)*1
		lims = [map(lambda x:(x[0]-1,x[-1]+1),
			[np.where(np.any(nonblank!=0,axis=j))[0]])[0] for j in range(2)]
		return lims

	def sidestack_images(fns):
		"""Combine snapshots into one image with same zoom and minimal whitespace."""
		#---preassemble images for a particular row
		imgs = [mpimg.imread(fn) for fn in fns]
		borders = np.array([get_blank_border(img) for img in imgs])
		lims = np.array([[i.min() for i in borders.T[0]]]+[[i.max() for i in borders.T[1]]]).T
		#---stack horizontally
		buffer_width = 100
		imgs_zoomed = [img[slice(*lims[1]),slice(*lims[0])] for img in imgs]
		imgs_cropped = [i[:,slice(*(get_blank_border(i)[0]))] for i in imgs_zoomed]
		buffer_strip = np.ones((imgs_cropped[0].shape[0],buffer_width,3))
		#---add a vertical buffer strip
		imgs_combo = np.concatenate([np.concatenate((i,buffer_strip),axis=1) 
			for i in imgs_cropped[:-1]]+[imgs_cropped[-1]],axis=1)
		return imgs_combo

	layouts = {
		'panel_n4':{
			#---previously panel_n1 rendered with fig.hydrogen_bonding.v8_ptdins_ptdins
			#---on panel_n2 used full lipids from v10, not solvated
			#---panel_n3 uses the updated figures with neat color scheme
			#---panel_n4 unified the ion colors (green) and made some minor changes
			'figsize':(12,12),
			'style':'square_tiles','folder':'fig.hydrogen_bonding.v12_ptdins_solvated','files':[
				'fig.snapshot.membrane-v531.fr703.795_796_o0.png',
				'fig.snapshot.membrane-v531.fr50.780_782_o6.png',
				'fig.snapshot.membrane-v532.fr11.779_781_o0.png',
				'fig.snapshot.membrane-v532.fr1.793_790_o1.png',
				'fig.snapshot.membrane-v533.fr13.798_794_o0.png',
				'fig.snapshot.membrane-v534.fr344.766_765_o8.png',
				]},
		'panel_n5':{
			#---previously panel_n1 rendered with fig.hydrogen_bonding.v8_ptdins_ptdins
			#---on panel_n2 used full lipids from v10, not solvated
			#---panel_n3 uses the updated figures with neat color scheme
			#---panel_n4 unified the ion colors (green) and made some minor changes
			'figsize':(12,10),
			'style':'square_tiles','folder':'fig.hydrogen_bonding.v12_ptdins_solvated','files':[
				'fig.snapshot.membrane-v531.fr703.795_796_o0.png',
				'fig.snapshot.membrane-v531.fr50.780_782_o6.png',
				'fig.snapshot.membrane-v532.fr11.779_781_o0.png',
				'fig.snapshot.membrane-v532.fr1.793_790_o1.png',
				]},
		'panel_n5_with_rdf':{
			'figsize':(6,12),'no_letters':True,
			'style':'square_tiles','folder':'.','files':[
				'fig.hydration_distribution.v1.png',
				'fig.snapshots.panel_n5.v1.png',
				]},}

	#---loop over layouts
	extras = []
	tagbox = dict(facecolor='w',lw=1,alpha=1.0,boxstyle="round,pad=0.5")
	for tag,layout in layouts.items():
		fns = layout['files']
		folder = layout['folder']
		if layout['style']=='square_tiles':
			ntiles = len(fns)
			axes,fig = square_tiles(ntiles=ntiles,figsize=layout['figsize'],wspace=0.1,hspace=0.1)
			for fnum,fn in enumerate(fns):
				#---second axis holds the exemplars
				ax = axes[fnum]
				#---use the side-stacker to remove whitespace, even for only one image
				image = sidestack_images([os.path.join(work.plotdir,folder,fn)])
				ax.imshow(image)
				ax.axis('off')
				if not layout.get('no_letters',False): 
					tb = ax.text(0.0,1.0,chr(ord('A')+fnum),fontsize=14,
						bbox=tagbox,rotation=0,ha="center",va="top",color='k',
						transform=ax.transAxes)
				extras.append(tb)
			picturesave('fig.snapshots.%s'%tag,work.plotdir,
				backup=False,version=True,meta={},extras=extras)
		else: raise Exception('invald layout style %s'%layout['style'])

if 'systematic_snapshotter' in routine:

	"""
	The above/previous hydrogen bonding plots (removed but available in git and the graveyard) have several 
	extremely useful features.
		1. identify relevant pairs.
		2. plot the pairs of lipids
		3. plot nearby ions, plot hydrogen bonds explicitly after orienting the lipids to show the bonds, etc
		4. plot mesh images so we can see a birdseye view of the configuration in case it is confusing
	We wish to part out these features into functions so we can make similar plots that show the ion-water 
	portion of this calculation.
	A single "experiment" is the union of:
		1. how to select pairs
		2. which bonds to draw
		3. where to save the file
	Note that we wish eventually to have this function replace the hydrogen bond rendering function but
	also render versions of the pairs with e.g. solvated water and other interesting aesthetic features.
	Deprecated image sets:
		1. everything before v7 is all development-only
		2. v7, v8 original set of images before discussions with david (v8 was overwritten with tails)
		3. v10 has useful images, but it was superceded
		4. v13 was used to develop the colored lipids
	Active image sets:
		1. v14 has hydrogen bonds formerly from v7 and v10 with lipids in a single color
		2. v12 has water colored by distance to the nearest ions
	NOTE that v536 still has incorrect meshes
	"""

	#---! this should match art_ptdins.py if plots are merged
	special_ion_color = 'green'
	lipid_colors = {'ptdins':'purple','PI2P':'purple','SAPI':'purple',
		'P35P':'purple','CHL1':'gray','DOPS':'red','DOPE':'blue2'}

	###---EXPERIMENTS
	methods = []
	#---overwriting is on by default so use "if False:" to avoid redoing things
	#---version 8: plot just the head groups for ptdins-ptdins pairs 
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v8_ptdins_solvated',
		render_modifiers=['no_tail'],
		lipid_style={'goodsell':False,'style':'Licorice 0.15 12.0 12.0','diffuse':True},
		show_water={'style':'licorice'},
		resname_filter=['PtdIns','PtdIns'])
		for sn in sns])
	#---version 14 (recapping v7,v10 with color changes) shows hydrogen bonds between different lipid types
	#---! this is completed but there are still problems with v536
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v14_ptdins_with_chl1',
		resname_filter=['PtdIns'],
		render_modifiers={'associates':['CHL1'],'barnacles':['CHL1'],'lipid_colors':lipid_colors})
		for sn in sns])
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v14_ptdins_dope_with_chl1',
		resname_filter=['PtdIns','DOPE'],
		render_modifiers={'associates':['CHL1'],'barnacles':['CHL1'],'lipid_colors':lipid_colors})
		for sn in sns])
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v14_ptdins_chl1_with_dope_chl1',
		resname_filter=['PtdIns','CHL1'],
		render_modifiers={'associates':['CHL1','DOPE'],'barnacles':['CHL1','DOPE'],
			'lipid_colors':lipid_colors})
		for sn in ['membrane-v532','membrane-v534']])
	#---! incomplete version 11: hydrogen bond plots extended to new lipids
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v11_ptdins_dops_assoc_chl1_barnacles_chl1_ptdins',
		resname_filter=['PtdIns','DOPS'],
		render_modifiers={'associates':['CHL1'],'barnacles':['CHL1','PtdIns']})
		for sn in sns])
	#---version 12: color-coded water distances
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v12_ptdins_solvated',
		render_modifiers={'no_tail':True,'lipid_color_specific':'black'},
		lipid_style={'goodsell':True,'style':'Licorice 0.15 12.0 12.0'},
		show_water={'style':'licorice'},
		resname_filter=['PtdIns'])
		for sn in sns])
	#---version 15: closer water without coloring
	if 1: methods.extend([dict(sn=sn,nranked=10,
		supername='v15_ptdins_solvated',
		render_modifiers={'no_tail':True,'lipid_color_specific':'black'},
		lipid_style={'goodsell':True,'style':'Licorice 0.15 12.0 12.0'},
		show_water={'style':'licorice','hydration_shell':True},
		resname_filter=['PtdIns'])
		for sn in sns])
	#---versions 10: recap the v7 plots from previous version of the code
	if 0: methods.extend([dict(sn=sn,nranked=10,
		supername='v13_ptdins_with_chl1_colors',
		resname_filter=['PtdIns'],
		render_modifiers={'associates':['CHL1'],'barnacles':['CHL1'],'lipid_colors':lipid_colors})
		for sn in sns])

	#---loop over methods, which also loops over simulations
	#---to debug just set `methods = methods[:1]` or pick one manually
	for mnum,method in enumerate(methods): 
		status('rendering method set',i=mnum,looplen=len(methods),tag='status')
		render_machine(**method)

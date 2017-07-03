#!/usr/bin/python

#---plot prep
if 'plotload' not in globals(): execfile('/etc/pythonstart')
execfile('./omni/base/header.py')
from plotter import *
from base.store import plotload
execfile('./calcs/specs/figures.py')
execfile('./calcs/specs/colors.py')
import numpy as np

#---settings
plotname = 'hydrogen_bonding'
withchar = ','
routine = ['detailed_plot','simple_plot','hbond_snaps','simple_plot_breakdown','pip2_only_breakdown'][-1:]
ns = next_script = sys.argv[0]
print_version = 'v20160923'

#---LOAD
#-------------------------------------------------------------------------------------------------------------

#---load everything
if 'data' not in globals(): data,calcs = plotload(plotname,work)

postdat = []
for sn in work.sns():
	dat = data[sn]['data']
	bl = dat['bond_list']
	nframes = len(dat['variations'])
	#---repetitive selections because I am too lazy to loop and we will probably do more
	#---both = 'not all SOL'
	subsel_both = np.where(~np.all((bl[:,0]=='SOL',bl[:,2]=='SOL'),axis=0))[0]
	counts_both = dat['variations'][:,subsel_both].sum(axis=1)
	#---lipid = 'not any SOL'
	subsel_lipid = np.where(~np.any((bl[:,0]=='SOL',bl[:,2]=='SOL'),axis=0))[0]
	counts_lipid = dat['variations'][:,subsel_lipid].sum(axis=1)
	#---lipid_water = 'exactly one SOL'
	subsel_lipid_water = np.where(np.sum((bl[:,0]=='SOL',bl[:,2]=='SOL'),axis=0)==1)[0]
	counts_lipid_water = dat['variations'][:,subsel_lipid_water].sum(axis=1)
	#---repeating some of the above for bonds involving PIP2
	#---any_pip2 = 'either or both PIP2'
	any_pip2 = np.any((bl[:,0]==work.meta[sn]['ptdins_resname'],
		bl[:,2]==work.meta[sn]['ptdins_resname']),axis=0)
	#---both_pip2 = '(not all SOL) and (either or both PIP2)' ... (both SOL or lipid) with PIP2
	subsel_both_pip2 = np.where(
		np.all((
			~np.all((bl[:,0]=='SOL',bl[:,2]=='SOL'),axis=0),
			any_pip2),axis=0))[0]
	counts_both_pip2 = dat['variations'][:,subsel_both_pip2].sum(axis=1)
	subsel_lipid_pip2 = np.where(np.all((~np.any((bl[:,0]=='SOL',
		bl[:,2]=='SOL'),axis=0),any_pip2),axis=0))[0]
	counts_lipid_pip2 = dat['variations'][:,subsel_lipid_pip2].sum(axis=1)
	subsel_lipid_water_pip2 = np.where(np.all((np.sum((bl[:,0]=='SOL',
		bl[:,2]=='SOL'),axis=0)==1,any_pip2),axis=0))[0]
	counts_lipid_water_pip2 = dat['variations'][:,subsel_lipid_water_pip2].sum(axis=1)
	#---repeating once more to isolate bonds involving PIP2 and cholesterol
	#any_pip2_chl = np.all((np.any((bl[:,0]==work.meta[sn]['ptdins_resname'],
	#	bl[:,2]==work.meta[sn]['ptdins_resname']),axis=0),
	#	np.any((bl[:,0]==work.meta[sn]['ptdins_resname'],
	#	bl[:,2]=='CHL1'),axis=0),
	#	),axis=0)
	#---! major screw-up above
	#---! another major screw-up below (this makes three, 
	#---! ...including the one on the atom-specfic items in the computation!)
	#any_pip2_chl = np.all((np.any((bl[:,0]==work.meta[sn]['ptdins_resname'],
	#	bl[:,2]=='CHL1'),axis=0),
	#	np.any((bl[:,2]==work.meta[sn]['ptdins_resname'],
	#	bl[:,0]=='CHL1'),axis=0),
	#	),axis=0)
	any_pip2_chl = np.any((np.all((bl[:,0]==work.meta[sn]['ptdins_resname'],
		bl[:,2]=='CHL1'),axis=0),
		np.all((bl[:,2]==work.meta[sn]['ptdins_resname'],
		bl[:,0]=='CHL1'),axis=0),
		),axis=0)
	counts_pip2_chl = dat['variations'][:,np.where(any_pip2_chl)].sum(axis=1).sum(axis=1)

	any_chl_lipid = np.any((np.all((bl[:,0]!='SOL',
		bl[:,2]=='CHL1'),axis=0),
		np.all((bl[:,2]!='SOL',
		bl[:,0]=='CHL1'),axis=0),
		),axis=0)
	counts_chl_lipid = dat['variations'][:,np.where(any_chl_lipid)].sum(axis=1).sum(axis=1)

	subsel_pip2_pip2 = np.where(np.all((bl[:,0]==work.meta[sn]['ptdins_resname'],
		bl[:,2]==work.meta[sn]['ptdins_resname']),axis=0))[0]
	counts_pip2_pip2 = dat['variations'][:,subsel_pip2_pip2].sum(axis=1)

	postdat.append((sn,{
		('lipid_water','mean'):counts_lipid_water.mean(),
		('lipid_water','std'):counts_lipid_water.std(),
		('lipid','mean'):counts_lipid.mean(),
		('lipid','std'):counts_lipid.std(),
		('both','mean'):counts_both.mean(),
		('both','std'):counts_both.std(),
		('lipid_water_pip2','mean'):counts_lipid_water_pip2.mean(),
		('lipid_water_pip2','std'):counts_lipid_water_pip2.std(),
		('lipid_pip2','mean'):counts_lipid_pip2.mean(),
		('lipid_pip2','std'):counts_lipid_pip2.std(),
		('both_pip2','mean'):counts_both_pip2.mean(),
		('both_pip2','std'):counts_both_pip2.std(),
		('pip2_chl1','mean'):counts_pip2_chl.mean(),
		('pip2_chl1','std'):counts_pip2_chl.std(),
		('pip2_pip2','mean'):counts_pip2_pip2.mean(),
		('pip2_pip2','std'):counts_pip2_pip2.std(),
		('lipid_chl1','mean'):counts_chl_lipid.mean(),
		('lipid_chl1','std'):counts_chl_lipid.std(),
		}))

n_most_popped_bonds = 3

#---precache identities in case they are already rendered
sns_group = work.sns()
if 'snap_plan' not in globals():
	snap_plan,bonds_tabbed = {},{}
	#---loop over simulations
	for sn_this in sns_group:
		#---print the top N most populated hydrogen bonds			

		result = data[sn_this]['data']
		#---saved the PIP2-lipid bonds specifically
		#---note that this is a single list over all frames, unpackable by pip2_non_sol_hbonds_bpf
		hbonds = result['pip2_non_sol_hbonds']
		resnum_index = 5

		if len(hbonds)==0: 
			print('[WARNING] no hbonds in %s ... skipping'%sn_this)
			bonds_tabbed[sn_this] = {'bonds':np.array([[]]),'counts':np.array([])}
			continue

		#---! intervene here to redefine hbonds by name, resid
		#hbonds = hbonds[:,np.array([0,1,4,5])]
		#resnum_index = 3

		#---unique-ify the hbonds and find the most prominent ones
		#---since it's unique, we don't care about frames -- there can be no multiples
		try: hbonds_coded = np.ascontiguousarray(hbonds).view(
			np.dtype((np.void,hbonds.dtype.itemsize*hbonds.shape[1])))
		#---some simulations have no bonds
		except: continue
		hbonds_unique_coded,counts = np.unique(hbonds_coded,return_counts=True)
		#---reconstitute the unique strings as tuples
		_,idx,counts_not_sorted = np.unique(hbonds_coded,return_index=True,return_counts=True)
		hbonds_unique_unsort = hbonds[idx]
		#---select the first most-populated hbond with two lipids
		hbonds_unique = hbonds_unique_unsort[np.argsort(counts_not_sorted)[::-1]]
		counts = counts_not_sorted[np.argsort(counts_not_sorted)[::-1]]

		bonds_tabbed[sn_this] = {'bonds':hbonds_unique,'counts':counts}

		for rank in range(n_most_popped_bonds):
			#---! ryan and david manually checked that the two above variables are indexed together
			index_of_type = [ii for ii,i in enumerate(hbonds_unique) if i[1]!=i[resnum_index]][rank]
			#---check the counts
			print('[NOTE] selecting bond type %s'%str(hbonds_unique[index_of_type]))
			print('[NOTE] this type is observed in %d frames'%counts[index_of_type])
			#---get the first occurence of this bond type in the trajectory
			index = [ii for ii,i in enumerate(hbonds) if 
				np.all(i==hbonds_unique[index_of_type])][rank]
			#---get the frame number from the list of bonds per frame
			frame = np.where(index<np.cumsum(result['pip2_non_sol_hbonds_bpf']))[0][0]
			#---get resid and resnames
			if not np.all(hbonds_unique[index_of_type]==hbonds[index]): raise Exception('problem')
			#rn1,ri1,an1,hname2,rn2,ri2,an2 = hbonds[index]
			snap_plan[(sn_this,rank)] = list(hbonds[index])+[frame]+[index]

def bonds_tabulate(sn,*args):
	"""
	Retrieve counts and lists of certain bonds.
	"""
	#---specify the format of the bond list
	items = np.array('resname resid atomname hydrogen resname resnum atomname'.split())
	bonds = bonds_tabbed[sn]['bonds']
	counts = bonds_tabbed[sn]['counts']
	assert len(args) in [1,2]
	rules = []
	#---for each arg we expand the selection criteria
	for arg in args:
		#---within each rule, we use "OR"
		for restrict in ['resname','atomname']:
			if restrict in arg:
				try: 
					rule = np.any([bonds[:,ii]==arg[restrict] for ii in np.where(items==restrict)[0]],axis=0)
					rules.append(rule)
				except: pass
	#---process all rules via "AND"
	if not rules: return {'inds':np.array([]),'bonds':np.array([]),'counts':np.array([]),'mean':0}
	inds = np.all(rules,axis=0)
	nframes = float(len(data[sn]['data']['variations']))
	return {'inds':inds,'bonds':bonds[inds],'counts':counts[inds],'mean':counts[inds].sum()/nframes}

sn = 'membrane-v532'
#---check that the counts on the bar graphs make sense
tab_pip2_chl1 = bonds_tabulate(sn,{'resname':'PI2P'},{'resname':'CHL1'})
tab_pip2 = bonds_tabulate(sn,{'resname':'PI2P'})
print((tab_pip2['mean'],tab_pip2_chl1['mean']))
#---continuing to debug
"""
STILL A DISCREPANCY
>>> dict(postdat)[sn][('pip2_chl1','mean')]
44.75530586766542
>>> dict(postdat)[sn][('lipid_pip2','mean')]
141.40074906367042
>>> print((tab_pip2['mean'],tab_pip2_chl1['mean']))
(141.40074906367042, 4.1285892634207242)
"""

#---check the counts for self bonds
sn = 'membrane-v532'
nframes = float(len(data[sn]['data']['variations']))
tabs = bonds_tabulate(sn,{'resname':work.meta[sn]['ptdins_resname']},{'resname':work.meta[sn]['ptdins_resname']})
print(tabs['counts'][np.where(np.all((tabs['bonds'][:,1]!=tabs['bonds'][:,5],tabs['bonds'][:,0]==work.meta[sn]['ptdins_resname'],tabs['bonds'][:,4]==work.meta[sn]['ptdins_resname']),axis=0))].sum()/nframes)
print(tabs['counts'][np.where(np.all((tabs['bonds'][:,1]==tabs['bonds'][:,5],tabs['bonds'][:,0]==work.meta[sn]['ptdins_resname'],tabs['bonds'][:,4]==work.meta[sn]['ptdins_resname']),axis=0))].sum()/nframes)
dict(postdat)[sn][('pip2_pip2','mean')]
#---so most of the PIP2-PIP2 bonds are within the lipid
#---this also confirms the basic/advanced tabulation matches
pct_pip2_pip2_hbonds = {}
for sn in sns_group:
	if sn in bonds_tabbed:
		tabs = bonds_tabulate(sn,
			{'resname':work.meta[sn]['ptdins_resname']},
			{'resname':work.meta[sn]['ptdins_resname']})
		if len(tabs['bonds'])!=0: 
			nframes = float(len(data[sn]['data']['variations']))
			self_other = tabs['counts'][np.where(np.all((tabs['bonds'][:,1]!=tabs['bonds'][:,5],
				tabs['bonds'][:,0]==work.meta[sn]['ptdins_resname'],
				tabs['bonds'][:,4]==work.meta[sn]['ptdins_resname']),axis=0))].sum()/nframes
			self_self = tabs['counts'][np.where(np.all((tabs['bonds'][:,1]==tabs['bonds'][:,5],
				tabs['bonds'][:,0]==work.meta[sn]['ptdins_resname'],
				tabs['bonds'][:,4]==work.meta[sn]['ptdins_resname']),axis=0))].sum()/nframes
			pct_pip2_pip2_hbonds[sn] = self_other/(self_self+self_other)
print(pct_pip2_pip2_hbonds)
print(max(pct_pip2_pip2_hbonds.values()))

###---PLOT ROUTINES

if 'detailed_plot' in routine:

	width = 1.2
	fsbase = 18
	patches = []
	#---re-order by cation
	cations_order = ['NA','K','MG','Cal']
	sns = [sn for cat in cations_order for sn in work.collection('all') 
		if work.meta[sn]['composition_name'] == 'asymmetric' and work.meta[sn]['cation']==cat]
	#---actually reorder by mean
	sns = [dict(postdat).keys()[i] for i in np.argsort([dict(postdat)[sn][('lipid','mean')] 
		for sn in dict(postdat)])[::-1]]
	allres = work.vars['selectors']['resnames_lipid_chol']
	resnames_pip2 = work.vars['selectors']['resnames_PIP2']
	tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.2")
	tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.2")
	axes,fig = panelplot(figsize=(30,15),layout={'out':{'grid':[2,1]},
		'ins':[{'grid':[1,3]} for j in range(2)]})
	for rnum,blist in enumerate([
		['both','lipid_water','lipid'],
		['both_pip2','lipid_water_pip2','lipid_pip2']]):
		sns = [dict(postdat).keys()[i] for i in 
			np.argsort([dict(postdat)[sn][('lipid' if rnum==0 else 'lipid_pip2','mean')] 
			for sn in dict(postdat)])[::-1]]
		for bnum,btype in enumerate(blist):
			ax = axes[rnum][bnum]
			ymin,ymax = 0,0
			counter,centers,scs = 0,[],[]
			for ss,sn in enumerate(sns):
				color = colorize(work.meta[sn],comparison='protonation')
				dc = dict(postdat)[sn][(btype,'mean')]
				dc_std = dict(postdat)[sn][(btype,'std')]
				ymax = dc+dc_std if dc+dc_std > ymax else ymax
				ymin = dc+dc_std if dc+dc_std < ymin else ymin
				if work.meta[sn]['composition_name']=='symmetric': alpha = 0.25
				else: alpha = 0.65
				bar = ax.bar(ss*width,dc,width=width,alpha=alpha,zorder=2,lw=0,color=color,edgecolor=color)
				#---also plot standard deviations
				std_lines = [([(ss+0.25)*width,(ss+0.75)*width],[dc-dc_std,dc-dc_std]),
					([(ss+0.25)*width,(ss+0.75)*width],[dc+dc_std,dc+dc_std]),
					([(ss+0.5)*width,(ss+0.5)*width],[dc-dc_std,dc+dc_std])]
				for l in std_lines: 
					ax.plot(l[0],l[1],color=color,lw=4,solid_capstyle='round',alpha=1,zorder=3)
					ymax = l[1] if dc+dc_std > ymax else ymax
				centers.append(np.mean(np.arange(counter,counter+width)+1))
				counter += width+0
			ax.set_xticks(centers)
			outer_labels = []
			ax.set_xticklabels([])
			for snum,sn in enumerate(sns):
				text = work.meta[sn]['ion_label']
				tb = ax.text(width*(snum+0.5),(ymax-ymin)*-0.02,text,
					bbox=tagbox_ion,ha='center',va='top',rotation=90 if len(text)>25 else 0,
					color='k',fontsize=fsbase-8)
				outer_labels.append(tb)
				text = work.meta[sn]['ptdins_label']
				tb = ax.text(width*(snum+0.5),(ymax-ymin)*0.02,text,
					bbox=tagbox_ptdins,rotation=90,ha="center",va="bottom",color='k',fontsize=fsbase-2)
			ax.set_ylim(ymin,ymax*1.1)
			ax.set_yticklabels(['%.f'%i for i in ax.get_yticks()],fontsize=fsbase+2)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_xlim(-1,counter+1)
			ax.set_title({
				'both':'lipid hydrogen bonds','lipid_water':'lipid-water hydrogen bonds',
				'lipid':'lipid-lipid hydrogen bonds',
				'both_pip2':'lipid hydrogen bonds (with PIP2)',
				'lipid_water_pip2':'lipid-water hydrogen bonds (with PIP2)',
				'lipid_pip2':'lipid-lipid hydrogen bonds (with PIP2)',
				}[btype],fontsize=fsbase+4)
			ax.set_ylabel('observed number of bonds',fontsize=fsbase+2)
			patches.extend([(mpl.lines.Line2D([],[],marker='o',markersize=15,markeredgecolor='w',mew=2,lw=0,
				markerfacecolor=colorize(work.meta,resname=(r if r!='PtdIns' else 'PI2P'),),label=r),r) 
				for r in [a for a in allres if a not in resnames_pip2]+['PtdIns']])
	picturesave('fig.%s.simple'%plotname,work.plotdir,backup=False,
		version=True,meta={},extras=outer_labels)
	plt.close()

if 'simple_plot' in routine:

	"""
	new plotting method whereby we prepare the data ahead of time
	"""

	figsize = (10,10)
	fsbase = 20

	#---previously we were looking at pip2-chl1 bonds but those were overcounted 
	#---...(because it included chl1-chl1 accidentally). since the new counts were low, 
	#---...let's look at pip2-pip2
	alt_rep = ['pip2_chl1','pip2_pip2'][-1]

	#---selection
	sns = work.vars['orders']['canon']['asymmetric']
	btypes = ['lipid_water','lipid_pip2',alt_rep]
	#---switch to PIP2-water bonds
	btypes = ['lipid_water_pip2','lipid_pip2',alt_rep]

	axes,fig = panelplot(figsize=(16,12),layout={'out':{'grid':[1,2]},
		'ins':{'grid':[1,1]}})
	patches = []

	ylabels = {
		'lipid_water':'lipid-solvent hydrogen bonds',
		'lipid_pip2':'$\mathrm{{PIP}_{2}}$-lipid hydrogen bonds',
		'lipid_water_pip2':'$\mathrm{{PIP}_{2}}$-water hydrogen bonds',
		}

	extra_legends = []
	#---turned off the gray legend to simplify the plot
	if False:
		extra_legends.append([
			mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,alpha=0.5,fc="gray",lw=1,edgecolor='k',hatch='//'),
			#---after switching the alt_rep from CHOL to PIP2, we changed the text
			'$\mathrm{{PIP}_{2}-\mathrm{PIP}_{2}}$',
			])

	#---NEED A MORE PHOTOGENIC PLOT

	for bb,btype in enumerate(btypes):

		pnum = [(0,0),(1,0),(1,0)][bb]
		prepped = dict([(sn,{}) for sn in sns])
		for ss,sn in enumerate(sns):
			prepped[sn]['mean'] = dict(postdat)[sn][(btype,'mean')]
			prepped[sn]['std'] = dict(postdat)[sn][(btype,'std')]

		ax = axes[pnum[0]][pnum[1]]
		outs = ppi_bar_stylized_panel(ax,
			name='fig.%s.redux1'%plotname,
			data=prepped,
			sns=sns,
			zero_low=False,
			ylabel_format=lambda x:'%d'%x,
			ylabel='hydrogen bonds',
			legend=(bb==len(btypes)-1),
			altback=None,#"gray" if btype==alt_rep else "None",
			extra_legends=extra_legends if btype==alt_rep else None,
			std_dev_ec='k',#'k' if btype!='pip2_pip2' else 'gray',
			errorbars_only=True if btype==alt_rep else False,
			)
		patches.extend(outs.get('patches',[]))
		if outs.get('legend',None): patches.append(outs['legend'])
		name = ylabels.get(btype,None)
		if name: ax.set_title(name,fontsize=fsbase+2)

	picturesave('fig.%s.redux2'%plotname,work.plotdir,backup=False,
		version=True,meta={},extras=patches,pdf=True)
	plt.close()

#---show the good side as long as partner_1_bead, partner_2_bead and me (everything) are defined
#---! add the vector projection into vmdmake
good_side_deprecated = """
$me moveby [vecscale -1.0 [expr [$partner_2_bead get {x y z}]]]
display resetview
mouse stoprotation
rotate x to -90
set vec_view_target [expr [$partner_2_bead get {x y z}]] [expr [$partner_1_bead get {x y z}]]
set vec_up {0.0 0.0 1.0}
set vec_x {1.0 0.0 0.0}
set vec_look [veccross $vec_view_target $vec_up]
set tmat [transvecinv $vec_look]
$me move $tmat
set vec_look2 [veccross $vec_look $vec_x]
set tmat2 [transvecinv $vec_look2]
$me move $tmat2
"""

good_side = """
# project a vector to a certain plane
# Usage: projection $v1 x ==> project vector v1 to y-z plane
proc projection { vec axis } {
if {$axis == "x"} {
set vec "0 [lindex $vec 1] [lindex $vec 2]"
} elseif {$axis == "y"} {
set vec "[lindex $vec 0] 0 [lindex $vec 2]"
} else {
set vec "[lindex $vec 0] [lindex $vec 1] 0"
}
return $vec
}

#---move to partner_2_bead, which will be the axis
$me moveby [vecscale -1.0 [expr [$partner_2_bead get {x y z}]]]

#---reset the view
display resetview
mouse stoprotation
rotate x to -90

#---project the vector onto the plane and perform the rotation
set vec_main [vecsub [expr [$partner_2_bead get {x y z}]] [expr [$partner_1_bead get {x y z}]]]
set vec_main_flat [projection $vec_main z]
set vec_up {0.0 0.0 1.0}
set vec_look [veccross $vec_main_flat $vec_up]
set tmat [transvecinv $vec_look]
$me move $tmat

#---tie the whole thing together because we just rotated into the wrong direction
rotate y by 90
"""

if 'hbond_snaps' in routine:

	###---BEGIN SNAPSHOT PROCEDURE

	"""
	note that the main view = vmdmake block is taken from plot-head_angle.py
	so is the beginning directory creation and imports section, and the form of the snapshot namer 
	pseudocode:
		prepare a directory -- make sure it doesn't exist
		decide which snapshots to make (overwrite decisions, etc)
		run the snapshotter once
		...
		loop-ify it
	"""

	#---versioning the layouts
	tag = 'v1'
	sns_group = work.sns()

	def snapshot_namer(**kwargs):
		"""Name the snapshots in the plot folder."""
		canon_order = ['tag','sn','index','code','rank']
		mark_out,mark_in = '.','_'
		suffix = mark_out.join(
			'%s%s%s'%(str(key),mark_in,str(val))
			for key,val in [(k,kwargs[k]) for k in canon_order])
		return 'fig.hydrogen_bonding.snap.'+suffix

	#---precache identities in case they are already rendered
	if 'snap_plan' not in globals() and False:
		snap_plan = {}

		#---loop over simulations
		for sn_this in sns_group:
			#---print the top N most populated hydrogen bonds			
			for rank in range(n_most_popped_bonds):
				result = data[sn_this]['data']
				#---saved the PIP2-lipid bonds specifically
				#---note that this is a single list over all frames, unpackable by pip2_non_sol_hbonds_bpf
				hbonds = result['pip2_non_sol_hbonds']
				resnum_index = 5

				if len(hbonds)==0: continue

				#---! intervene here to redefine hbonds by name, resid
				#hbonds = hbonds[:,np.array([0,1,4,5])]
				#resnum_index = 3

				#---unique-ify the hbonds and find the most prominent ones
				#---since it's unique, we don't care about frames -- there can be no multiples
				try: hbonds_coded = np.ascontiguousarray(hbonds).view(
					np.dtype((np.void,hbonds.dtype.itemsize*hbonds.shape[1])))
				#---some simulations have no bonds
				except: continue
				hbonds_unique_coded,counts = np.unique(hbonds_coded,return_counts=True)
				#---reconstitute the unique strings as tuples
				_,idx,counts_not_sorted = np.unique(hbonds_coded,return_index=True,return_counts=True)
				hbonds_unique = hbonds[idx]
				#---select the first most-populated hbond with two lipids
				hbonds_unique_sort = hbonds_unique[np.argsort(counts_not_sorted)[::-1]]
				counts = counts_not_sorted[np.argsort(counts_not_sorted)[::-1]]
				#---! ryan and david manually checked that the two above variables are indexed together
				index_of_type = [ii for ii,i in enumerate(hbonds_unique_sort) if i[1]!=i[resnum_index]][rank]
				#---check the counts
				print('[NOTE] selecting bond type %s'%str(hbonds_unique_sort[index_of_type]))
				print('[NOTE] this type is observed in %d frames'%counts[index_of_type])
				#---get the first occurence of this bond type in the trajectory
				index = [ii for ii,i in enumerate(hbonds) if 
					np.all(i==hbonds_unique_sort[index_of_type])][rank]
				#---get the frame number from the list of bonds per frame
				frame = np.where(index<np.cumsum(result['pip2_non_sol_hbonds_bpf']))[0][0]
				#---get resid and resnames
				if not np.all(hbonds_unique_sort[index_of_type]==hbonds[index]): raise Exception('problem')
				#rn1,ri1,an1,hname2,rn2,ri2,an2 = hbonds[index]
				snap_plan[(sn_this,rank)] = list(hbonds[index])+[frame]+[index]
				
	#---single tests
	if False:
		singleout = ('membrane-v530',0)
		snap_plan = {singleout:snap_plan[singleout]}

	#---store the snapshots in the post_plot_spot according to the tag
	tempdir = os.path.join(work.paths['post_plot_spot'],'fig.hydrogen_bonding.%s'%tag)
	#---only make figures if the folder is empty. it only takes a few minutes -- no more coding logic
	if not os.path.isdir(tempdir): 
		os.mkdir(tempdir)
		status('snapshots dropping to %s (delete them if you want to re-make them)'%tempdir,tag='note')

		from codes.vmdmake import vmdmake
		import matplotlib.image as mpimg
		from matplotlib.path import Path
		import matplotlib.patches as patches
		import tempfile

		for sn_this,rank in snap_plan:

			rn1,ri1,an1,hname2,rn2,ri2,an2,frame,index = snap_plan[(sn_this,rank)]
			code = rn1,ri1,an1,hname2,rn2,ri2,an2
			snapshot_fn = snapshot_namer(tag=tag,sn=sn_this,index=index,
				code='_'.join(code),rank=rank)
			#---move the cursor back
			work.cursor = ('ptdins','xtc')
			gro = work.slice(sn_this)['current']['all']['gro']
			xtc = work.slice(sn_this)['current']['all']['xtc']
			tpr = work.get_last_start_structure(sn_this,part_name='tpr')
			kwargs_vmdmake = {'CGBONDSPATH':'~/libs/cg_bonds.tcl','GMXDUMP':'~/libs/gmxdump'}

			#---run VMD to make the snapshots via vmdmake
			view = vmdmake.VMDWrap(site=tempdir,
				gro=work.paths['post_data_spot']+gro,xtc=work.paths['post_data_spot']+xtc,tpr=tpr,
				frames='',res=(4000,4000),**kwargs_vmdmake)
			view.do('load_dynamic','standard','bonder')
			view.command("animate goto %d"%frame)
			#---show both lipids
			selections = []
			selections.append("(resname %s and resid %s) and not hydrogen"%(rn1,ri1))
			view.select(**{'partner_1':selections[-1],'style':
				'Licorice 0.200000 12.000000 12.000000','goodsell':True})
			selections.append("(resname %s and resid %s) and not hydrogen"%(rn2,ri2))
			view.select(**{'partner_2':selections[-1],
				'style':'Licorice 0.200000 12.000000 12.000000','goodsell':True})
			#---beads for the acceptor/donor
			selections.append("((resname %s and resid %s and name %s) and not hydrogen)"%(rn1,ri1,an1))
			view.select(**{
				'partner_1_bead':selections[-1],
				'style':'VDW 0.400000 12.000000','goodsell':True})
			#---beads for the acceptor/donor
			selections.append("((resname %s and resid %s and name %s) and not hydrogen)"%(rn2,ri2,an2))
			view.select(**{
				'partner_2_bead':selections[-1],
				'style':'VDW 0.400000 12.000000','goodsell':True})
			#---bead for the hydrogen (I think residue 2 is always the donor but this makes it sure)
			selections.append("(resid %s or resid %s) and name %s"%(ri1,ri2,hname2))
			view.command('set me [atomselect top "%s"]'%' or '.join(['(%s)'%i for i in selections]))
			view.select(**{
				'hydrogen_bead':selections[-1],
				'style':'VDW 0.400000 12.000000','goodsell':True})
			#---move everything
			#---draw a line for the hydrogen bond because the hydrogen bond rep is too permissive 
			#---...and dashed is hard to see
			#---below is a method for hydrogen bonding but VMD is too permissive 
			#---...and allows carbon to be the donor !?
			if False:
				#---show both in case you want to do a hydrogen bond representation
				view.select(**{'both_partners':
					"((resname %s and resid %s)) or ((resname %s and resid %s))"%(rn1,ri1,rn2,ri2),
					'style':'HBonds 4.000000 30.000000 3.000000'})
				view.command("mol modselect 2 0 ((resname %s and resid %s)) or ((resname %s and resid %s))"%
					(rn1,ri1,rn2,ri2))
				view.command("mol modstyle 2 0 HBonds 3.000000 20.000000 4.000000")
				view.command("mol modcolor 2 top ColorID 16")
			view['snapshot_filename'] = snapshot_fn
			view.command(good_side)
			view.command("draw cylinder [expr [$partner_2_bead get {x y z}]] "+
				"[expr [$partner_1_bead get {x y z}]] radius 0.1")
			view.do('snapshot')
			view.show(quit=True)

	#---MAKE A PANEL PLOT OF ALL OF THE THUMBNAILS
	#---several functions poached from plot-head_angle.py

	import matplotlib.image as mpimg
	from matplotlib.path import Path
	import matplotlib.patches as patches
	import tempfile

	def get_blank_border(img):
		"""Return the border limits for removing whitespace."""
		nonblank = np.any(img!=1.0,axis=2)*1
		lims = [map(lambda x:(x[0]-1,x[-1]+1),
			[np.where(np.any(nonblank!=0,axis=j))[0]])[0] for j in range(2)]
		return lims

	def sidestack_images(fns):
		"""Combine snapshots into one image with same zoom and minimal whitespace."""
		imgs = [mpimg.imread(os.path.join(tempdir,fn+'.png')) for fn in fns]
		borders = np.array([get_blank_border(img) for img in imgs])
		lims = np.array([[i.min() for i in borders.T[0]]]+[[i.max() for i in borders.T[1]]]).T
		#---stack horizontally (so the second column of centers is the only relevant one)
		buffer_width = 100
		imgs_zoomed = [img[slice(*lims[1]),slice(*lims[0])] for img in imgs]
		imgs_cropped = [i[:,slice(*(get_blank_border(i)[0]))] for i in imgs_zoomed]
		buffer_strip = np.ones((imgs_cropped[0].shape[0],buffer_width,3))
		#---add a vertical buffer strip
		imgs_combo = np.concatenate([np.concatenate((i,buffer_strip),axis=1) 
			for i in imgs_cropped[:-1]]+[imgs_cropped[-1]],axis=1)
		#---return the centers as well (assume horizontal here)
		widths = np.array([i.shape[:2] for i in imgs_cropped])[:,1]
		seps = np.cumsum(np.concatenate([[0]]+[[i,buffer_width] for i in widths]))[:-1]
		centers = ((seps[1:]+seps[:-1])/2.0)[::2]
		return imgs_combo,centers

	fsbase = 18
	stagger_labels = False
	extra_stagger = 250
	tag_drop = 100

	for sn_this in list(set(zip(*snap_plan.keys())[0])):

		snapshot_fns = [snapshot_namer(tag=tag,sn=sn_this,index=snap_plan[(sn_this,r)][-1],
			code='_'.join(snap_plan[(sn_this,r)][:-2]),rank=r) for r in range(n_most_popped_bonds)]

		axes,fig = panelplot(figsize=(12,12),layout={'out':{'grid':[1,1]},
			'ins':[{'grid':[1,1],'wspace':0.1}]})

		patchlist = []
		ax = axes[0]

		img,centers = sidestack_images(snapshot_fns)

		tagbox = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.3")
		for cc,c in enumerate(centers):
			text = '%s-%s'%(
				work.vars['names']['proper_residue_names_long'][snap_plan[(sn_this,cc)][0]],
				snap_plan[(sn_this,cc)][4])
			#---use the modulo operation to stagger but this is deprecated in favor of a newline
			patchlist.append(ax.text(c,img.shape[0]+tag_drop+stagger_labels*extra_stagger*(cc%2==0),
				text,bbox=tagbox,ha='center',va='top',rotation=0,
				color='k',fontsize=fsbase+2))

		ax.imshow(img)

		#---clean up snapshot axes
		for aa,ax in enumerate(axes):
			#-----------------------vvvvvvv
			ax.set_ylabel(work.meta[sn_this]['ion_label'],fontsize=fsbase+8,labelpad=60,rotation=0)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_xticks([]);ax.set_yticks([])
			[ax.spines[i].set_visible(False) for i in ['top','bottom','left','right']]

		picturesave('fig.hydrogen_bonding_snap.temp_summary.%s'%sn_this,
			work.plotdir,backup=False,version=False,meta={},extras=patchlist+[])
		plt.close()

"""
final notes: 
	need to remove stray hydrogen
	need to figure out where the cholesterol are ... do the unique lipids by residue name and not resname/id combos
confirming pip2-pip2 works the same in advanced vs basic indexing whatever
	>>> counts_pip2_pip2.mean()
	46.97627965043695
	>>> counts[np.where(np.sum((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,2]=='PI2P'),axis=0)==2)].sum()/float(nframes)
	46.97627965043695
more
	>>> counts[np.where(np.sum((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,4]=='PI2P'),axis=0)==1)].sum()/float(nframes)
	10.414481897627965
	>>> counts_lipid_pip2.mean()
	157.2796504369538
	>>> 
INCONSISTENCY
	>>> counts[np.where(np.sum((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,4]=='PI2P'),axis=0)==1)].sum()/float(nframes)
	10.414481897627965
	>>> counts_lipid_pip2.mean()
	157.2796504369538
	>>> 
INCONSISTENCY
	>>> bl[any_pip2][0]
	array(['CHL1', 'O3', 'PI2P', 'O11'], 
	      dtype='|S4')
	>>> np.where(np.all((hbonds[:,4]=='CHL1',hbonds[:,6]=='O3',hbonds[:,0]=='PI2P',hbonds[:,2]=='O11'),axis=0))
	(array([], dtype=int64),)
	>>> np.where(np.all((hbonds[:,0]=='CHL1',hbonds[:,2]=='O3',hbonds[:,4]=='PI2P',hbonds[:,6]=='O11'),axis=0))
	(array([], dtype=int64),)

DEBUGGING
	the problem above demonstrates that the bond list for items with pip2 has a bond that is not present in the more detailed data dump. we need to locate the CHL1 to PIP2-O11 bond
	NOTES FROM hydrogen_bonding.py
		reduced_tabs[np.where(np.all((np.any((reduced_tabs[:,0]=='CHL1',reduced_tabs[:,2]=='CHL1'),axis=0),np.any((reduced_tabs[:,0]=='PI2P',reduced_tabs[:,2]=='PI2P'),axis=0)),axis=0))[0]]
	PROVING THE BUG WAS FIXED
		>>> counts[np.where(np.sum((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,4]=='PI2P'),axis=0)==1)].sum()/float(nframes)+counts[np.where(np.sum((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,4]=='PI2P'),axis=0)==2)].sum()/float(nframes)
		157.2796504369538
		>>> counts_lipid_pip2.mean()
		157.2796504369538
comparing chl1 to PI2P
	counts[np.where(np.all((np.any((hbonds_unique_sort[:,0]=='PI2P',hbonds_unique_sort[:,4]=='PI2P'),axis=0),np.any((hbonds_unique_sort[:,0]=='CHL1',hbonds_unique_sort[:,4]=='CHL1
'),axis=0)),axis=0))[0]].sum()/float(nframes)
"""

### ABOVE DETRITUS saved for posterity after ryan and david debug on 26th and finished by 28th
### ca 2017.2.24 before meeting ryan tries to make the figure he said he would make

if 'simple_plot_breakdown' in routine:
	"""
	new plotting method whereby we prepare the data ahead of time
	"""
	#---selection
	sns = work.vars['orders']['canon']['asymmetric']
	sn = sns[0]
	bonds_tabbed[sn]

	figsize = (10,10)
	fsbase = 20
	#---previously we were looking at pip2-chl1 bonds but those were overcounted 
	#---...(because it included chl1-chl1 accidentally). since the new counts were low, 
	#---...let's look at pip2-pip2
	alt_rep = ['pip2_chl1','pip2_pip2'][-1]
	#---selection
	sns = work.vars['orders']['canon']['asymmetric']
	btypes = ['lipid_water','lipid_pip2',alt_rep]
	#---switch to PIP2-water bonds
	btypes = ['lipid_water_pip2','lipid_pip2',alt_rep]
	axes,fig = panelplot(figsize=(12,12),layout={'out':{'grid':[1,1]},
		'ins':{'grid':[1,1]}})
	patches = []
	ylabels = {
		'lipid_water':'lipid-solvent hydrogen bonds',
		'lipid_pip2':'$\mathrm{{PIP}_{2}}$-lipid hydrogen bonds',
		'lipid_water_pip2':'$\mathrm{{PIP}_{2}}$-water hydrogen bonds'}
	extra_legends = []
	#---NEED A MORE PHOTOGENIC PLOT
	for bb,btype in enumerate(btypes):
		pnum = [(0,0),(1,0),(1,0)][bb]
		prepped = dict([(sn,{}) for sn in sns])
		for ss,sn in enumerate(sns):
			prepped[sn]['mean'] = dict(postdat)[sn][(btype,'mean')]
			prepped[sn]['std'] = dict(postdat)[sn][(btype,'std')]
		ax = axes[pnum[0]][pnum[1]]
		outs = ppi_bar_stylized_panel(ax,
			name='fig.%s.redux1'%plotname,
			data=prepped,
			sns=sns,
			zero_low=False,
			ylabel_format=lambda x:'%d'%x,
			ylabel='hydrogen bonds',
			legend=(bb==len(btypes)-1),
			altback=None,#"gray" if btype==alt_rep else "None",
			extra_legends=extra_legends if btype==alt_rep else None,
			std_dev_ec='k',#'k' if btype!='pip2_pip2' else 'gray',
			errorbars_only=True if btype==alt_rep else False)
		patches.extend(outs.get('patches',[]))
		if outs.get('legend',None): patches.append(outs['legend'])
		name = ylabels.get(btype,None)
		if name: ax.set_title(name,fontsize=fsbase+2)
	picturesave('fig.%s.redux3'%plotname,work.plotdir,backup=False,
		version=True,meta={},extras=patches,pdf=True)
	plt.close()

if 'pip2_only_breakdown' in routine:

	"""
	In previous versions of the plot we were overcounting pip2-other hbonds because we included chl1-chl1
	bonds. This was fixed in the redux v2 and v3 plots. But these are somewhat awkward because they included
	a top bar for pip2 and lower bars that included an alt_rep that was pip2-chl2 below it.
	One one of the plots I switched the alternate representation to show pip2-pip2 bonds instead of the 
	much sparser pip2-chl1 bonds.
	New idea is to just do groups of bars to handle inter vs intra bonds, and pip2-chl1, etc.
	"""

	fsbase = 20
	alt_rep = ['pip2_chl1','pip2_pip2'][-1]
	#---selection
	sns = work.vars['orders']['canon']['asymmetric_pi_first']
	btypes = ['lipid_water','lipid_pip2',alt_rep]
	#---switch to PIP2-water bonds
	btypes = ['lipid_water_pip2','lipid_pip2',alt_rep]
	btypes = ['lipid_pip2',alt_rep,'lipid_water_pip2','lipid_chl1','pip2_chl1']
	axes,fig = panelplot(figsize=(16,16),layout={'out':{'grid':[3,3],'hspace':0.6},
		'ins':{'grid':[1,1]}})
	axes = [i for j in axes for i in j]
	patches = []
	ylabels = {
		'lipid_water':'lipid-solvent',
		'lipid_pip2':r'$\mathrm{{PIP}_{2}-\text{lipid}}$',
		'lipid_water_pip2':r'$\mathrm{{PIP}_{2}-\text{water}}$',
		'pip2_pip2':r'$\mathrm{{PIP}_{2}}-\mathrm{{PIP}_{2}}$',
		'pip2_chl1':r'$\mathrm{{PIP}_{2}}-\mathrm{CHOL}$',
		'lipid_chl1':r'$\mathrm{CHOL}-\mathrm{lipid}$'}
	extra_legends = []
	for bb,btype in enumerate(btypes):
		prepped = dict([(sn,{}) for sn in sns])
		for ss,sn in enumerate(sns):
			try: 
				prepped[sn]['mean'] = dict(postdat)[sn][(btype,'mean')]
				prepped[sn]['std'] = dict(postdat)[sn][(btype,'std')]
			except:
				prepped[sn]['mean'] = prepped[sn]['std'] = 0
		ax = axes[bb]
		outs = ppi_bar_stylized_panel(ax,
			name='fig.%s.redux1'%plotname,
			data=prepped,
			sns=sns,
			zero_low=False,
			ylabel_format=lambda x:'%d'%x,
			ylabel='hydrogen bonds' if bb in [0,3] else None,
			legend=bb==4,ncol=2,
			extra_legends=extra_legends if bb==2 else None)
		patches.extend(outs.get('patches',[]))
		if outs.get('legend',None): patches.append(outs['legend'])
		name = ylabels.get(btype,None)
		if name: ax.set_title(name,fontsize=fsbase+2)

	#---! the following wastes lines of code
	#---last plot is the self-other bonds, but we have no variation on that ...
	ax = axes[len(btypes)+1]
	pip2_self_other = {}
	for sn in sns:
		tabs = bonds_tabulate(sn,{'resname':work.meta[sn]['ptdins_resname']},
			{'resname':work.meta[sn]['ptdins_resname']})
		bonds,counts = [tabs[j] for j in ['bonds','counts']]
		if len(bonds)>0:
			nframes = float(len(data[sn]['data']['variations']))
			count = np.sum(counts[np.where(bonds[:,1]!=bonds[:,5])])/nframes
			pip2_self_other[sn] = {'mean':count,'std':0}
		else: pip2_self_other[sn] = {'mean':0,'std':0}
	outs = ppi_bar_stylized_panel(ax,
		name='fig.%s.redux1'%plotname,
		data=pip2_self_other,
		sns=sns,
		zero_low=False,
		ylabel_format=lambda x:'%d'%x,
		ylabel='hydrogen bonds',
		legend=False,
		extra_legends=None)
	ax.set_title('$\mathrm{{PIP}_{2}}-\mathrm{lipid}$ (inter)',fontsize=fsbase+2)

	#---! the following wastes lines of code
	#---last plot is the self-other bonds, but we have no variation on that ...
	ax = axes[len(btypes)+2]
	pip2_self_other = {}
	for sn in sns:
		tabs = bonds_tabulate(sn,{'resname':work.meta[sn]['ptdins_resname']},
			{'resname':work.meta[sn]['ptdins_resname']})
		bonds,counts = [tabs[j] for j in ['bonds','counts']]
		if len(bonds)>0:
			nframes = float(len(data[sn]['data']['variations']))
			count = np.sum(counts[np.where(np.all((bonds[:,1]!=bonds[:,5],
				np.all((bonds[:,0]==work.meta[sn]['ptdins_resname'],
					bonds[:,4]==work.meta[sn]['ptdins_resname']),axis=0)),axis=0))])/nframes
			pip2_self_other[sn] = {'mean':count,'std':0}
		else: pip2_self_other[sn] = {'mean':0,'std':0}
	outs = ppi_bar_stylized_panel(ax,
		name='fig.%s.redux1'%plotname,
		data=pip2_self_other,
		sns=sns,
		zero_low=False,
		ylabel_format=lambda x:'%.1f'%x,
		ylabel=None,
		legend=False,
		extra_legends=None)
	ax.set_title('$\mathrm{{PIP}_{2}}-\mathrm{{PIP}_{2}}$ (inter)',fontsize=fsbase+2)

	fig.delaxes(axes[5])
	fig.delaxes(axes[8])

	picturesave('fig.%s.redux4'%plotname,work.plotdir,backup=False,
		version=True,meta={},extras=patches,pdf=True)
	plt.close()

	#---CHECKING SOME STUFF
	#---use bonds_tabbed to get PIP2-other hydrogen bonds that excludes internal ones
	sn = 'membrane-v530'
	counts,bonds = [bonds_tabbed[sn][j] for j in ['counts','bonds']]
	#---both PIP2 and not the same residue number
	np.sum(counts[np.where(np.all((bonds[:,1]!=bonds[:,5],bonds[:,0]=='PI2P',bonds[:,4]=='PI2P'),axis=0))])/nframes
	dict(postdat)[sn][(btype,'mean')]
	"""
	ANOTHER SERIOUS FUCKING PROBLEM.
	>>> dict(postdat)[sn][(btype,'mean')]
	46.97627965043695
	>>> counts,bonds = [bonds_tabbed[sn][j] for j in ['counts','bonds']]
	>>> np.sum(counts[np.where(np.all((bonds[:,1]!=bonds[:,5],bonds[:,0]=='PI2P',bonds[:,4]=='PI2P'),axis=0))])/nframes
	2.3250000000000002
	>>> np.sum(counts[np.where(np.all((bonds[:,1]==bonds[:,5],bonds[:,0]=='PI2P',bonds[:,4]=='PI2P'),axis=0))])/nframes
	76.066666666666663
	"""

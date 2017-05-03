#!/usr/bin/python -i

#---plot prep
if False:
	if 'plotload' not in globals(): execfile('/etc/pythonstart')
	execfile('./omni/base/header.py')
	from plotter import *
	from base.store import plotload
	execfile('./calcs/specs/figures.py')
	execfile('./calcs/specs/colors.py')

art = {
	'fs':{'legend':10},
	}

#---settings
plotname = 'ion_binding'
######### 
if 'data' not in globals(): data,calcs = plotload(plotname,work)
if False:
	if 'data' not in globals():
		calc_name = 'ion_binding'
		data = {}
		data[calc_name] = {}
		sns = work.sns()
		sn = sns[0]	
		data[calcname][sn] = work.plotload_manual(calcname=calc_name,specs={})
		print('DONE')
		import ipdb;ipdb.set_trace()

sns = data['ion_binding'].keys()
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib import patches
import brewer2mpl
clrs = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors
import numpy as np
ns = sys.argv[0]

#---settings
max_by_ion = True
onemax = 700
bars = False

#---tasks
routine = ['plot','ion_binding_redux'][-1:]

#---PLOTS
#-------------------------------------------------------------------------------------------------------------

if 'plot' in routine:

	print_version = 'v20160929'

	layout_name = 'summary2'
	bar_suffix = '_farbar' if bars else ''
	zone_cutoffs = sorted([c['specs']['zonecut'] for c in calcs['ion_binding_combinator']])
	for cnum,zonecut in enumerate(zone_cutoffs):
		axes,fig = panelplot(layout=figlayout[layout_name],figsize=(12,16))
		countmax = [0 for i in figlayout[layout_name]['ins']]
		for pnum,sn in enumerate(sns):
			axrow,axcol = figplacer(sn,figplace[layout_name])
			ax = axes[axrow][axcol]
			ax.set_ylabel('bound ions')
			ax.set_xlabel('time (ns)')
			ax.set_title(' '.join([
				work.meta[sn]['composition_name'],
				work.meta[sn]['ptdins_label'],
				work.meta[sn]['ion_label'],
				]))
			ts_skip = work.slice(sn)[calcs['ion_binding']['slice_name']][
				calcs['ion_binding']['group']]['skip']		
			datb = data['ion_binding_combinator'][cnum][sn]['data']
			nframes = len(datb[i2s(0,'wcs')])
			base = np.zeros(nframes)
			for nn in range(3):
				totals = np.sum(datb[i2s2(nn,'wcs')],axis=1)
				ax.fill_between(np.arange(nframes)*ts_skip/1000.,base,totals+base,
					color=clrs[nn],alpha=0.5,lw=0)
				involves_ptdins = np.sum(datb[i2s2(nn,'wcs')][:,np.array([ii 
					for ii,i in enumerate(datb[i2s2(nn,'combonames')]) 
					if work.meta[sn]['ptdins_resname'] in i])],axis=1)
				ax.fill_between(np.arange(nframes)*ts_skip/1000.,base,involves_ptdins+base,
					color=clrs[nn],alpha=1.,lw=0)
				if nn == 2: ax.plot(np.arange(nframes)*ts_skip/1000.,totals+base,lw=1,c='k')
				base += totals
				countmax[axrow] = max([countmax[axrow],max(totals+base)])
			if max_by_ion: 
				ax.set_ylim(0,700 if work.meta[sn]['cation']=='NA' or 
				sn in ['membrane-v542','membrane-v543'] else 400)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_xlim(0,nframes*ts_skip/1000.)
		if not max_by_ion:
			for ii,axrow in enumerate(axes):
				for ax in axrow: ax.set_ylim(0,topcount)
		for rr,cc in [(ii,jj) for ii,i in enumerate(figplace[layout_name]) 
			for jj,j in enumerate(i) if j==None]:
			fig.delaxes(axes[rr][cc])
		descriptor = 'cutoff%.1f'%zonecut
		ax = axes[0][1]
		status_bar_patches = []
		if bars:
			status_bar_patches.append(
				patches.Rectangle((0.0,1.6),1.0,0.08,transform=ax.transAxes,
				facecolor='w',edgecolor='k',clip_on=False,lw=1.5))
			status_bar_patches.append(
				patches.Rectangle((0.0,1.6),
				float(cnum)/len(zone_cutoffs),0.08,transform=ax.transAxes,
				facecolor='k',edgecolor='k',clip_on=False,lw=1.5))
			for i in status_bar_patches: ax.add_patch(i)
			ax = axes[0][0]
			at = AnchoredText('$\mathrm{cutoff=%.1f\AA}$'%zonecut,
				prop=dict(size=20),frameon=True,loc=3,bbox_to_anchor=(0.5,1.3),
				bbox_transform=ax.transAxes)
			at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
			ax.add_artist(at)		
		picturesave('fig.%s.%s%s'%(plotname,descriptor,bar_suffix),work.plotdir,backup=False,
			version=True,meta={'zonecut':zonecut,'sns':sns,'print_version':print_version},
			extras=status_bar_patches)
		plt.close()

if 'debug_visually' in routine:

	import itertools
	sn = 'membrane-v532'
	dat = data['ion_binding']['membrane-v532']['data']
	pas = dat['partners_atoms']
	ld = lipid_distances = dat['lipid_distances']
	resnames = dat['resnames']
	lipids = array(list(resnames[sort(unique(resnames,return_index=True)[1])]))
	#---! big difference at frame 503
	#---! [data['ion_binding_combinator'][0]['membrane-v532']['data']['%d.wcs'%nn][503] for nn in range(3)]

	from codes.plotter3d import *
	from base.store import load
	mesh = load('v538.10000-90000-100.lipid_mesh.n1.dat',work.paths['post_data_spot'])

	meshpoints(mesh['%d.%d.points'%(0,502)],color=(1,0,1))
	meshpoints(mesh['%d.%d.points'%(0,503)])
	pbcwire(mesh['%d.%d.vec'%(0,502)])

	zonecut = 2.2
	nn = 2

	combos = array([''.join(j) for j in 
		itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) 
		if sum([int(k) for k in j])==nn+1])
	cind = lambda a : where(combos==''.join([str(sum(array(a)==i)) for i in lipids]))[0][0]
	for fr in [502,503]:
		parts = resnames[pas[fr,where(sum(lipid_distances[fr]<zonecut/10.,axis=1)==nn+1)[0]]][:,:nn+1]
		print len(parts)
		ans = array([sum(array([cind(j) for j in parts])==i) for i in range(len(combos))])
		print all(ans==array([data['ion_binding_combinator'][0]['membrane-v532']['data']['%d.wcs'%nn][fr] 
			for nn in range(3)][2]).astype(int))

	ax = plt.subplot(121)
	nums = sum(ld[...,0]<2.2/10.,axis=1)
	colors = 'rbg'
	eg = zeros((801))
	for nn in range(3):
		eg += sum(data['ion_binding_combinator'][0]['membrane-v532']['data']['%d.wcs'%nn],axis=1)	
		ax.plot(eg,c=colors[nn])		
	ax = plt.subplot(122)
	alls = sum(lipid_distances[...,0]<zonecut/10.,axis=1)
	ax.plot(alls,c='k')
	ax.set_ylim(0,alls.max()*1.1)
	plt.show()

if 'debug_ion_binding' in routine:

	import MDAnalysis
	from codes.binding import *
	from codes.plotter3d import *
	from base.store import load

	if 'ans' not in globals():

		sn = 'membrane-v532'
		lan = 'v532.20000-100000-100.all.pbcmol.lipid_abstractor.n0.dat'
		abstract = load(lan,work.paths['post_data_spot'])

		#---prepare universe	
		if 'uni' not in globals():

			grofile,trajfile = [work.slice(sn)['current']['all'][i] for i in ['gro','xtc']]
			uni = MDAnalysis.Universe(work.postdir+grofile,work.postdir+trajfile)
			nframes = len(uni.trajectory)
			#---MDAnalysis uses Angstroms not nm
			lenscale = 10.

			#---compute masses by atoms within the selection
			sel_lipids = uni.select_atoms(' or '.join('resname %s'%r 
				for r in work.vars['selectors']['resnames_lipid_chol']))
			sel_ions = uni.select_atoms(work.vars['selectors']['cations'])

			#---load lipid points into memory
			trajectory_ions = zeros((nframes,len(sel_ions),3))
			trajectory_lipids = zeros((nframes,len(sel_lipids),3))
			vecs = zeros((nframes,3))
			for fr in range(nframes):
				status('loading frame',tag='load',i=fr,looplen=nframes)
				uni.trajectory[fr]
				trajectory_lipids[fr] = sel_lipids.coordinates()/lenscale
				trajectory_ions[fr] = sel_ions.coordinates()/lenscale
				vecs[fr] = sel_lipids.dimensions[:3]/lenscale

			#monolayer_indices = kwargs['upstream']['lipid_abstractor']['monolayer_indices']
			monolayer_indices = abstract['monolayer_indices']
			#resids = kwargs['upstream']['lipid_abstractor']['resids']
			resids = abstract['resids']
			monolayer_residues = [resids[where(monolayer_indices==mn)[0]] for mn in range(2)]
			group_lipid = uni.select_atoms(' or '.join(['resid '+str(i) for mononum in range(2) 
				for i in monolayer_residues[mononum]]))
			lipid_resids = array([i.resid for i in group_lipid])
			if work.meta[sn]['composition_name'] != 'asymmetric': lipid_resid_subselect = slice(None,None)
			#---hack to account for asymmetric bilayer by analyzing only the first (top) monolayer 
			else: lipid_resid_subselect = where([i.resid in monolayer_residues[0] for i in group_lipid])[0]

	comparison = [502,503]
	comparison = [501,502]

	if 'err' not in globals(): err = {}
	nrank = 3
	for fr in comparison:
		if fr not in err:
			status('frame %d'%fr,tag='debug')
			ld,pas = partnerfinder(
				trajectory_lipids[fr],trajectory_ions[fr],vecs[fr],
				lipid_resids,nrank,includes=lipid_resid_subselect)
			err[fr] = {'ld':ld,'pas':pas}

	for fr in []:
		pas,ld = err[fr]['pas'],err[fr]['ld']
		meshpoints(trajectory_lipids[0][pas.T[0][where(ld<=0.22)[0]]],color=(1,1,1),scale_factor=1.)
		meshpoints(trajectory_lipids[0][pas.T[0][where(all([ld<0.46,ld>0.22],axis=0))[0]]],
			color=(0,0,1),scale_factor=1.)

	#---everything looks fine at lipid_distances
	zonecut = 0.22
	print "LIPID DISTANCES LOOK LIKE THEY HAVE CONSISTENT COUNTS"
	for fr in comparison:
		print cumsum([len(np.where(np.sum(err[fr]['ld']<zonecut,axis=1)==nn+1)[0]) for nn in range(3)])

	#---moving on to check binding combinator
	cnum,zonecut = list(enumerate([c['specs']['zonecut'] for c in calcs['ion_binding_combinator']]))[0]
	zonecut = zonecut/10.
	datb = data['ion_binding_combinator'][cnum][sn]['data']
	print "HERE IS THE PROBLEM FROM ION BINDING COMBINATOR"
	print [[sum(datb['%d.wcs'%nn],axis=1)[fr] for nn in [0,1,2]] for fr in comparison]
	print [sum([sum(datb['%d.wcs'%nn],axis=1)[fr] for nn in [0,1,2]]) for fr in comparison]

	#---reconstructing the combinator code here
	nn = 0
	all_combos = [datb['%d.combos'%nn] for nn in range(3)]
	resnames = data['ion_binding'][sn]['data']['resnames']
	lipids = array(list(resnames[sort(unique(resnames,return_index=True)[1])]))
	for fr in comparison:
		for nn in range(3):
			lipid_distances = err[fr]['ld']
			pas = err[fr]['pas']
			combos = all_combos[nn]
			#cind = lambda a : where(combos==''.join([str(sum(array(a)==i)) for i in lipids]))[0][0]
			#[where(combos==''.join([str(sum(array(a)==i)) for i in lipids]))[0][0] for a in parts]
			#print sum(array([sum(array([cind(j) for j in parts])==i) for i in range(len(combos))]))
			parts = resnames[pas[where(sum(lipid_distances<zonecut,axis=1)==nn+1)[0]]][:,:nn+1]
			count = sum(array([sum(array([where(combos==''.join([str(sum(array(a)==j)) 
				for j in lipids]))[0][0] for a in parts])==i) for i in range(len(combos))]))
			err[fr][('parts',nn)] = parts
			err[fr][('count',nn)] = count

	for fr in comparison:
		print [len(err[fr][('parts',nn)]) for nn in range(3)]
		print [err[fr][('count',nn)] for nn in range(3)]

	"""
	nn=1;fr=comparison[1];
	parts = resnames[pas[where(sum(err[fr]['ld']<zonecut,axis=1)==nn+1)[0]]][:,:nn+1];parts.T[0];
	len(parts);
	count = sum(array([sum(array([where(all_combos[nn]==''.join([str(sum(array(a)==j)) 
		for j in lipids]))[0][0] for a in parts])==i) for i in range(len(combos))]));count
	"""

	#---note an aborted attempt to "plot the new lipid distances" here possibly during 
	#---...debugging but incomplete
	#---construct a GIF with: "convert -delay 20 -loop 0 fig.ion_binding.cutoff*.png myimage.gif"

#---! superceded below for plotting different sets of simulations on one plot
if 'ion_binding_redux_superceded' in routine:

	"""
	Better looking plots for the charging curves.
	"""

	do_single = [False,2.2][0]

	#---from plot-lipid_areas2d
	colors_ions = {
		'NA':'green',
		'Na,Cal':'green',
		'MG':'red',
		'Cal':'blue',
		'K':'gray',
		}

	print_version = 'v20161110'

	#---sweep parameters
	zone_cutoffs = sorted([c['specs']['zonecut'] for c in calcs['ion_binding_combinator']])

	#---loop over all plots
	if do_single: zone_cutoffs = [do_single]
	for cnum,zonecut in enumerate(zone_cutoffs):

		#---aesthetics
		descriptor = 'cutoff%.1f'%zonecut
		picname = 'fig.ion_binding2.%s'%descriptor
		time_pts = 3
		ymin = 0
		fsbase = 16
		patches = []

		#---selections
		sns_sep = [work.vars['orders']['canon']['symmetric'],work.vars['orders']['canon']['asymmetric']]
		
		#---figure
		axes,fig = panelplot(figsize=(16,10),
			layout={'out':{'grid':[1,2],'wratios':[len(s) for s in sns_sep],'wspace':0.1},
			'ins':{'grid':[1,1]}})
		axes = [a[0] for a in axes]
		
		#---format count data up to the maximum number of nearest-neighbors (hardcoded to 3)
		ncats = 3
		postdat_sep = [{} for s in sns_sep] 
		for pnum,sns in enumerate(sns_sep):
			postdat = dict([(sn,{}) for sn in sns])
			for sn in sns:
				datb = data['ion_binding_combinator'][cnum][sn]['data']
				total = np.sum([datb['%d.wcs'%i].sum(axis=1) for i in range(ncats)],axis=0)
				postdat[sn]['total'] = total
				sel_frames = np.linspace(0,total.shape[0]-1,time_pts).astype(int)
				postdat[sn]['total_few'] = [total[i] for i in sel_frames]
				postdat[sn]['nn_count_few'] = [list(np.array([datb['%d.wcs'%i].sum(axis=1) 
					for i in range(ncats)]).T[j]) for j in sel_frames]
				involves_ptdins = [np.sum(datb[i2s2(nn,'wcs')][:,np.array([ii for ii,i in enumerate(
					datb[i2s2(nn,'combonames')]) if work.meta[sn]['ptdins_resname'] in i])],axis=1) 
					for nn in range(ncats)]
				postdat[sn]['nn_count_few_ptdins'] = [[involves_ptdins[nn][j] for nn in range(ncats)] for j in sel_frames]
				postdat[sn]['nn_count_few_not_ptdins'] = [[
					postdat[sn]['nn_count_few'][j][i]-postdat[sn]['nn_count_few_ptdins'][j][i] 
					for i in range(ncats)] for j in range(len(sel_frames))]
			postdat_sep[pnum] = postdat

		for pnum,sns in enumerate(sns_sep):

			postdat = postdat_sep[pnum]
			ax = axes[pnum]

			#---plot the total ions at a few time points
			barspecs_all = dict(lw=0)
			#---if we just set yvals to the totals, we get a standard bar plot
			yvals = [postdat[sn]['total_few'] for sn in sns]
			#---stacked bar plot from the totals by number of neighbors
			yvals = [postdat[sn]['nn_count_few'] for sn in sns]
			#---specify by pip2 representation in the nearest neighbors
			yvals = []
			for ss,sn in enumerate(sns):
				ysn = []
				for t in range(time_pts):
					ysnt = []
					for l in range(ncats):
						for key in ['nn_count_few_not_ptdins','nn_count_few_ptdins']:
							ysnt.append(postdat[sn][key][t][l])
					ysn.append(ysnt)
				yvals.append(ysn)
			#---bars are colored by ion
			barspecs_base = [{'color':colors_ions[work.meta[sn]['cation']]} for sn in sns]
			alpha_sweep = dict([(ii,i) for ii,i in enumerate(np.linspace(0.35,1,3))])
			#---could apply different patterns to each of the time_pts by adding e.g. alpha_sweep to the dict below
			barspecs = [[dict(**barspecs_base[ss]) for j in range(time_pts)] 
				for ss,sn in enumerate(sns)]
			color_blender = lambda c,alpha : tuple(list(mpl.colors.ColorConverter().to_rgb(c))+[alpha]) 
			#---add additional granular settings for stacked bar plots
			for ii in range(len(barspecs_base)): barspecs_base[ii].pop('color')
			#---TUNE THE STACK FOR MAXIMUM COMPREHENSIBILITY !!!
			barspecs_stack = [[[dict(
				alpha=alpha_sweep[k/2],
				#lw=2 if k%2==1 else 0,
				lw=0,
				#hatch='--' if k%2==1 else None,
				#edgecolor=color_blender(colors_ions[work.meta[sn]['cation']],alpha_sweep[k/2]),
				#color='gray' if k%2==0 else colors_ions[work.meta[sn]['cation']],
				edgecolor='w',
				hatch='///' if k%2==0 else None,
				**barspecs_base[ss]) 
				for k in range(ncats*2)] 
				for j in range(time_pts)] for ss,sn in enumerate(sns)]

			out_specs = barmaker(ax,yvals,gap=1.0,
				barspecs_all=barspecs_all,barspecs=barspecs,barspecs_stack=barspecs_stack)

			ymax = out_specs['ymax']
			#---if we are making a batch we fix the ymin and ymax
			if not do_single: ymax = [400,250][pnum]
			ax.set_ylim(ymin,ymax*1.1)

			#---copied from elsewhere
			tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.5")
			tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.5")
			tagbox_ptdins = tagbox_ion
		
			ax.set_xticks([])
			label_drop = (ymax-ymin)*0.2
			width = out_specs['width']
			xcenters = [np.mean(np.array(i)+out_specs['width']/2.0) for i in out_specs['xvals']]
			for snum,sn in enumerate(sns):
				text = work.meta[sn]['ion_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*-0.02-label_drop,text,
					bbox=tagbox_ion,ha='center',va='top',rotation=90 if len(text)>25 else 0,
					color='k',fontsize=fsbase)
				patches.append(tb)
				text = work.meta[sn]['ptdins_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*0.02-label_drop,text,
					bbox=tagbox_ptdins,rotation=90,ha="center",va="bottom",color='k',fontsize=fsbase)
				patches.append(tb)

			#---the caption should explain the panels
			if False: ax.set_title(work.meta[sn]['composition_name'],fontsize=fsbase+4)

		#---annotations after the last bar
		rtx = out_specs['xvals'][-1][-1]+out_specs['width']
		cat_tops = np.array([0]+list(np.cumsum(yvals[-1][-1])[1::2]))
		midys = (cat_tops[1:]-cat_tops[:-1])/2.0+cat_tops[:-1]
		top_ys = np.cumsum(yvals[-1][-1])[1::2]
		el = mpl.patches.Ellipse((2, -1), 0.5, 0.5)
		for ncat,y in enumerate(top_ys):
			color = '#bdbdbd'
			ann = ax.annotate('ions bound\nto %d lipid%s'%(ncat+1,{0:''}.get(ncat,'s')),
				xy=(rtx,y),xycoords='data',
				xytext=(35,0),textcoords='offset points',
				size=16,va="center",
				bbox=dict(boxstyle="round",fc=color,ec="none"),
				arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					fc=color,ec="none",patchA=None,patchB=el,relpos=(0.2, 0.5),))
			patches.append(ann)

		#---spine formatting
		for ax in axes:
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('none')

		for ax in axes:
			ax.set_yticklabels(['%d'%i for i in ax.get_yticks()],fontsize=fsbase+2)
		axes[0].set_ylabel('bound ions (within $\mathrm{\mathbf{%.1f} \AA}$)'%zonecut,fontsize=fsbase+4)

		#---custom legend on the second plot
		legend_labels,legend_patches = [],[]
		ion_names = ['K','NA','MG','Cal']	
		marks_sns = [sn for sn in [m[0] for m in [[sn for sn in sns if work.meta[sn]['cation']==i] 
			for i in ion_names] if m]]
		for nn,name in enumerate(ion_names):
			legend_labels.append(work.meta[marks_sns[nn]]['ion_label'])
			legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,
				facecolor=colors_ions[name],edgecolor='k',lw=2))
		legend_labels.append('bound to\nat least\n1 PtdIns')
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,color='gray'))
		legend_labels.append('not bound\nto PtdIns')
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,fc='w',edgecolor='gray',hatch='///'))
		legend = ax.legend(legend_patches,legend_labels,loc='upper left',fontsize=art['fs']['legend']+6,
			ncol=3,labelspacing=1.2,handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)

		#---render
		picturesave(picname,work.plotdir,backup=False,extras=patches,
			version=True,meta={'zonecut':zonecut,'sns':sns_sep[0]+sns_sep[1],'print_version':print_version},
			#---pad inches because the tagboxes get cut off, pdf because hatches must be thicker
			pad_inches=0.5,pdf=True)

	#---make the gif
	os.remove(work.plotdir+'/fig.ion_binding2.gif')
	os.system('convert -delay 20 -loop 0 %s %s'%(
		work.plotdir+'/fig.ion_binding2.cutoff*v1.png',
		work.plotdir+'/fig.ion_binding2.gif'))

if 'ion_binding_redux' in routine:

	"""
	Better looking plots for the charging curves.
	redeveloped from above
	"""

	#---NO MORE HACKING MAHHHH
	mpl.rcParams['hatch.linewidth'] = 3.0

	from plotter.extras import i2s2

	do_single = [False,2.2][-1]

	#---from plot-lipid_areas2d
	colors_ions = {
		'NA':'green',
		'Na,Cal':'green',
		'MG':'red',
		'Cal':'blue',
		'K':'gray',
		}

	print_version = 'v20161110'

	#---loop over all plots
	if do_single: zone_cutoffs = [do_single]
	else:
		#---sweep parameters
		zone_cutoffs = sorted([c['specs']['zonecut'] for c in calcs['ion_binding_combinator']])
		#---! above as a full loop in the plot spec and currently that is disabled FOR REDEVELOPMENT

	for cnum,zonecut in enumerate(zone_cutoffs):

		#---aesthetics
		descriptor = 'cutoff%.1f'%zonecut
		picname = 'fig.ion_binding3.%s'%descriptor
		time_pts = 3
		ymin = 0
		fsbase = 16
		patches = []

		#---selections
		sns_sep = [work.vars['orders']['canon']['symmetric'],work.vars['orders']['canon']['asymmetric']]
		sns = sns_sep[0]+sns_sep[1]
		figsize = (16,10)
		sns_sep = [work.vars['orders']['canon']['asymmetric']]
		sns = sns_sep[0]
		figsize = (12,12)

		#---figure
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,len(sns_sep)],'wratios':[len(s) for s in sns_sep],'wspace':0.1},
			'ins':{'grid':[1,1]}})
		#axes = [a[0] for a in axes]
		axes = axes

		
		#---format count data up to the maximum number of nearest-neighbors (hardcoded to 3)
		ncats = 3
		postdat_sep = [{} for s in sns_sep] 
		for pnum,sns in enumerate(sns_sep):
			postdat = dict([(sn,{}) for sn in sns])
			for sn in sns:
				#datb = data['ion_binding_combinator'][cnum][sn]['data']
				datb = data['ion_binding_combinator'][sn]['data']
				total = np.sum([datb['%d.wcs'%i].sum(axis=1) for i in range(ncats)],axis=0)
				postdat[sn]['total'] = total
				sel_frames = np.linspace(0,total.shape[0]-1,time_pts).astype(int)
				postdat[sn]['total_few'] = [total[i] for i in sel_frames]
				postdat[sn]['nn_count_few'] = [list(np.array([datb['%d.wcs'%i].sum(axis=1) 
					for i in range(ncats)]).T[j]) for j in sel_frames]
				involves_ptdins = [np.sum(datb[i2s2(nn,'wcs')][:,np.array([ii for ii,i in enumerate(
					datb[i2s2(nn,'combonames')]) if work.meta[sn]['ptdins_resname'] in i])],axis=1) 
					for nn in range(ncats)]
				postdat[sn]['nn_count_few_ptdins'] = [[involves_ptdins[nn][j] for nn in range(ncats)] for j in sel_frames]
				postdat[sn]['nn_count_few_not_ptdins'] = [[
					postdat[sn]['nn_count_few'][j][i]-postdat[sn]['nn_count_few_ptdins'][j][i] 
					for i in range(ncats)] for j in range(len(sel_frames))]
			postdat_sep[pnum] = postdat

		for pnum,sns in enumerate(sns_sep):

			postdat = postdat_sep[pnum]
			ax = axes[pnum]

			#---plot the total ions at a few time points
			barspecs_all = dict(lw=0)
			#---if we just set yvals to the totals, we get a standard bar plot
			yvals = [postdat[sn]['total_few'] for sn in sns]
			#---stacked bar plot from the totals by number of neighbors
			yvals = [postdat[sn]['nn_count_few'] for sn in sns]
			#---specify by pip2 representation in the nearest neighbors
			yvals = []
			for ss,sn in enumerate(sns):
				ysn = []
				for t in range(time_pts):
					ysnt = []
					for l in range(ncats):
						for key in ['nn_count_few_not_ptdins','nn_count_few_ptdins']:
							ysnt.append(postdat[sn][key][t][l])
					ysn.append(ysnt)
				yvals.append(ysn)
			#---bars are colored by ion
			barspecs_base = [{'color':colors_ions[work.meta[sn]['cation']]} for sn in sns]
			alpha_sweep = dict([(ii,i) for ii,i in enumerate(np.linspace(0.35,1,3))])
			#---could apply different patterns to each of the time_pts by adding e.g. alpha_sweep to the dict below
			barspecs = [[dict(**barspecs_base[ss]) for j in range(time_pts)] 
				for ss,sn in enumerate(sns)]
			color_blender = lambda c,alpha : tuple(list(mpl.colors.ColorConverter().to_rgb(c))+[alpha]) 
			#---add additional granular settings for stacked bar plots
			for ii in range(len(barspecs_base)): barspecs_base[ii].pop('color')
			#---TUNE THE STACK FOR MAXIMUM COMPREHENSIBILITY !!!
			barspecs_stack = [[[dict(
				#---drop the alpha on the hatches because now I cannot change the hatch color and it's black
				alpha=alpha_sweep[k/2] if not k%2==0 else alpha_sweep[k/2]*0.65,
				#lw=2 if k%2==1 else 0,
				lw=0,
				#hatch='--' if k%2==1 else None,
				#edgecolor=color_blender(colors_ions[work.meta[sn]['cation']],alpha_sweep[k/2]),
				color='w' if k%2==0 else colors_ions[work.meta[sn]['cation']],
				#edgecolor='k' if k%2==0 else colors_ions[work.meta[sn]['cation']],
				#---!ahhh they are fucking around with hatches again!
				edgecolor=color_blender(colors_ions[work.meta[sn]['cation']],alpha_sweep[k/2]),
				hatch='///' if k%2==0 else None,
				**barspecs_base[ss]) 
				for k in range(ncats*2)] 
				for j in range(time_pts)] for ss,sn in enumerate(sns)]

			out_specs = barmaker(ax,yvals,gap=1.0,
				barspecs_all=barspecs_all,barspecs=barspecs,barspecs_stack=barspecs_stack)

			ymax = out_specs['ymax']
			#---if we are making a batch we fix the ymin and ymax
			if not do_single: ymax = [400,250][pnum]
			ax.set_ylim(ymin,ymax*1.1)

			#---copied from elsewhere
			tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.5")
			tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.5")
			tagbox_ptdins = tagbox_ion
		
			ax.set_xticks([])
			label_drop = (ymax-ymin)*0.2
			width = out_specs['width']
			xcenters = [np.mean(np.array(i)+out_specs['width']/2.0) for i in out_specs['xvals']]
			for snum,sn in enumerate(sns):
				text = work.meta[sn]['ion_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*-0.02-label_drop,text,
					bbox=tagbox_ion,ha='center',va='top',rotation=90 if len(text)>25 else 0,
					color='k',fontsize=fsbase)
				patches.append(tb)
				text = work.meta[sn]['ptdins_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*0.02-label_drop,text,
					bbox=tagbox_ptdins,rotation=90,ha="center",va="bottom",color='k',fontsize=fsbase)
				patches.append(tb)

			#---the caption should explain the panels
			if False: ax.set_title(work.meta[sn]['composition_name'],fontsize=fsbase+4)

		#---annotations after the last bar
		rtx = out_specs['xvals'][-1][-1]+out_specs['width']
		cat_tops = np.array([0]+list(np.cumsum(yvals[-1][-1])[1::2]))
		midys = (cat_tops[1:]-cat_tops[:-1])/2.0+cat_tops[:-1]
		top_ys = np.cumsum(yvals[-1][-1])[1::2]
		el = mpl.patches.Ellipse((2, -1), 0.5, 0.5)
		for ncat,y in enumerate(top_ys):
			color = '#bdbdbd'
			ann = ax.annotate('ions bound\nto %d lipid%s'%(ncat+1,{0:''}.get(ncat,'s')),
				xy=(rtx,y),xycoords='data',
				xytext=(35,0),textcoords='offset points',
				size=16,va="center",
				bbox=dict(boxstyle="round",fc=color,ec="none"),
				arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					fc=color,ec="none",patchA=None,patchB=el,relpos=(0.2, 0.5),))
			patches.append(ann)

		#---spine formatting
		for ax in axes:
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)
			#ax.spines['bottom'].set_visible(False)
			ax.yaxis.set_ticks_position('left')
			ax.xaxis.set_ticks_position('none')

		for ax in axes:
			ax.set_yticklabels(['%d'%i for i in ax.get_yticks()],fontsize=fsbase+2)
		axes[0].set_ylabel('bound ions (within $\mathrm{\mathbf{%.1f} \AA}$)'%zonecut,fontsize=fsbase+4)

		#---custom legend on the second plot
		legend_labels,legend_patches = [],[]
		ion_names = ['K','NA','MG','Cal']	
		marks_sns = [sn for sn in [m[0] for m in [[sn for sn in sns if work.meta[sn]['cation']==i] 
			for i in ion_names] if m]]
		for nn,name in enumerate(ion_names):
			legend_labels.append(work.meta[marks_sns[nn]]['ion_label'])
			legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,
				facecolor=colors_ions[name],edgecolor='k',lw=2))
		legend_labels.append('bound to\nat least\n1 PtdIns')
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,color='gray'))
		legend_labels.append('not bound\nto PtdIns')
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,fc='w',edgecolor='gray',hatch='///'))
		legend = ax.legend(legend_patches,legend_labels,loc='upper left',fontsize=art['fs']['legend']+6,
			ncol=3,labelspacing=1.2,handleheight=2.0,markerscale=0.5,shadow=True,fancybox=True)

		#---render
		picturesave(picname,work.plotdir,backup=False,extras=patches,
			version=True,meta={'zonecut':zonecut,'sns':sns,'print_version':print_version},
			#---pad inches because the tagboxes get cut off, pdf because hatches must be thicker
			pad_inches=0.5,pdf=True)

	#---make the gif
	if False:
		os.remove(work.plotdir+'/fig.ion_binding2.gif')
		os.system('convert -delay 20 -loop 0 %s %s'%(
			work.plotdir+'/fig.ion_binding2.cutoff*v1.png',
			work.plotdir+'/fig.ion_binding2.gif'))
		
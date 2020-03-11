#!/usr/bin/env python

import numpy as np
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib import patches
import brewer2mpl

@autoload(plotrun)
def load():
	sns = work.sns()
	data = plotload('ion_binding')
	clrs = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors
	# settings
	max_by_ion = True
	onemax = 700
	bars = False
	art = {}

def plot_ion_binding_waterfall(ax,zonecut,sns_sep,do_single=True,
	multiplexed_cutoffs=None,charge_scaled=False,do_legend=True,extras=None,norm_leaflet=False,
	do_tag_bridges=True,ion_names=None,no_vert_labels=False,ymax_adjust=True,
	labeler=None,do_time_ticks=False,ptdins_name='PtdIns'):
	"""
	Single-panel plot of the cool ion-binding waterfall visualization.
	"""
	def set_time_pts(sn):
		global time_pts
		# handle variable number of timepoints here
		if time_pts_abs!=-1: time_pts = time_pts_abs
		else: 
			#! no: time_pts = int(len(data.this[sn]['0.wcs'])/time_pts_interval)
			duration = data.extras[sn]['end']-data.extras[sn]['start']
			time_pts = int(int(data.extras[sn]['end']-data.extras[sn]['start'])/time_pts_interval)
		return time_pts
	if extras is None: extras = []
	# format count data up to the maximum number of nearest-neighbors (hardcoded to 3)
	ncats = 3
	ymax_all = 0.
	global postdat_sep
	postdat_sep = [{} for s in sns_sep] 
	for pnum,sns in enumerate(sns_sep):
		postdat = dict([(sn,{}) for sn in sns])
		for sn in sns:
			try: data.set(name='ion_binding_combinator')
			except:
				if not multiplexed_cutoffs: 
					raise Exception('requested no multiplex cutoffs but we have too much data')
				else:
					# MULTIPLEXED CUTOFFS HERE
					# from the draft: $3.0, 3.0, 2.3, 2.6 \,\AA$ for K+, Na+, Mg2+, and Ca2+,
					zc = multiplexed_cutoffs[work.meta[sn]['cation']]
					data.set(name='ion_binding_combinator',select={('zonecut',):zc})
			set_time_pts(sn)
			if charge_scaled: factor = {'NA':1.,'Na,Cal':1.,'MG':2.,'K':1.,'Cal':2.}[work.meta[sn]['cation']]
			else: factor = 1.0
			if norm_leaflet and work.meta[sn]['composition_name'] == 'symmetric': factor = factor/2.0
			else: pass
			datb = data.this[sn]
			total = np.sum([datb['%d.wcs'%i].sum(axis=1) for i in range(ncats)],axis=0)*factor
			postdat[sn]['total'] = total
			sel_frames = np.linspace(0,total.shape[0]-1,time_pts).astype(int)
			postdat[sn]['total_few'] = [total[i]*factor for i in sel_frames]
			postdat[sn]['nn_count_few'] = [list(np.array([datb['%d.wcs'%i].sum(axis=1) 
				for i in range(ncats)]).T[j]*factor) for j in sel_frames]
			involves_ptdins = [np.sum(datb[i2s2(nn,'wcs')][:,np.array([ii for ii,i in enumerate(
				datb[i2s2(nn,'combonames')]) if work.meta[sn]['ptdins_resname'] in i])],axis=1)*factor
				for nn in range(ncats)]
			postdat[sn]['nn_count_few_ptdins'] = [[involves_ptdins[nn][j]
				for nn in range(ncats)] for j in sel_frames]
			postdat[sn]['nn_count_few_not_ptdins'] = [[
				postdat[sn]['nn_count_few'][j][i]-postdat[sn]['nn_count_few_ptdins'][j][i] 
				for i in range(ncats)] for j in range(len(sel_frames))]
		postdat_sep[pnum] = postdat
	global postdat # debugging
	postdat = postdat_sep[pnum]
	for pnum,sns in enumerate(sns_sep):
		# plot the total ions at a few time points
		barspecs_all = dict(lw=0)
		# if we just set yvals to the totals, we get a standard bar plot
		yvals = [postdat[sn]['total_few'] for sn in sns]
		# stacked bar plot from the totals by number of neighbors
		yvals = [postdat[sn]['nn_count_few'] for sn in sns]
		# specify by pip2 representation in the nearest neighbors
		yvals = []
		for ss,sn in enumerate(sns):
			set_time_pts(sn)
			ysn = []
			for t in range(time_pts):
				ysnt = []
				for l in range(ncats):
					for key in ['nn_count_few_not_ptdins','nn_count_few_ptdins']:
						ysnt.append(postdat[sn][key][t][l])
				ysn.append(ysnt)
			yvals.append(ysn)
		# bars are colored by ion
		barspecs_base = [{'color':colors_ions[work.meta[sn]['cation']]} for sn in sns]
		alpha_sweep = dict([(ii,i) for ii,i in enumerate(np.linspace(0.35,1,3))])
		# apply different patterns to each of the time_pts by adding e.g. alpha_sweep to the dict below
		barspecs = [[dict(**barspecs_base[ss]) for j in range(set_time_pts(sn))] 
			for ss,sn in enumerate(sns)]
		color_blender = lambda c,alpha : tuple(list(mpl.colors.ColorConverter().to_rgb(c))+[alpha]) 
		# add additional granular settings for stacked bar plots
		for ii in range(len(barspecs_base)): barspecs_base[ii].pop('color')
		# tune the stack so the plot is comprehensible
		barspecs_stack = [[[dict(
			# drop the alpha on the hatches because now I cannot change the hatch color and it's black
			alpha=alpha_sweep[k/2] if not k%2==0 else alpha_sweep[k/2]*0.65,
			#! lw=2 if k%2==1 else 0,
			lw=0,
			#! every time I remake this plot there is new hatch behavior
			#! hatch='--' if k%2==1 else None,
			#! edgecolor=color_blender(colors_ions[work.meta[sn]['cation']],alpha_sweep[k/2]),
			color='w' if k%2==0 else colors_ions[work.meta[sn]['cation']],
			#! edgecolor='k' if k%2==0 else colors_ions[work.meta[sn]['cation']],
			#! modifications for changing hatch handling in mpl
			edgecolor='k' if k%2==0 else color_blender(
				colors_ions[work.meta[sn]['cation']],alpha_sweep[k/2]),
			hatch='///' if k%2==0 else None,
			**barspecs_base[ss]) 
			for k in range(ncats*2)] 
			for j in range(set_time_pts(sn))] for ss,sn in enumerate(sns)]
		global out_specs
		out_specs = barmaker(ax,yvals,gap=1.0,
			barspecs_all=barspecs_all,barspecs=barspecs,barspecs_stack=barspecs_stack)
		ymax = out_specs['ymax']
		ymax_all = max([ymax,ymax_all])
		# if we are making a batch we fix the ymin and ymax
		if not do_single: ymax = [400,250][pnum]
		if ymax_adjust: ax.set_ylim(ymin,ymax*1.5) # changed 1.2 to 1.3 for BPS 
		# copied from elsewhere
		tagbox_ion = dict(facecolor='w',alpha=1.0,boxstyle="round,pad=0.3")
		tagbox_ptdins = dict(facecolor='w',lw=0,alpha=0.0,boxstyle="round,pad=0.4")
		tagbox_ptdins = tagbox_ion
		if do_time_ticks==None: ax.set_xticks([])
		elif do_time_ticks=='physiological_long':
			ax.set_xticks([0,50,100,151,201])
			ax.set_xticklabels([0,250,500,250,500])

		label_drop = (ymax-ymin)*0.2
		width = out_specs['width']
		xcenters = [np.mean(np.array(i)+out_specs['width']/2.0) for i in out_specs['xvals']]
		#! no_vert_labels is for the toc entry, a combination plot
		for snum,sn in enumerate(sns):
			text = work.meta[sn]['ion_label']
			#downdrop = -0.02 if not no_vert_labels else -0.04
			#tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*downdrop-label_drop,text,
			#	bbox=tagbox_ion,ha='center',va='top',rotation=-90 if not no_vert_labels and len(text)>25 else 0,
			#	color='k',fontsize=art['fs']['tags'])
			#extras.append(tb)
			#text = work.meta[sn]['ptdins_label']
			#tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*0.02-label_drop,text,
			#	bbox=tagbox_ptdins,rotation=-90 if not no_vert_labels and len(text)>15 else 0,ha="center",va="bottom",
			#	color='k',fontsize=art['fs']['tags'])
			#extras.append(tb)
			if labeler==None:
				#! hacking to move labels on the symmetric. tuned it manually
				if 'membrane-v509' in sns: extra_drop = 20.0
				else: extra_drop = 0.0
				print('extra drop is %s'%str(extra_drop))
				## HACKING FOR BPS
				# downdrop = -0.02 if not no_vert_labels else -0.04
				downdrop = 0.
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*downdrop-label_drop-extra_drop,text,
					bbox=tagbox_ion,ha='center',va='top',rotation=-90,
					color='k',fontsize=art['fs']['tags'])
				extras.append(tb)
				text = work.meta[sn]['ptdins_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*downdrop-label_drop/4.,text,
					bbox=tagbox_ptdins,rotation=-90,ha='center',va='top',
					color='k',fontsize=art['fs']['tags'])
				extras.append(tb)
				if do_time_ticks!=None: ax.set_xlabel('time',fontsize=16)
			# added extra label handling here. acquiring technical debt but things are too urgent!
			elif labeler=='physiological_long':
				label_drop = (ymax-ymin)*0.2
				downdrop = 0.
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*downdrop-label_drop,text,
					bbox=tagbox_ion,ha='center',va='top',rotation=0,
					color='k',fontsize=art['fs']['tags']+4)
				extras.append(tb)
				text = work.meta[sn]['ptdins_label']
				tb = ax.text(width*(xcenters[snum]),(ymax-ymin)*downdrop-label_drop/2.,text,
					bbox=tagbox_ptdins,rotation=-0,ha='center',va='top',
					color='k',fontsize=art['fs']['tags']+4)
				extras.append(tb)
				ax.set_xlabel('time',fontsize=16)
			else: raise Exception	

	# annotations after the last bar
	rtx = out_specs['xvals'][-1][-1]+out_specs['width']
	cat_tops = np.array([0]+list(np.cumsum(yvals[-1][-1])[1::2]))
	midys = (cat_tops[1:]-cat_tops[:-1])/2.0+cat_tops[:-1]
	top_ys = np.cumsum(yvals[-1][-1])[1::2]
	el = mpl.patches.Ellipse((2, -1), 0.5, 0.5)
	for ncat,y in enumerate(top_ys):
		color = '#bdbdbd'
		tag_bridge = 'ions bound\nto %d lipid%s'%(ncat+1,{0:''}.get(ncat,'s'))
		tag_bridge = '%d-bridge'%(ncat+1)
		if do_tag_bridges:
			ann = ax.annotate(tag_bridge,
				xy=(rtx,y),xycoords='data',
				xytext=(35,0),textcoords='offset points',
				size=16,va="center",
				bbox=dict(boxstyle="round",fc=color,ec="none"),
				arrowprops=dict(arrowstyle="wedge,tail_width=1.",
					fc=color,ec="none",patchA=None,patchB=el,relpos=(0.2, 0.5),))
			extras.append(ann)
	# spine formatting
	for ax in [ax]:
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		# ax.spines['bottom'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('none')
	ax.tick_params(axis='both',which='major',labelsize=art['fs']['axlabel'])
	if charge_scaled: ylabel = r'bound ions $\times$ charge ($e$)'
	else: ylabel = 'bound ions'
	if multiplexed_cutoffs:
		ax.set_ylabel(r'%s'%ylabel,fontsize=art['fs']['axlabel'])
	else: ax.set_ylabel(r'%s ($\vec{r}\leq\mathrm{\mathbf{%.1f} \AA}$)'%(ylabel,zonecut),
		fontsize=art['fs']['axlabel'])
	# custom legend on the second plot
	if do_legend:
		legend_labels,legend_patches = [],[]
		if not ion_names: ion_names = ['K','NA','MG','Cal']
		marks_sns = [sn for sn in [m[0] for m in [[sn for sn in sns 
			if work.meta[sn]['cation']==i and work.meta[sn].get('cations',i)==i] for i in ion_names] if m]]
		for nn,name in enumerate(ion_names):
			legend_labels.append(work.meta[marks_sns[nn]]['ion_label'])
			legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,
				facecolor=colors_ions[name],edgecolor=colors_ions[name],lw=2))
		legend_labels.append('bound to\nat least\n1 %s'%ptdins_name)
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,color='gray'))
		legend_labels.append('not bound\nto %s'%ptdins_name)
		legend_patches.append(mpl.patches.Rectangle((0,0),1.0,1.0,fc='w',edgecolor='gray',hatch='///'))
		if multiplexed_cutoffs:
			# mild offset for the legend so it clears the tops
			legend = ax.legend(legend_patches,legend_labels,loc='upper left',fontsize=art['fs']['legend'],
			ncol=3,labelspacing=1.2,handleheight=2.0,markerscale=0.5,
			bbox_to_anchor=(0.0,1.15 if not charge_scaled else 0.95))
		else:
			legend = ax.legend(legend_patches,legend_labels,loc='upper left',fontsize=art['fs']['legend'],
				ncol=2,labelspacing=1.2,handleheight=2.0,markerscale=0.5) 
				#! shadow=True,fancybox=True)
		frame = legend.get_frame()
		frame.set_edgecolor('w')
		frame.set_facecolor('w')
	return dict(ymax=ymax_all)

if __name__=='__main__':

	# several layouts (iterate manually for now)
	layout = ['comprehensive','physiological_with_EDPs','simple_snapshot',
		'physiological_long'][-1]

	# hatching shennanigans in various matplotlib versions
	mpl.rcParams['hatch.linewidth'] = 3.0
	from plotter.extras import i2s2

	# set multiplexed cutoffs to None to use zonecut instead
	multiplexed_cutoffs = {'NA':3.0,'K':3.0,'MG':2.3,'Cal':2.6,'Na,Cal':3.0}
	zonecut = 2.2

	charge_scaled = True
	norm_leaflet = True
	print_version = 'v20180426'
	colors_ions = {'NA':'green','Na,Cal':'green','MG':'red','Cal':'blue','K':'gray',}
	descriptor = '.cutoff%.1f'%zonecut if not multiplexed_cutoffs else ''
	# binding5 was the latest "comprehensive" plot
	picname = 'fig.ion_binding5%s%s'%(descriptor,'' if not multiplexed_cutoffs else '.special')
	# binding6 below is the new physiological-only plot
	picname = 'fig.ion_binding6%s%s'%(descriptor,'' if not multiplexed_cutoffs else '.special')
	if layout=='simple_snapshot': picname = 'fig.ion_binding_toc'
	if layout=='comprehensive': picname = 'fig.ion_binding7'

	# aesthetics
	global art
	time_pts_abs = -1 # previously called time_pts, set to 16, now variable
	time_pts_interval = 5000. # 5 ns per bar
	ymin = 0

	# selections
	if layout=='comprehensive':
		art = {'fs':{'legend':14,'title':20,'tags':14,'axlabel':14},}
		sns_sep = [work.vars['orders']['canon']['symmetric'],work.vars['orders']['canon']['asymmetric']]
		titles = ['symmetric','physiological']
		sns = sns_sep[0]+sns_sep[1]
		do_time_ticks = None
	elif layout=='physiological_with_EDPs':
		art = {'fs':{'legend':14,'title':20,'tags':14,'axlabel':16},}
		sns_sep = [work.vars['orders']['canon']['asymmetric_no_dilute']]
		sns = sns_sep[0]
	elif layout=='simple_snapshot':
		art = {'fs':{'legend':14,'title':20,'tags':16,'axlabel':16},}
		sns_sep = [['membrane-v531','membrane-v532']]
		sns = sns_sep[0]
	elif layout=='physiological_long':
		art = {'fs':{'legend':14,'title':20,'tags':16,'axlabel':16},}
		sns_sep = [['membrane-v563','membrane-v565']]
		sns = sns_sep[0]
		do_time_ticks = 'physiological_long'
	else: raise Exception

	extras = []
	if layout=='comprehensive':
		figsize = (20,10)
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,len(sns_sep)],'wratios':[len(s) for s in sns_sep],
			'wspace':0.1},'ins':[{'grid':[1,1]},{'grid':[1,1]}]})
		axes = [i for j in axes for i in j]
		ymax = 0.
	elif layout=='physiological_with_EDPs':
		figsize = (16,14)
		# dial in the wratio below to balance things
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,1+len(sns_sep)],'wratios':[4]+[len(s) for s in sns_sep],
			'wspace':0.3},'ins':[{'grid':[2,1]},{'grid':[2,1],'hratios':[7,1]}]})
	elif layout=='simple_snapshot':
		figsize = (16,10)
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,1+len(sns_sep)],'wratios':[3,1],'hspace':0.0,
			'wspace':0.1},'ins':[{'grid':[1,1]},{'grid':[4,1],'hratios':[2,5,1,1]}]})
		ymax = 0.
	elif layout=='physiological_long':
		figsize = (12,10)
		axes,fig = panelplot(figsize=figsize,
			layout={'out':{'grid':[1,1],'hspace':0.0,
			'wspace':0.1},'ins':[{'grid':[1,1]}]})
		ymax = 0.
	else: raise Exception

	# main plot loop
	if layout=='comprehensive':
		for anum,sns_this in enumerate(sns_sep):
			ax = axes[anum]
			detail = plot_ion_binding_waterfall(ax,sns_sep=[sns_this],zonecut=zonecut,extras=extras,
				multiplexed_cutoffs=multiplexed_cutoffs,charge_scaled=charge_scaled,do_legend=anum==1,
				norm_leaflet=norm_leaflet,do_tag_bridges=anum==1,labeler=None,do_time_ticks=do_time_ticks)
			print([len(i) for i in out_specs['xvals']])
			ymax = max(detail['ymax'],ymax)
			ax.set_title(titles[anum],fontsize=art['fs']['title'])
			# get the average totals for comparison between compositions
			if layout=='comprehensive':
				for sn in ['membrane-v511','membrane-v532']:
					if sn in postdat_sep[0]:
						status('%s totals average: %.2f'%(sn,postdat_sep[0][sn]['total'].mean()))
				"""
				RESULTS
				[STATUS] membrane-v511 totals average: 364.92
				[STATUS] membrane-v532 totals average: 435.48
				"""
	elif layout=='physiological_long':
		ax = axes[0]
		plot_ion_binding_waterfall(ax,sns_sep=sns_sep,zonecut=zonecut,
			multiplexed_cutoffs=multiplexed_cutoffs,charge_scaled=charge_scaled,extras=extras,
			ion_names=['MG','Cal'],do_legend=True,no_vert_labels=True,ymax_adjust=False,
			labeler='physiological_long',do_time_ticks=do_time_ticks,
			ptdins_name=r"$\mathrm{{PIP}_2}$")
		ax.tick_params(axis=u'both',which=u'both',length=0)
	elif layout=='physiological_with_EDPs':
		ax = axes[1][0]
		fig.delaxes(axes[1][1])
		plot_ion_binding_waterfall(ax,sns_sep=sns_sep,zonecut=zonecut,
			multiplexed_cutoffs=multiplexed_cutoffs,charge_scaled=charge_scaled,extras=extras)

		# add EDPs (depends strongly on ptdins_electron_density_profiles.py)
		#! weird import
		plotrun.routine = []
		outgoing = dict(autoload=autoload,plotrun=plotrun,autoplot=autoplot)
		#!!! why are there three copies of the EDPs data?
		execfile('calcs/ptdins_electron_density_profiles.py',globals(),outgoing)
		post_process_edps = outgoing['post_process_edps']
		kwargs = dict(sns=sns,data=data,bilayer_width=outgoing['bilayer_width'])
		if multiplexed_cutoffs:
			#! redundant data for some reason! needs fixed in omnicalc
			reds = [ii for ii,i in enumerate(data.names) if i[0]=='electron_density_profiles']
			data.this = data.data[reds[0]]
		else:
			data.set(name='electron_density_profiles')
		# get the post for EDP and export it to globals
		outgoing['post'] = post = post_process_edps(**kwargs)
		plot_interesting = outgoing['plot_interesting']
		# if you just run plot_interesting here, you would get the original EDPs from that code
		def get_panel(sn,ion_name): return {'membrane-v531':0,'membrane-v532':1}[sn]
		def get_color(sn,ion_name):
			if ion_name=='anion': return 'k'
			else: return colorize(work.meta[sn],comparison='protonation')
		sns_this = ['membrane-v531','membrane-v532']
		extras.extend(plot_interesting(axes=axes[0],get_panel=get_panel,
			do_external=True,sns_this=sns_this,layout_style='by_sn',
			xlim=(1.5,4.5),mark_peaks=False,get_color=get_color,
			nsegs=5,early_color='#A9A9A9',fs_labels=14,fs_legend=14,show_title=False))

	elif layout=='simple_snapshot':
		ax = axes[1][1]
		for j in [0,2,3]: fig.delaxes(axes[1][j])
		plot_ion_binding_waterfall(ax,sns_sep=sns_sep,zonecut=zonecut,
			multiplexed_cutoffs=multiplexed_cutoffs,charge_scaled=charge_scaled,extras=extras,
			ion_names=['MG','Cal'],do_legend=False,no_vert_labels=True,
			ptdins_name=r"$\mathrm{{PIP}_2}$")
		ax.tick_params(axis=u'both',which=u'both',length=0)
		# hijacked from lipid_lipid_bonds
		ax = axes[0][0]
		import matplotlib.image as mpimg
		# previously sent a single image to sidestack but the following routine makes it square
		image = mpimg.imread(os.path.join(work.plotdir,'snap.vmd.bilayer.png'))
		border = get_blank_border(image)
		sizes = np.array([np.ptp(i) for i in border])
		padding = np.array([(max(sizes)-s)/2. for s in sizes]).astype(int)
		zooms = np.array([[border[i][0]-1*padding[i],border[i][1]+padding[i]] 
			for i in range(2)])
		# note that you have to reverse the order of dimensions here
		image_zoomed = image[zooms[1][0]:zooms[1][1],zooms[0][0]:zooms[0][1]]
		ax.imshow(image_zoomed)
		ax.axis('off')

	else: raise Exception
	# note that the weird plateau is legitimate; the plot is not cut off for v565
	# leaflet norming also matches the axes for a direct comparison between compositions
	if norm_leaflet and layout=='comprehensive':
		for ax in axes: ax.set_ylim((0,ymax))
	# PDF method previously used to make hatches thicker
	meta = {'zonecut':zonecut if not multiplexed_cutoffs else multiplexed_cutoffs,
		'charge_scaled':charge_scaled,'sns':sns,'print_version':print_version,'time_pts':time_pts_abs},
	fn = picturesave(picname,work.plotdir,backup=False,extras=extras,dpi=500,
		version=True,meta=meta,pad_inches=0.5,pdf=False,form='pdf')
	#! abandoned
	if False:
		# trim white space
		#! this function is useful. formalize it and put it somewhere
		image = mpimg.imread(fn)
		border = get_blank_border(image)
		sizes = np.array([np.ptp(i) for i in border])
		padding = np.array([(max(sizes)-s)/2. for s in sizes]).astype(int)
		zooms = np.array([[border[i][0]-1*padding[i],border[i][1]+padding[i]] 
			for i in range(2)])
		# note that you have to reverse the order of dimensions here
		image_zoomed = image[zooms[1][0]:zooms[1][1],zooms[0][0]:zooms[0][1]]
		fig = plt.figure(figsize=figsize)
		ax = plt.subplot(111)
		ax.imshow(image)
		picturesave(picname,work.plotdir,backup=False,extras=extras,
			version=True,meta=meta,pad_inches=0.5,pdf=False)
	conv_cmd = '-trim -bordercolor white -border 25 +repage'
	os.system('convert %s %s %s'%(fn,conv_cmd,re.sub(r'\.pdf$','.trim.pdf',fn)))

"""
# checking the number of bars
# [len(i) for i in out_specs['xvals']]
[16, 16, 16, 9, 16, 16]
[STATUS] membrane-v511 totals average: 364.92
[16, 16, 16, 9, 16, 16, 100, 16, 16, 16, 16, 100, 16]
[STATUS] membrane-v532 totals average: 435.48
[STORE] searching pictures
[STORE] saving picture to fig.ion_binding7.v1.pdf
[WARNING] you are saving as pdf and only png allows metadata-versioned pictures

NOTE THAT hatching looked different in firefox vs a rendered latex document even for the same file! no idea why!
"""


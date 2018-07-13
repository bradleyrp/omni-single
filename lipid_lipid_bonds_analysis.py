#!/usr/bin/env python

"""
LIPID-LIPID hydrogen bond and salt-bridge analysis
"""

import itertools

### LOAD

def load_actinlink(data):
	"""
	Plot settings for the "actinlink" project. Called by load() and sent to globals.
	"""
	sns = actinlink_sns_mdia2
	color_by_simulation = actinlink_color_by_simulation
	replicate_mapping = actinlink_replicate_mapping
	extra_labels = actinlink_extra_labels
	if set(sns)!=set(work.sns()): raise Exception
	nmol_counts = dict([(sn,dict(zip(*[data['hydrogen_bonding'][sn]['data'][k] 
		for k in ['resnames','nmols']]))) for sn in sns])
	residue_renamer = None
	outgoing = dict(sns=sns,color_by_simulation=color_by_simulation,
		replicate_mapping=replicate_mapping,extra_labels=extra_labels,nmol_counts=nmol_counts,
		residue_renamer=residue_renamer)
	return outgoing

def load_ptdins(data):
	"""
	Plot settings for the "ptdins" project. Called by load() and sent to globals.
	"""
	collection = 'asymmetric_all_no_mixed'
	sns = work.metadata.collections[collection]
	if set(sns)!=set(work.sns()): raise Exception('plot is set up for %s'%collection)
	def color_by_simulation(sn):
		return colorize(work.meta[sn],comparison='asymmetric_all')
	replicate_mapping = dict([(sn,[sn]) for sn in sns])
	extra_labels = dict([(sn,'%s, %s'%(work.meta.get(sn,{}).get('ptdins_label','ptdins'),
		work.meta.get(sn,{}).get('ion_label','ion'))) 
		for sn in sns])
	nmol_counts = dict([(sn,dict(zip(*[data['hydrogen_bonding'][sn]['data'][k] 
		for k in ['resnames','nmols']]))) for sn in sns])
	def residue_renamer(resname):
		if resname in ['PI2P','P35P','PIPP','PIPU','SAPI']: return 'PtdIns'
		else: return resname
	outgoing = dict(sns=sns,color_by_simulation=color_by_simulation,
		replicate_mapping=replicate_mapping,extra_labels=extra_labels,nmol_counts=nmol_counts,
		residue_renamer=residue_renamer)
	return outgoing

@autoload(plotrun)
def load():
	"""
	Load once and export to globals. Uses project-specific load-functions.
	"""
	data,calc = work.plotload('lipid_lipid_bonds_analysis')
	# format of the bonds data
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']
	# load switch by project name
	project_name = os.path.basename(os.getcwd())
	# write to locals so the Observer can export to globals
	#! note that previous uses of the Observer may have forgotten to get _locals from the function
	_locals = globals()['load_%s'%project_name](data=data)

### PLOT

def make_legend(ax,replicate_mapping_this,ncol=1,keys=None):
	"""Make a legend. Note that this is redundant with plot-lipid_rdfs.py."""
	legendspec = []
	for sn_general,sns in replicate_mapping_this:
		if not keys or sn_general in keys:
			legendspec.append(dict(
				#! absurd one-liners below
				name=extra_labels.get(sn_general,work.meta.get(sn_general,{}).get('label','')),
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.),ncol=ncol)
	frame = legend.get_frame()
	frame.set_edgecolor('k')
	frame.set_facecolor('w')
	return legend

def plot_bonds(name,kinds,**kwargs):
	"""
	Plot the distributions of observed bond counts.
	"""
	merged = kwargs.get('merged',False)
	normed = kwargs.get('normed',False)
	symmetrize = kwargs.get('symmetrize',False)
	# fancy imshow method is a great substitute for error bars
	style = kwargs.get('style',['bars','imshow'][1])
	# count scaling applies to bars and imshow counts (not normalized scores)
	count_scaling = kwargs.get('count_scaling',['absolute','global','local'][1])
	count_scaling_norm = kwargs.get('count_scaling_norm',['global','local'][0])
	# aesthetics
	wspace = kwargs.get('wspace',0.4)
	hspace = kwargs.get('hspace',0.2)
	legend_ncols = kwargs.get('legend_ncols',1)
	#! ignoring POPC for now
	#! farm this out to the specs
	resnames_exclude = []#['POPC']
	debug_black_background = kwargs.get('debug_black_background',False)
	def get_norm_factor(sn,combo):
		if normed:
			if 'PtdIns' in combo:
				combo = tuple([c if c!='PtdIns' else work.meta[sn]['ptdins_resname'] for c in combo])
			nmols = [nmol_counts.get(sn,nmol_counts[dict(replicate_mapping).get(sn,[sn])[0]])[r] 
				for r in combo]
			norm_factor = nmols[0]*nmols[1]
		else: norm_factor = 1.0
		return norm_factor
	global post,post_norm
	post = {}
	global replicate_mapping_this
	if merged: replicate_mapping_this = replicate_mapping 
	else: replicate_mapping_this = [(sn,[sn]) for sn in sns]
	# tabulate nframes by simulation
	nframes_by_sn = dict([(sn,len(data['salt_bridges'][sn]['data']['observations'])) for sn in sns])
	valid_frames_by_kind = dict([(sn,dict([(k,data[k][sn]['data']['valid_frames']) for k in kinds])) 
		for sn in sns])
	# tabulate the bonds
	for sn_group,sns_this in replicate_mapping_this:
		for snum,sn in enumerate(sns_this):
			for knum,kind in enumerate(kinds):
				# collect bonds
				bonds,obs = [data[kind][sn]['data'][k] for k in ['bonds','observations']]
				# collect pairs
				if len(bonds)==0: continue
				subjects = np.unique(bonds[:,rowspec.index('subject_resname')])
				targets = np.unique(bonds[:,rowspec.index('target_resname')])
				resnames = [i for i in np.unique(np.concatenate((subjects,targets)))
					if i not in resnames_exclude]
				combos = list(itertools.product(resnames,repeat=2))
				if residue_renamer!=None: 
					combos_renamed = [tuple([residue_renamer(i) for i in j]) for j in combos]
				else: combos_renamed = combos
				for combo,combo_renamed in zip(combos,combos_renamed):
					rows = np.where(np.all((
						bonds[:,rowspec.index('subject_resname')]==combo[0],
						bonds[:,rowspec.index('target_resname')]==combo[1],
						bonds[:,rowspec.index('subject_resid')]!=bonds[:,rowspec.index('target_resid')],
						),axis=0))[0]
					if len(rows)>0:
						if sn_group not in post: post[sn_group] = {}
						if combo_renamed not in post[sn_group]: 
							#! we have to load with zeros here because one simulation has a single CHL1-CHL1
							#! ... bond while the other replicate lacks it, and later we are summing and 
							#! ... concatenating hence we need proper placeholders
							post[sn_group][combo_renamed] = [[np.zeros(nframes_by_sn[ii]) 
								for ii in sns_this] for jj in kinds]
						counts = obs[:,rows].sum(axis=1)
						# we have to be careful about counting. we want to concatenate simulations
						# ... and sum different kinds of bonds. hence we construct a two-dimensional 
						# ... list here, and then process it when we are done. the first dimension gets the 
						# ... while the second dimension is contatenated 
						#! post[sn_group][combo] = np.append(post[sn_group][combo],counts)
						post[sn_group][combo_renamed][knum][snum] = counts
	# post-processing to sum bond types (first dimension) then concatenate replicates
	for sn_group in post:
		for combo in post[sn_group]:
			try:
				post[sn_group][combo] = np.concatenate(np.array(post[sn_group][combo]).sum(axis=0))
			# quarantining the hacking for PtdIns
			except: 
				#! dark is stuck on scipy 0.19.1 while the actinlink machine has 1.0.0 where they fixed the
				#! ... bug with kdtree where you needed two copies of the box vectors and some frames broke
				#! ... however rather than recompute with new scipy I will just use the standard fix for when
				#! ... there is a mismatch in the frames: valid_frames tells us which ones work
				if len(dict(replicate_mapping_this)[sn_group])>1: 
					raise Exception('the valid_frames fix does not work for replicates')
				intersect = np.intersect1d(*[data[k][sn_group]['data']['valid_frames'] for k in kinds])
				valid_frame_map = [np.array([np.where(np.arange(nframes_by_sn[sn_group])==i)[0][0] 
					for i in intersect]) for k in kinds]
				#! running out of names. hacking through this swiftly!
				valid_frame_remap = [np.array([np.where(i==v)[0][0] 
					for i in intersect]) for v in valid_frame_map]
				try:
					post[sn_group][combo] = [[post[sn_group][combo][kk][ss][valid_frame_remap[kk]] 
						for ss,s in enumerate(dict(replicate_mapping_this)[sn_group])] 
						for kk,kind in enumerate(kinds)]
				# some bond types have an empty bond list hence an empty array so we hack through this
				# ... by inferring the length of the trajectory and using zeros when the array is empty. 
				# ... note that this only happens when you drop the salt bridge cutoff down from 3.4 to 2,2
				except:
					post[sn_group][combo] = [[i if len(i)>0 else [np.zeros(max([len(k[0]) 
						for k in post[sn_group][combo]]))] for i in j] for j in post[sn_group][combo]]					
				post[sn_group][combo] = np.concatenate(np.array(post[sn_group][combo]).sum(axis=0))
	# diagnostic
	try:
		status('mean counts v532, PtdIns-PtdIns: %s'%
			post['membrane-v532'][('PtdIns','PtdIns')].mean(),tag='diagnostic')
		status('mean counts v532, DOPE-DOPE: %s'%
			post['membrane-v532'][('DOPS','DOPE')].mean(),tag='diagnostic')
	except: pass
	if symmetrize:
		post_symmetric = {}
		for sn in post:
			post_symmetric[sn] = {}
			combos_sym = list(set([tuple(sorted(i)) for i in post[sn].keys()]))
			for combo in combos_sym:
				post_symmetric[sn][combo] = np.sum([post[sn][c] for c in post[sn] 
					if set(c)==set(combo)],axis=0)
		post = post_symmetric
	# reformulate post into images
	combos_u = list(set([k for j in [i.keys() for i in post.values()] for k in j]))
	sns_by_mapping = list(zip(*replicate_mapping_this))[0]
	try: max_count = int(np.concatenate([np.concatenate(post[sn].values()) for sn in post]).max())
	except: 
		try: max_count = int(np.concatenate([np.array(post[sn].values()).reshape(-1) for sn in post]).max())
		except: 
			max_count = int(np.array([i for j in post[sn].values() for i in j for sn in post]).max())
	if style=='imshow':
		images = dict([(combo,{}) for combo in combos_u])
		for combo in combos_u:
			sns_this = [sn for sn in sns_by_mapping if sn in post and combo in post[sn]]
			images[combo] = np.zeros((len(sns_this),max_count))
			for snum,sn in enumerate(sns_this):
				if combo not in post[sn]: continue
				# adding counts here for salt bridges and hydrogen bonds
				counts = post[sn][combo]
				# bins also get the normalization factor along with the counts when we normalize by comp
				if not normed: 
					hist,_ = np.histogram(counts,bins=(np.arange(max_count+1)-0.5),normed=True)
					images[combo][snum] = hist
		# each image of the lipid-count-normalized histograms has a maximum value which we get from
		# ... the maximum of each simulation's max bond counts divided by the normalization factor
		post_norm,post_norm_interp = {},{}
		for combo in combos_u:
			sns_this = [sn for sn in sns_by_mapping if sn in post and combo in post[sn]]
			obs_by_snum = []
			post_norm_interp[combo] = {}
			for sn in sns_this:
				counts = post[sn][combo]
				norm_factor = get_norm_factor(sn,combo)
				counts_normed = counts/norm_factor
				# diagnostic
				if normed and sn=='membrane-v532' and combo==('PtdIns','PtdIns'):
					status('normed mean counts v532, PtdIns-PtdIns: %s'%
						counts_normed.mean(),tag='diagnostic')
				if normed and sn=='membrane-v532' and combo==('DOPE','DOPE'):
					status('normed mean counts v532, DOPE-DOPE: %s'%
						counts_normed.mean(),tag='diagnostic')
				#! alternate count 
				vals_u,idx,counts_u = np.unique(counts_normed,return_counts=True,return_index=True)
				obs_by_snum.append(counts_normed)
				post_norm_interp[combo][sn] = dict(counts_normed=counts_normed,
					vals_u=vals_u,counts_u=counts_u,idx=idx,norm_factor=norm_factor,counts=counts)
			post_norm[combo] = obs_by_snum
		# loop over normed data and interpolate. note that iterpolation is *essential* to rendering
		# ... this data correctly because we have integer counts normalized by floats (see note above)
		images,images_extra = {},{}
		for combo in combos_u:
			sns_this = [sn for sn in sns_by_mapping if sn in post and combo in post[sn]]
			# the number of bins is the maximum of the discrete count of bins. we could modify this later to 
			# ... exclude some whitespace if necessary
			nbins = float(max_count)
			# if we are norming we use a standard resolution and then round into the bins
			if normed: base_res = 1000
			# if we are not norming we keep absolute counts
			else: base_res = int(nbins)
			# we construct a high resolution image and just round things onto it
			image = np.zeros((base_res,base_res,4))
			#! debugging below
			if debug_black_background:
				image[...,:3] = 0.0
				image[...,3] = 1.0
			shrink_max = float(min([v['norm_factor'] for k,v in post_norm_interp[combo].items()]))
			# after normalization, the top of the bar is the maximum count by the maximum shrink
			max_factor = nbins/shrink_max
			# assume we want a square tile so the width of the simulation with the maximum shrink should be
			# one over the number of simulations
			#! sns_this is wrong!
			width_base = 1./len(sns_this)
			# the base area is the area of one tile, corresponding to one bin in the simulation with the max
			# ... shrink, hence the one that runs all the way to the top, which will be larger than the others
			# ... because it has the smallest normalization factor.
			area_base = 1.0/nbins*width_base
			for sn in sns_this:
				# each simulation has a maximum value given by the shrink factor
				vals_u,counts_u,idx,norm_factor,counts_normed,counts = [post_norm_interp[combo][sn][k] 
					for k in ['vals_u','counts_u','idx','norm_factor','counts_normed','counts']]
				# this is the maximum for this simulation after the normalization factor
				post_norm_interp[combo][sn]['max_this_abs'] = max_this_abs = float(nbins)/norm_factor
				post_norm_interp[combo][sn]['max_this_rel'] = max_this_rel = max_this_abs/max_factor
				# bin heights are given in relative units
				post_norm_interp[combo][sn]['height_bins'] = height_bins = max_this_rel/nbins
				# bin widths depend on the base area to maintain constant area for each bin
				post_norm_interp[combo][sn]['width_bins'] = width_bins = area_base/height_bins
			# second pass to make the image is necessary to get the width sums
			total_width = sum([post_norm_interp[combo][s]['width_bins'] for s in post_norm_interp[combo]])
			xpos = 0
			images_extra[combo] = dict(xticks=[])
			for sn in sns_this:
				vals_u,counts_u,idx,norm_factor,counts_normed,counts = [post_norm_interp[combo][sn][k] 
					for k in ['vals_u','counts_u','idx','norm_factor','counts_normed','counts']]
				width_bins = post_norm_interp[combo][sn]['width_bins']/total_width
				height_bins = post_norm_interp[combo][sn]['height_bins']
				# at this point we have the relative height, the bin heights, and the unique counts
				# ... however the observations do not need to span the whole set so we have to slot them 
				# ... into place. we use the following trick: multiply by norm factor and you get ints.
				# ... note that rounding is essential otherwise you get repeats
				bins_which = np.round(vals_u*norm_factor).astype(int)
				# normalize the intensity by max instead of sum
				bins_intensity = counts_u.astype(float)/counts_u.max()
				for bin_num,intensity in zip(bins_which[bins_which>0],bins_intensity[bins_which>0]):
					# now we are ready to assemble the squares
					lrbt = np.array([xpos,xpos+width_bins,
						bin_num*height_bins,(bin_num+1)*height_bins])
					inds = np.round(lrbt*base_res).astype(int)
					# if we multiply the rgba value by intensity then it fades to grey
					# ... so instead we multiply only the alpha by intensity. this fades to background
					# ... and there is nothing behind it so no overlays. this looks better.
					image[inds[0]:inds[1],inds[2]:inds[3]] = np.array(mpl.colors.to_rgba(
						color_by_simulation(sn)))*np.array([1.,1.,1.,intensity])
				images_extra[combo]['xticks'].append((inds[1]-0.5)/total_width)
				xpos += width_bins
			# transpose so the image looks correct on imshow using the 3-dimensional color values
			images[combo] = np.transpose(image,(1,0,2))/image.max()
		max_row_global = max([max(np.where(np.sum(v[:,:,:3].sum(axis=2),axis=1)>0)[0]) 
			for v in images.values()])
	# plot
	collected_means = {}
	axes,fig = square_tiles(len(combos_u),figsize=16,wspace=wspace,hspace=hspace)
	for cnum,combo in enumerate(combos_u):
		collected_means[combo] = {}
		ax = axes[cnum]
		# always tabulate with bars to get the global max which is more convenient then the elaborate
		# ... imshow manipulations required to standardize the y-axes for the imshow version
		counter = 0
		sns_this = [sn for sn in sns_by_mapping if sn in post and combo in post[sn]]
		for sn in sns_this:
			collected_means[combo][sn] = []
			if combo in post[sn]:
				norm_factor = get_norm_factor(sn,combo)
				mean = post[sn][combo].mean()/norm_factor
				std = post[sn][combo].std()/norm_factor
				collected_means[combo][sn].append(mean+std)
				if style=='bars':
					ax.bar([counter],mean,width=1,color=color_by_simulation(sn))
					#! error bar is really huge hence the distribution is not normal. previously:
					ax.errorbar([counter],mean,std,xerr=0,
						zorder=3,fmt='none',ecolor='k',lw=0.5)
				counter += 1
		if style=='imshow':
			cmap_name = 'binary'
			if not normed:
				# when normed, this option is not extremely useful because it trims all images by their 
				# ... common lack of whitespcae, but since they have different maximum scores, then we cannot 
				# ... compare their scores directly. hence it does the opposite of its name!
				if count_scaling=='global': 
					raw = images[combo][:max_row_global]
					max_row_this = max_row_global
				# local also trims the most whitespace possible for each simulation but leaves you with 
				# ... different scaling on the y-axis when the data are normed by composition
				elif count_scaling=='local':
					# find the maximum row with a nonzero element (excluding opacity)
					# ... note that this line was written quickly using a sum trick over columns
					max_row_this = max(np.where(np.sum(image[:,:,:3].sum(axis=2),axis=1)>0)[0])
					raw = images[combo][:max_row_this]
				# when normed the following option creates different y-axis and a lot of whitespace because 
				# ... it assumes that each combo could have the maximum discrete count for all bonds, 
				# ... divided by the normalization factor for that specific combination (by composition) 
				# ... hence it is not useful (when normed)
				elif count_scaling=='absolute':
					raw = images[combo]
					max_row_this = base_res-1
				else: raise Exception
			# we ignore the other methods given above for norming because 
			else:
				max_row_this = max(np.where(np.sum(image[:,:,:3].sum(axis=2),axis=1)>0)[0])
				raw = images[combo][:max_row_this]
			#! interesting note: had to do the transpose when the image is made after switching to 4-d
			#! ... image values and then I had to change axis commands and array slicing and the shape 
			#! ... indexing in the aspect ratio in five places to correct this single upstream transpose\
			kwargs_to_imshow = {}
			max_y = max([post[sn][combo].max()/get_norm_factor(sn,combo) for sn in post if combo in post[sn]])
			if normed: kwargs_to_imshow['extent'] = [0.,1.,0.,max_y]
			ax.imshow(raw,origin='lower',interpolation='nearest',
				cmap=mpl.cm.__dict__[cmap_name],zorder=1,**kwargs_to_imshow)
		if style=='imshow' and not normed: ax.set_aspect(float(raw.shape[1])/raw.shape[0])
		ax.set_title('%s-%s'%tuple([work.vars.get('names',{}).get('short',{}).get(p,p) for p in combo]))
		if normed: ax.set_ylabel('score')
		else: ax.set_ylabel('bonds')
		#! have not added grid lines on normed yet
		if not normed and style=='imshow':
			ax.set_xticks(images_extra[combo]['xticks'])
			ax.xaxis.grid(True,which='major',zorder=0)
			ax.set_axisbelow(True)
		else: ax.set_xticks([])
		ax.set_xticklabels([])
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		if normed and count_scaling_norm=='local':
			y_max_this = max([max(i) for i in collected_means[combo].values()])
			ax.set_ylim((0,y_max_this))
			if style=='imshow':
				ax.set_aspect(1./float(y_max_this))
			else: ax.set_aspect('auto')
	# maximum from the bars is useful as a global maximum on the imshow
	max_y_bars = max([max([max(i) for i in m.values()]) for m in collected_means.values()])
	if style=='bars':
		if (normed and count_scaling_norm=='global') or (not normed and count_scaling=='global'):
			for ax in axes: ax.set_ylim((0,max_y_bars))
	#! normed-max is currently handled at the project level
	elif style=='imshow' and normed and count_scaling_norm=='global':
		for ax in axes:
			normed_y_max_global = max_y_bars
			ax.set_ylim((0,normed_y_max_global))
			ax.set_aspect('auto')
	extras = [make_legend(axes[-1],replicate_mapping_this=replicate_mapping_this,
		keys=[i for i,j in replicate_mapping_this],ncol=legend_ncols)]
	title_names = {'hydrogen_bonding':'hydrogen bonds','salt_bridges':'salt bridges'}
	extras.append(plt.suptitle(' and '.join(title_names[k] for k in kinds)+
		(' (normalized)' if normed else ''),fontsize=20,y=0.92))
	picturesave('fig.lipid_lipid_bonds.%s'%(name),work.plotdir,
		backup=False,version=True,meta={},extras=extras)
	post_debugger[tuple(kinds)] = post

def plots_actinlink():
	"""
	Plot actinlink bonds analysis including merged replicates.
	"""
	figspec = {}
	global post_debugger
	post_debugger = {}
	kind_map = {'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	symmetrize = True
	for normed in [True,False]:
		for merged in [True,False]:
			for key,kind in kind_map.items():
				for style in ['imshow','bars']:
					for count_scaling in ['global','local']:
						name = '%s.%s'%(style,key)+('.merged' if merged else '')+\
							('.normed' if normed else '')+('.local' if count_scaling=='local' else '')
						figspec[name] = {'merged':merged,'kinds':kind,'normed':normed,
							'symmetrize':symmetrize,'style':style,
							'count_scaling_norm':count_scaling,'count_scaling':count_scaling}
	test_key = ['imshow.hbonds_salt.normed.local',None][-1]
	if test_key: figspec = {test_key: figspec[test_key]}
	for key,spec in figspec.items(): plot_bonds(name=key,**spec)

def plots_ptdins():
	"""
	Plot actinlink bonds analysis including merged replicates.
	"""
	figspec = {}
	global post_debugger
	post_debugger = {}
	kind_map = {
		'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	merged = False
	symmetrize = True
	for normed in [True,False][:1]:
		for style in ['imshow','bars'][-1:]:
			for key,kind in kind_map.items():
				for count_scaling in ['global','local']:
					name = '%s.%s'%(style,key)+('.merged' if merged else '')+\
						('.normed' if normed else '')+('.local' if count_scaling=='local' else '')
					figspec[name] = {'merged':merged,'kinds':kind,'normed':normed,
						'symmetrize':symmetrize,'style':style,
						'count_scaling_norm':count_scaling,'count_scaling':count_scaling}
	test_key = ['imshow.hbonds_salt.normed.local',None][-1]
	if test_key: figspec = {test_key: figspec[test_key]}
	for key,spec in figspec.items(): plot_bonds(name=key,**spec)

@autoplot(plotrun)
def plots():
	"""Project-specific plotting sets."""
	project_name = os.path.basename(os.getcwd())
	globals()['plots_%s'%project_name]()	

def plots():
	"""
	Plot *many* versions of this visualization.
	Deprecated plot function (see custom versions above).
	"""
	kind_map = {
		'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	figspec = {
		('%s.%s'%(style,key)+
		{1:'.merged',0:''}[merged]+
		{1:'.normed',0:''}[normed]+
		{'local':'.local','global':''}[count_scaling]):{
			'merged':merged,'kinds':kind,'normed':normed,'symmetrize':symmetrize,
			'style':style,'count_scaling_norm':count_scaling,'count_scaling':count_scaling}
			for normed in [True,False]
			for merged in [True,False]
			for key,kind in kind_map.items()
			for style in ['imshow','bars']
			for count_scaling in ['global','local']
			for symmetrize in [True,False]}
	figspecs = {}
	test_key = ['imshow.hbonds.normed.local',None][-1]
	if test_key: figspec = {test_key: figspec[test_key]}
	for key,spec in figspec.items(): plot_bonds(name=key,**spec)

plotrun.routine = None
if __name__=='__main__': 

	sn = 'membrane-v532'
	kind = ['hydrogen_bonding','salt_bridges'][-1]
	bonds,obs = [data[kind][sn]['data'][k] for k in ['bonds','observations']]
	red_idx = obs.mean(axis=0).argsort()[::-1]
	counts_top = obs.mean(axis=0)[red_idx]
	red = bonds[red_idx][:,np.array([0,2,3,5])]
	_,idx = np.unique(red,return_index=True,axis=0)
	# ranked_bonds tells us the most common hbonds while counts_top gives the prevalence
	ranked_bonds = red[np.sort(idx)]
	# atom_ranks gives the prevalence ranking for a specific atom on the target lipid
	atom_ranks = {'acceptor':{},'donor':{}}
	# select bonds where one participant is a PIP2
	# the donor is the first set of columns on the bond list
	for name,i,j in [('donor',0,1),('acceptor',2,3)]:
		these_inds = np.where(ranked_bonds[:,0]=='PI2P')[0]
		for t in these_inds: 
			if ranked_bonds[t,j] not in atom_ranks[name]: 
				atom_ranks[name][ranked_bonds[t,j]] = counts_top[t]
			else: atom_ranks[name][ranked_bonds[t,j]] += counts_top[t]
	pass

	#!!! do not use the variable script
	this_script = {}
	scale_factor = {'hydrogen_bonding':0.2,'salt_bridges':0.1}[kind]
	for color,name in [(9,'donor'),(10,'acceptor')]:
		this_script[name] = ('\n'.join(["set i %(ind)d;set size %(size)f;mol addrep $molno;mol modcolor $i $molno ColorID %(color)d;mol modselect $i $molno name %(name)s;mol modstyle $i $molno Points $size;mol modstyle $i $molno Beads $size 12.0;mol modmaterial $i $molno EdgyGlass"%{'size':scale_factor*j,'ind':ii+2,'name':i,'color':color} for ii,(i,j) in enumerate(atom_ranks[name].items())]))

	this_script['header'] = """

color Display Background white
display projection Orthographic

set molno 0
mol new PI2P.gro
mol delrep 0 $molno

mol addrep $molno
mol modstyle 0 $molno DynamicBonds 1.700000 0.100000 6.000000
mol addrep $molno
mol modstyle 1 $molno Licorice 0.100000 10.000000 10.000000

"""
	this_script['middle'] = """

set molno 1
mol new PI2P.gro
mol delrep 0 $molno

mol addrep $molno
mol modstyle 0 $molno DynamicBonds 1.700000 0.100000 6.000000
mol addrep $molno
mol modstyle 1 $molno Licorice 0.100000 10.000000 10.000000

"""

	with open('/home/rpb/omicron/factory/calc/ptdins/aside-pip2-hbond-sites/run_%s.tcl'%kind,'w') as fp:
		for key in 'header donor middle acceptor'.split():
			fp.write(this_script[key]+'\n\n')

"""
ABSOLUTELY DO NOT USE THIS FILE UNLESS YOU ARE DOING THE STUFF ABOVE WITH VMD
SEE lipid_lipid_bonds.py instead!!!
"""

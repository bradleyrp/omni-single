#!/usr/bin/env python

"""
LIPID-LIPID hydrogen bond and salt-bridge analysis
"""

import itertools

def load_actinlink(data):
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
	collection = 'asymmetric_all'
	sns = work.metadata.collections[collection]
	if set(sns)!=set(work.sns()): raise Exception('plot is set up for %s'%collection)
	def color_by_simulation(sn):
		return colorize(work.meta[sn],comparison='asymmetric_all')
	replicate_mapping = dict([(sn,[sn]) for sn in sns])
	extra_labels = dict([(sn,'%s, %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label'])) 
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
	data,calc = work.plotload('lipid_lipid_bonds_analysis')
	# load switch by project name
	project_name = os.path.basename(os.getcwd())
	# load to locals then automatically to globals
	locals().update(**globals()['load_%s'%project_name](data=data))
	# format of the bonds data
	rowspec = ['subject_resname','subject_resid','subject_atom',
		'target_resname','target_resid','target_atom']

def make_legend(ax,replicate_mapping_this,ncol=1,keys=None):
	"""Make a legend."""
	# redundant with plot-lipid_rdfs.pf
	legendspec = []
	for sn_general,sns in replicate_mapping_this:
		if not keys or sn_general in keys:
			legendspec.append(dict(
				#! preference for one-liners is getting the best of me!
				name=extra_labels.get(sn_general,work.meta.get(sn_general,{}).get('label','')),
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.),ncol=ncol)
	frame = legend.get_frame()
	frame.set_edgecolor('k')
	frame.set_facecolor('w')
	return legend

def plot_bonds(name,kinds,**kwargs):
	"""Plot the distributions of observed bond counts."""
	merged = kwargs.get('merged',False)
	normed = kwargs.get('normed',False)
	symmetrize = kwargs.get('symmetrize',False)
	# fancy imshow method is a great substitute for error bars
	style = kwargs.get('style',['bars','imshow'][-1])
	#! hardcoded settings
	do_fancy_interp = True
	fancy_interp_max = ['absolute','global','local'][1]
	#! only 10 bins look good if you do not do the advanced interpolation method
	nbins_normed = 10.
	norm_global_max = False
	# aesthetics
	wspace = kwargs.get('wspace',0.4)
	hspace = kwargs.get('hspace',0.2)
	legend_ncols = kwargs.get('legend_ncols',1)
	#! ignoring POPC for now
	resnames_exclude = ['POPC']
	def get_norm_factor(sn,combo):
		if normed:
			nmols = [nmol_counts.get(sn,nmol_counts[dict(replicate_mapping).get(sn,[sn])[0]])[r] 
				for r in combo]
			norm_factor = nmols[0]*nmols[1]
		else: norm_factor = 1.0
		return norm_factor
	global post
	post = {}
	global replicate_mapping_this
	if merged: replicate_mapping_this = replicate_mapping 
	else: replicate_mapping_this = [(sn,[sn]) for sn in sns]
	# tabulate nframes by simulation
	nframes_by_sn = dict([(sn,len(data['salt_bridges'][sn]['data']['observations'])) for sn in sns])
	valid_frames_by_kind = dict([(sn,dict([(k,data[k][sn]['data']['valid_frames']) for k in kinds])) for sn in sns])
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
	if symmetrize:
		post_symmetric = {}
		for sn in post:
			post_symmetric[sn] = {}
			combos_sym = list(set([tuple(sorted(i)) for i in post[sn].keys()]))
			for combo in combos_sym:
				post_symmetric[sn][combo] = np.sum([post[sn][c] for c in post[sn] if set(c)==set(combo)],axis=0)
		import ipdb;ipdb.set_trace()
		post = post_symmetric
	# reformulate post into images
	combos_u = list(set([k for j in [i.keys() for i in post.values()] for k in j]))
	sns_by_mapping = list(zip(*replicate_mapping_this))[0]
	try: max_count = int(np.concatenate([np.concatenate(post[sn].values()) for sn in post]).max())
	except: max_count = int(np.concatenate([np.array(post[sn].values()).reshape(-1) for sn in post]).max())
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
				#! alternate count 
				vals_u,idx,counts_u = np.unique(counts_normed,return_counts=True,return_index=True)
				obs_by_snum.append(counts_normed)
				post_norm_interp[combo][sn] = dict(counts_normed=counts_normed,
					vals_u=vals_u,counts_u=counts_u,idx=idx,norm_factor=norm_factor,counts=counts)
			post_norm[combo] = obs_by_snum
			#!np.concatenate(obs_by_snum).max()
		#! simple image method which may be deprecated by the superior interpolation method
		#! THIS SECTION IS DEPRECATED. Use the "fancy" method even when not doing normalization
		if not do_fancy_interp:
			#! changes at the end of the following are incorrect.
			# get the global maximum of the normed counts
			max_norm = max([np.concatenate(i).max() for i in post_norm.values()])
			bins_global = np.linspace(0,max_norm,nbins_normed)
			images = {}
			for combo in combos_u:
				sns_this = [sn for sn in sns_by_mapping if combo in post[sn]]
				images[combo] = []
				norm_max_local = np.concatenate(post_norm[combo]).max()
				bins_local = np.linspace(0,norm_max_local,nbins_normed)
				for item,sn in zip(post_norm[combo],sns_this):
					"""
					in this standard method we use uniform bins to make a uniform image however this fails
					to look reasonable because the distributions with the largest normalization factor
					end up looking like steps with nothing between. to put it another way, when we take a few 
					columns of integer values and divide them by a different float to normalize them, they still 
					take a very small number of values, but placing them all in uniform bins makes the ones that 
					have shrunk the most or have the largest normalization factors to occupy discrete steps that 
					the others do not. they then end up taking discrete values on a continuum, leaving a lot of 
					whitespace in between. this is especially jarring when we use one maximum per combination.
					current solution is to use a small number of bins to handle the problem. proposed best 
					solution is to do an elaborate interpolation which maintains the area density of the counts
					to avoid being misleading by making thinner columns for the simulations with the highest
					amounts of normalization, and then eliminating whitespace by mapping the discrete normalized 
					values onto the full number of bins somehow. the advanced method will really help with
					the clarity of the plot. 20 bins looks better than 10, but starts to have the whitespace 
					problem while 10 looks indistinct and more like a rorschach. this motivated the advanced 
					interpolation method
					"""
					# a second problem is that the global max shows lots of high-likelihood PIP2-PIP2 bonds
					# ... for the 10% system so all the other tiles are mostly white. it would be good to plot 
					# ... both with and without the global max for comparison
					if norm_global_max: bins = bins_global
					else: bins = bins_local
					hist,_ = np.histogram(item,bins=bins,normed=True)
					images[combo].append(hist)
				images[combo] = np.array(images[combo])
				# reformulate by simulation color (sorry this is messy. written quickly but works)
				images[combo] = (np.concatenate((np.tile(np.ones(images[combo].shape),(3,1,1)),
					[images[combo]]))*np.tile(np.array(mpl.colors.to_rgba(color_by_simulation(sn))),
					(images[combo].shape[0],images[combo].shape[1],1)).transpose((2,0,1))).transpose((2,1,0))
		# loop over normed data and interpolate. note that iterpolation is *essential* to rendering
		# ... this data correctly because we have integer counts normalized by floats (see note above)
		else:
			images = {}
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
					xpos += width_bins
				# transpose so the image looks correct on imshow using the 3-dimensional color values
				images[combo] = np.transpose(image,(1,0,2))/image.max()
		max_row_global = max([max(np.where(np.sum(v[:,:,:3].sum(axis=2),axis=1)>0)[0]) for v in images.values()])
	# plot
	axes,fig = square_tiles(len(combos_u),figsize=16,wspace=wspace,hspace=hspace)
	for cnum,combo in enumerate(combos_u):
		ax = axes[cnum]
		if style=='bars':
			counter = 0
			for sn in sns_this:
				if combo in post[sn]:
					norm_factor = get_norm_factor(sn,combo)
					mean = post[sn][combo].mean()/norm_factor
					std = post[sn][combo].std()/norm_factor
					ax.bar([counter],mean,width=1,color=color_by_simulation(sn))
					#! error bar is really huge hence the distribution is not normal. previously:
					#! ... ax.errorbar([counter],mean-std,mean+std,xerr=0,
					#! ... zorder=3,fmt='none',ecolor='k',lw=0.5)
					counter += 1
		elif style=='imshow':
			cmap_name = 'binary'
			if do_fancy_interp and fancy_interp_max=='global': 
				raw = images[combo][:max_row_global]
			elif do_fancy_interp and fancy_interp_max=='local':
				# find the maximum row with a nonzero element (excluding opacity)
				# ... note that this line was written quickly using a sum trick over columns
				max_row_this = max(np.where(np.sum(image[:,:,:3].sum(axis=2),axis=1)>0)[0])
				raw = images[combo][:max_row_this]
			else: raw = images[combo]
			#! interesting note: had to do the transpose when the image is made after switching to 4-d
			#! ... image values and then I had to change axis commands and array slicing and the shape 
			#! ... indexing in the aspect ratio in five places to correct this single upstream transpose\
			kwargs_to_imshow = {}
			if normed: kwargs_to_imshow['extent'] = [0.,1.,0.,1.]
			ax.imshow(raw,origin='lower',interpolation='nearest',
				cmap=mpl.cm.__dict__[cmap_name],**kwargs_to_imshow)
			if style=='imshow' and do_fancy_interp: 
				if normed: ax.set_aspect(1.0)
				else: ax.set_aspect(float(raw.shape[1])/raw.shape[0])
			else: ax.set_aspect(float(len(sns_this))/(max_count if not normed else nbins_normed))
		else: raise Exception
		ax.set_title('%s-%s'%tuple([work.vars['names']['short'].get(p,p) for p in combo]))
		if normed: ax.set_ylabel('score')
		else: ax.set_ylabel('bonds')
		ax.set_xticks([])
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	extras = [make_legend(axes[-1],replicate_mapping_this=replicate_mapping_this,
		keys=[i for i,j in replicate_mapping_this],ncol=legend_ncols)]
	title_names = {'hydrogen_bonding':'hydrogen bonds','salt_bridges':'salt bridges'}
	extras.append(plt.suptitle(' and '.join(title_names[k] for k in kinds)+
		(' (normalized)' if normed else ''),fontsize=20,y=0.92))
	picturesave('fig.lipid_lipid_bonds.%s.%s'%(style,name),work.plotdir,
		backup=False,version=True,meta={},extras=extras)
	post_debugger[tuple(kinds)] = post

def plots_actinlink():
	"""Plot actinlink bonds analysis including merged replicates."""
	figspec = {}
	global post_debugger
	post_debugger = {}
	kind_map = {'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],
		'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	for normed in [True,False]:
		for merged in [True,False]:
			for key,kind in kind_map.items():
				name = key+('.merged' if merged else '')+('.normed' if normed else '')
				figspec[name] = {'merged':merged,'kinds':kind,'normed':normed}
	#! testing figspec = {'salt.normed':figspec['salt.normed']}
	for key,spec in figspec.items(): plot_bonds(name=key,**spec)

def plots_ptdins():
	"""Plot actinlink bonds analysis including merged replicates."""
	figspec = {}
	global post_debugger
	post_debugger = {}
	kind_map = {'hbonds':['hydrogen_bonding'],'salt':['salt_bridges'],'hbonds_salt':['hydrogen_bonding','salt_bridges']}
	# no need to norm because the concentrations are identical
	normed = False
	merged = False
	symmetrize = True
	for key,kind in kind_map.items():
		name = key+('.normed' if normed else '')
		figspec[name] = {'merged':merged,'kinds':kind,'normed':normed,'symmetrize':symmetrize}
	#! testing figspec = {'salt.normed':figspec['salt.normed']}
	for key,spec in figspec.items(): plot_bonds(name=key,**spec)

@autoplot(plotrun)
def plots():
	"""Project-specific plotting sets."""
	project_name = os.path.basename(os.getcwd())
	globals()['plots_%s'%project_name]()	

plotrun.routine = None
if __name__=='__main__': pass

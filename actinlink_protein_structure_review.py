#!/usr/bin/env python

if 'data' not in globals():

	import MDAnalysis
	import brewer2mpl
	import tempfile,subprocess,re,time

	def align(r0,r1):
		r1 -= np.mean(r1,axis=0)
		U,s,Vt = np.linalg.svd(np.dot(r0.T,r1))
		signer = np.identity(3)
		signer[2,2] = np.sign(np.linalg.det(np.dot(Vt.T,U)))
		RM = np.dot(np.dot(U,signer),Vt)
		return RM

	sns,(data,calc) = work.sns(),plotload('protein_rmsd',work)
	post = {}
	for sn in sns:
		if sn in post: continue
		fn = 'rep.%s.gro'%sn
		gro,xtc = [os.path.join(work.postdir,'%s.%s'%(calc['extras'][sn]['slice_path'],s)) 
			for s in ['gro','xtc']]
		uni = MDAnalysis.Universe(gro,xtc)
		sel = uni.select_atoms('name CA')
		nframes = len(uni.trajectory)
		# load coordinates
		pos = np.zeros((nframes,len(sel),3))
		for fr in range(nframes):
			status('reading %s'%sn,looplen=nframes,i=fr,tag='load')
			uni.trajectory[fr]	
			sel = uni.select_atoms('name CA')
			pos[fr] = sel.positions
		post[sn] = dict(coords=pos,uni=uni)

	# deprecated method. it prints the structures but fails to align them properly, so just use VMD
	if False:

		for sn in sns:
			status('aligning within simulation: %s'%sn)
			nframes = len(post[sn]['coords'])
			post[sn]['coords'][0] -= post[sn]['coords'][0].mean(axis=0)
			for fr in range(nframes)[1:]:
				post[sn]['coords'][fr] -= post[sn]['coords'][fr].mean(axis=0)
				rm = align(post[sn]['coords'][0],post[sn]['coords'][fr])
				post[sn]['coords'][fr] = np.dot(rm,post[sn]['coords'][fr].T).T
		ref_key = post.keys()[0]
		for sn in post.keys():
			fn = 'rep.%s.pdb'%sn
			coords = post[sn]['coords']
			fr_near = (np.linalg.norm(coords-coords.mean(axis=0),axis=2)**2).sum(axis=1).argmin()
			if sn==ref_key: 
				fr_near_ref = fr_near
			uni = post[sn]['uni']
			sel = uni.select_atoms('protein')
			uni.trajectory[fr_near]
			sel.positions -= sel.positions.mean(axis=0)
			if sn!=ref_key:
				sel.positions = np.dot(align(
					post[ref_key]['coords'][fr_near_ref],post[sn]['coords'][fr_near]),
					sel.positions.T).T
			sel.write(fn)

def compute_secondary_structure(sn,uni,fn_out):
	"""Call DSSP to get secondary structure."""
	sel = uni.select_atoms('protein')
	nframes = len(uni.trajectory)
	structure = []
	start_time = time.time()
	coder = dict(zip('GHITEBSC ',range(9)))
	coder[''] = 9
	"""
	G = 3-turn helix (310 helix). Min length 3 residues.
	H = 4-turn helix (α helix). Minimum length 4 residues.
	I = 5-turn helix (π helix). Minimum length 5 residues.
	T = hydrogen bonded turn (3, 4 or 5 turn)
	E = extended strand in parallel and/or anti-parallel β-sheet conformation. Min length 2 residues.
	B = residue in isolated β-bridge (single pair β-sheet hydrogen bond formation)
	S = bend (the only non-hydrogen-bond based assignment).
	C = coil (residues which are not in any of the above conformations).
	"""
	for fr in range(nframes):
		uni.trajectory[fr]
		status('reading secondary structure %s'%sn,i=fr,looplen=nframes,start=start_time)
		temp_file = tempfile.NamedTemporaryFile(suffix='.pdb',delete=False)
		pdb = MDAnalysis.Writer(temp_file.name,multiframe=False)
		pdb.write(sel)
		pdb.close()
		output = subprocess.check_output(['dssp',temp_file.name])
		parsed = re.findall('#  RESIDUE(.*?)\n(.+)',output,re.M+re.DOTALL)[0][1]
		ss_frame = [coder[i[16:17]] for i in parsed.split('\n')]
		structure.append(ss_frame)
	return np.array(structure)

def make_legend(**kwargs):
	letters = kwargs.get('letters','GHITEBSC ')
	legendspec = []
	for key in letters:
		num = coder[key]
		legendspec.append(dict(name=coder_to_name[key],
			patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_list[num],lw=1,linestyle='-',ec='k')))
	patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
	legend = ax.legend(patches,labels,loc='lower left',
		bbox_to_anchor=(1.0,0.0,1.,1.),ncol=1,fontsize=14)
	frame = legend.get_frame()
	frame.set_edgecolor('k')
	frame.set_facecolor('white')
	return legend

# compute secondary structure
for sn in post:
	fn = '%s.out'%sn
	if os.path.isfile(os.path.join(work.plotdir,fn)): continue
	structure = compute_secondary_structure(sn,post[sn]['uni'],'%s.out'%sn)
	np.savetxt(os.path.join(work.plotdir,fn),np.array(structure).astype(int))

sns_groups = [
	['mdia2bilayer_nochl2','mdia2bilayer_nochl3'],
	['mdia2bilayer10','mdia2bilayer10_2',],
	['mdia2bilayerphys','mdia2bilayerphys2'],
	['mdia2bilayer30','mdia2bilayer30_2']]

routine = ['structure_maps','bars','bars_breakdown'][:]
panel_kwargs = dict(figsize=(16,8),layout={'out':{'grid':[1,4]},'ins':[{'grid':[2,1]} 
	for i in range(4)]})
color_list = brewer2mpl.get_map('Set1','qualitative',9).mpl_colors
color_list[-1] = (1.,1.,1.)
cmap_custom = mpl.colors.ListedColormap(color_list)

if 'structure_maps' in routine:

	extras = []
	axes,fig = panelplot(**panel_kwargs)
	for gnum,grp in enumerate(sns_groups):
		for cnum,sn in enumerate(grp):
			ax = axes[gnum][cnum]
			dat = np.loadtxt('%s.out'%sn)
			ax.set_yticks(range(dat.shape[1]))
			ax.imshow(dat.T,
				interpolation='nearest',aspect=dat.shape[0]/dat.shape[1],cmap=cmap_custom)
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_title(work.meta[sn].get('label',sn))
			resnames = post[sn]['uni'].select_atoms('protein').residues.resnames
			ax.set_yticklabels([residue_codes[j] for j in resnames])
	coder_to_name = {'G':'helix (3-10)','H':'helix','I':'helix (pi)','T':'turn','E':'Beta sheet (extended)',
		'B':'Beta sheet (bridge)','S':'bend','C':'coil',' ':'N/A','':'N/A'}
	coder = dict(zip('GHITEBSC',range(8)),**{' ':8,'':8})
	coder[''] = 8
	extras.append(make_legend())
	picturesave('fig.protein_structure.timeseries',directory=work.plotdir,extras=extras)

if 'bars' in routine:

	extras,kinds_all,y_max = [],[],0.0
	axes,fig = panelplot(**panel_kwargs)
	for gnum,grp in enumerate(sns_groups):
		for cnum,sn in enumerate(grp):
			ax = axes[gnum][cnum]
			# load and reduce
			dat = np.loadtxt('%s.out'%sn)
			dat[np.where(dat==9)] = 8
			kinds = np.unique(dat).astype(int)
			kinds_all = list(set(list(kinds)+kinds_all))
			kinds_indexer = dict(zip(kinds,range(len(kinds))))
			dat_reduce = np.zeros((len(dat),len(kinds)))
			for rnum,row in enumerate(dat):
				vals,counts = np.unique(row,return_counts=True)
				inds,counts = tuple(np.array([(kinds_indexer[i],j) for i,j in zip(vals,counts)]).T)
				dat_reduce[rnum,inds] = counts
			ax.bar(np.arange(len(kinds)),dat_reduce.mean(axis=0),width=1.0,
				color=[color_list[i] for i in kinds])
			y_max = max([y_max,dat_reduce.mean(axis=0).max()])
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_title(work.meta[sn].get('label',sn))
			ax.set_ylabel('residues per frame')
			letters = dict([i[::-1] for i in coder.items()])
			letters[8] = 'N/A'
			ax.set_xticks(np.arange(len(kinds)))
			ax.set_xticklabels([letters[i] for i in kinds])
			ax.set_xlim((-0.5,len(kinds)-0.5))
	for ax in [i for j in axes for i in j]: ax.set_ylim(0,y_max*1.05)
	extras.append(make_legend(
		letters=''.join([dict([(j,i) for i,j in coder.items()])[k] for k in kinds_all])))
	picturesave('fig.protein_structure.breakdown_types',directory=work.plotdir,extras=extras)

if 'bars_breakdown' in routine:

	extras,kinds_all = [],[]
	axes,fig = panelplot(**panel_kwargs)
	for gnum,grp in enumerate(sns_groups):
		for cnum,sn in enumerate(grp):
			status('reducing for %s'%sn)
			ax = axes[gnum][cnum]
			# load and reduce
			dat = np.loadtxt('%s.out'%sn)
			dat[np.where(dat==9)] = 8
			kinds = np.unique(dat).astype(int)
			kinds_all = list(set(list(kinds)+kinds_all))
			kinds_indexer = dict(zip(kinds,range(len(kinds))))
			nres = dat.shape[1]
			dat_reduce = np.zeros((len(kinds),dat.shape[0],nres))
			for rnum,row in enumerate(dat):
				resid_to_kind = [kinds_indexer[i] for i in row.astype(int)]
				dat_reduce[tuple(
					np.array([np.array(resid_to_kind),
					np.ones(nres)*rnum,
					np.arange(nres)]).astype(int))] += 1
			img_dat = dat_reduce.mean(axis=1)
			#ax.imshow(img_dat,aspect=img_dat.shape[1]/img_dat.shape[0])
			vals_prev = np.zeros(nres)
			for knum,vals in enumerate(img_dat.cumsum(axis=0)):
				ax.bar(np.arange(nres),vals-vals_prev,
					bottom=vals_prev,color=color_list[kinds[knum]],width=1.0)
				vals_prev = vals
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_title(work.meta[sn].get('label',sn))
			ax.set_ylabel('frequency')
			ax.set_xlim((-0.5,nres-0.5))
			ax.set_xticks(range(nres))
			resnames = post[sn]['uni'].select_atoms('protein').residues.resnames
			ax.set_xticklabels([residue_codes[j] for j in resnames])
	extras.append(make_legend(
		letters=''.join([dict([(j,i) for i,j in coder.items()])[k] for k in kinds_all])))
	picturesave('fig.protein_structure.breakdown_residues',directory=work.plotdir,extras=extras)

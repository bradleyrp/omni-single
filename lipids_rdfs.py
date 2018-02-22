#!/usr/bin/env python

"""
Analyze the lipid-lipid RDFs.
"""

import scipy
import scipy.spatial
from base.compute_loop import basic_compute_loop
import MDAnalysis
import itertools

def compute_rdf(sn,fr,out=[0,0]):
	"""Parallel compute function."""
	global dat,pts2d,dists,vecs_lib
	#vec = dat[sn]['data']['vecs'][fr]
	vec = vecs_lib[sn][fr]
	boxstuff = lambda pts,vec,d=3 : pts[...,:d]-(pts[...,:d]>vec[:d])*vec[:d]+\
		(pts[...,:2]<np.array([0.,0.,0.][:d]))*vec[:d]
	pts_stuff = boxstuff(pts2d[fr],vec,d=2)
	pd = np.array([scipy.spatial.distance.cdist(np.array([pts_stuff[:,d]]).T,
		np.array([pts_stuff[:,d]]).T) for d in range(2)])
	for d in range(2):
		pd[d] -= (pd[d]>vec[d]/2.)*vec[d]
		pd[d] += (pd[d]<-1*vec[d]/2.)*vec[d]
	if out=='xy':
		for d in [0,1]:
			pd[d] += vec[d]
	elif out=='x':
		for d in [0]:
			pd[d] += vec[d]
	elif out=='y':
		for d in [1]:
			pd[d] += vec[d]
	return np.linalg.norm(pd,axis=0)

@autoload(plotrun)
def load():
	data,calc = plotload('lipid_abstractor')
	def mesh_style(calc,abstractor): 
		"""Check if the incoming lipid abstractor data matches the abstractor key."""
		if abstractor=='lipid_chol_com':
			return 'CHL1' in calc['calcs']['specs']['selector']['resnames']
		elif abstractor=='lipid_com':
			return 'CHL1' not in calc['calcs']['specs']['selector']['resnames']
		else: raise Exception
	sns = work.sns()

#! debugging the long-range part
if __name__=='__main__' and False:

	cutoff = 20.0
	binsize = 0.05
	scanrange = np.arange(0,cutoff,binsize)

	sns_this = work.sns()[-1:]

	#! vector saving error!
	if 'vecs_lib' not in globals():
		vecs_lib = {}
		for sn in sns_this:
			uni = MDAnalysis.Universe(*[os.path.join(work.postdir,
				'%s.%s'%(calc['extras'][sn]['slice_path'],i)) for i in ['gro','xtc']])
			vecs_alt = []
			nframes = len(uni.trajectory)
			for fr in range(nframes):
				status('reading',i=fr,looplen=nframes)
				#! you have to recast this otherwise you get repeats! why?!
				vecs_alt.append(np.array(uni.trajectory[fr].dimensions[:3]))
			vecs_alt = np.array(vecs_alt)/10.
			vecs_lib[sn] = vecs_alt

	#! cannot get this into load for now
	if 'dists_all' not in globals():
		dists_all,counted_all = {},{}
		pairspec = {}
		for abstractor in ['lipid_com','lipid_chol_com'][:1]:
			if abstractor=='lipid_com': sns = sns_this
			else: sns = [sn for sn in sns_this if sn not in ['mdia2bilayer_nochl3','mdia2bilayer_nochl2']]
			# prepare a place to hold distances
			dists,resnames,imono = {},{},{}
			shifts = [[0,0],[0,1],[1,0],[1,1]]#! [:1] #!?
			shifts = ['','x','y','xy']
			#! note that we have to look things up by resname if we are not pulling an upstream
			#! ... calculation since the upstream is what allows us to look things up 
			#! ... by the internal loop name
			key_select_calc, = [key for key,val in calc.items() if 'calcs' in val 
				and mesh_style(calc=val,abstractor=abstractor)]
			dat = data[key_select_calc]
			for sn in sns_this:
				status('computing RDF for %s'%sn)
				imono = dat[sn]['data']['monolayer_indices']
				points = dat[sn]['data']['points']
				top_mn = actinlink_monolayer_indexer(sn=sn,abstractor=abstractor)
				resnames[sn] = dat[sn]['data']['resnames'][np.where(imono==top_mn)]
				pts2d = points[:,np.where(imono==top_mn)[0],:2]
				nframes,nmol = pts2d.shape[:2]
				looper = [dict(fr=fr,sn=sn,out=out) for fr in range(nframes) 
					for out in shifts]
				dists[sn] = np.zeros((nframes,nmol,nmol))
				incoming = basic_compute_loop(compute_rdf,looper)
				dists[sn] = np.array(incoming)
			# counting sequence starts here
			resnames_all = np.unique(np.concatenate([val for key,val in resnames.items()]))
			pairs = ([([r],[r]) for r in resnames_all]+
				[([i],[j]) for i,j in itertools.combinations(resnames_all,2)]+
				[(resnames_all,resnames_all)])
			pairnames = [(i[0],j[0]) for i,j in pairs[:-1]]+[('all lipids','all lipids')]
			pairspec[abstractor] = dict(pairs=pairs,pairnames=pairnames)

			#! override for this debug seq
			pairnames = pairnames[-1:]
			pairs = pairs[-1:]

			def count_rdf_reduce(vals):
				global scanrange
				return np.histogram(vals,bins=scanrange)[0] 
			counted = dict([(sn,{}) for sn in sns_this])
			for sn in sns_this:
				status('reducing %s'%sn)
				for pairname,pair in zip(pairnames,pairs):
					group_1 = np.where(np.in1d(resnames[sn],pair[0]))[0]
					group_2 = np.where(np.in1d(resnames[sn],pair[1]))[0]
					incoming = basic_compute_loop(count_rdf_reduce,[
						dict(vals=i) for i in np.concatenate(dists[sn][:,group_1][...,group_2])])
					counted[sn][pairname] = np.array(incoming)
			dists_all[abstractor] = dists
			counted_all[abstractor] = counted

	axes,fig = square_tiles(1,figsize=(8,8),hspace=0.4,wspace=0.4)
	sn = sns_this[0]

	middles = (scanrange[1:]+scanrange[:-1])/2
	areas = np.array([np.pi*(binsize)*middles[i]*2 for i in range(len(middles))])

	raw = np.array(counted_all[abstractor][sn].values()[0])
	raw[np.where(raw[:,0]==1),0] -= 1
	counts = raw.mean(axis=0)
	ax = axes[0]
	density = 1.0
	#density = nmol/np.mean(total_area)
	factor_nonlike = 1.0
	ax.plot(middles,counts/areas/density*factor_nonlike,color='k',lw=2)

	picturesave('fig.DEBUG',directory=work.plotdir,meta={},extras=[])

	pass

if __name__=='__main__': 

	cutoff = 20.0
	binsize = 0.05
	scanrange = np.arange(0,cutoff,binsize)

	#! vector saving error!
	if 'vecs_lib' not in globals():
		vecs_lib = {}
		for sn in sns:
			uni = MDAnalysis.Universe(*[os.path.join(work.postdir,
				'%s.%s'%(calc['extras'][sn]['slice_path'],i)) for i in ['gro','xtc']])
			vecs_alt = []
			nframes = len(uni.trajectory)
			for fr in range(nframes):
				status('reading',i=fr,looplen=nframes)
				#! you have to recast this otherwise you get repeats! why?!
				vecs_alt.append(np.array(uni.trajectory[fr].dimensions[:3]))
			vecs_alt = np.array(vecs_alt)/10.
			vecs_lib[sn] = vecs_alt

	#! cannot get this into load for now
	if 'dists_all' not in globals():
		dists_all,counted_all = {},{}
		pairspec = {}
		for abstractor in ['lipid_com','lipid_chol_com']:
			if abstractor=='lipid_com': sns = work.sns()
			else: sns = [sn for sn in work.sns() if sn not in ['mdia2bilayer_nochl3','mdia2bilayer_nochl2']]
			# prepare a place to hold distances
			dists,resnames,imono = {},{},{}
			shifts = [[0,0],[0,1],[1,0],[1,1]][:1] #!?
			#! note that we have to look things up by resname if we are not pulling an upstream
			#! ... calculation since the upstream is what allows us to look things up 
			#! ... by the internal loop name
			key_select_calc, = [key for key,val in calc.items() if 'calcs' in val 
				and mesh_style(calc=val,abstractor=abstractor)]
			dat = data[key_select_calc]
			for sn in sns:
				status('computing RDF for %s'%sn)
				imono = dat[sn]['data']['monolayer_indices']
				points = dat[sn]['data']['points']
				top_mn = actinlink_monolayer_indexer(sn=sn,abstractor=abstractor)
				resnames[sn] = dat[sn]['data']['resnames'][np.where(imono==top_mn)]
				pts2d = points[:,np.where(imono==top_mn)[0],:2]
				nframes,nmol = pts2d.shape[:2]
				looper = [dict(fr=fr,sn=sn,out=out) for fr in range(nframes) 
					for out in shifts]
				dists[sn] = np.zeros((nframes,nmol,nmol))
				incoming = basic_compute_loop(compute_rdf,looper)
				dists[sn] = np.array(incoming)
			# counting sequence starts here
			resnames_all = np.unique(np.concatenate([val for key,val in resnames.items()]))
			pairs = ([([r],[r]) for r in resnames_all]+
				[([i],[j]) for i,j in itertools.combinations(resnames_all,2)]+
				[(resnames_all,resnames_all)])
			pairnames = [(i[0],j[0]) for i,j in pairs[:-1]]+[('all lipids','all lipids')]
			pairspec[abstractor] = dict(pairs=pairs,pairnames=pairnames)
			def count_rdf_reduce(vals):
				global scanrange
				return np.histogram(vals,bins=scanrange)[0] 
			counted = dict([(sn,{}) for sn in sns])
			for sn in sns:
				status('reducing %s'%sn)
				for pairname,pair in zip(pairnames,pairs):
					group_1 = np.where(np.in1d(resnames[sn],pair[0]))[0]
					group_2 = np.where(np.in1d(resnames[sn],pair[1]))[0]
					incoming = basic_compute_loop(count_rdf_reduce,[
						dict(vals=i) for i in np.concatenate(dists[sn][:,group_1][...,group_2])])
					counted[sn][pairname] = np.array(incoming)
			dists_all[abstractor] = dists
			counted_all[abstractor] = counted

	#! needed to be changed
	def actinlink_monolayer_indexer(sn,abstractor):
		if abstractor=='lipid_chol_com':
			if sn in ['mdia2bilayer10_2','mdia2bilayer30_2']: return 1
			else: return 0
		elif abstractor=='lipid_com':
			if sn in ['mdia2bilayer10_2','mdia2bilayer30_2']: return 1
			else: return 0
		else: raise Exception

	# standard plots
	master_plots = ['lipid_com','lipid_chol_com']
	master_plots = []
	for abstractor in master_plots:
		pairs,pairnames = [pairspec[abstractor][k] for k in ['pairs','pairnames']]
		middles = (scanrange[1:]+scanrange[:-1])/2
		areas = np.array([np.pi*(binsize)*middles[i]*2 for i in range(len(middles))])
		#! move this stuff to the art file and collect from elsewhere!
		replicate_mapping = [('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
			('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
			('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
			('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])][1 if abstractor=='lipid_chol_com' else 0:]
		extra_labels = {
			'pip2_20_no_chol':r'mDia2, 20% $PIP_2$, no CHOL ($\times2$)',
			'pip2_20':r'mDia2, 20% $PIP_2$ ($\times2$)',
			'pip2_30':r'mDia2, 30% $PIP_2$ ($\times2$)',
			'pip2_10':r'mDia2, 10% $PIP_2$ ($\times2$)',}
		#! actinlink settings go in the art file please !!! this is verbatim from actinlink_bonds_analysis
		def color_by_simulation(sn):
			# tetradic colors via https://www.sessions.edu/color-calculator/
			colors = ['#ff6c28','#28c6ff','#a928ff','#ffd828']
			# tetradic colors via http://paletton.com/#uid=7030Z0kqbujggFvlnx6qcY+wDkl
			colors = ['#f2b52c','#22b93c','#2e47a4','#f23a2c']
			refs = ['^mdia2bilayer_nochl','^mdia2bilayer(?!(_nochl|[0-9]))',
				'^mdia2bilayer10','^mdia2bilayer30']
			ref_to_color = dict(zip(refs,colors))
			matches = [color for ref,color in ref_to_color.items() if re.match(ref,sn)]
			if len(matches)>1: raise Exception
			elif len(matches)==0:
				#! handling combos here with a minor hack see global replicate_mapping
				return dict(zip(zip(*replicate_mapping)[0],colors))[sn]
			else: return matches[0]
		#! from 3D: areas = np.array([4*np.pi*binsize*middles[i]**2 for i in range(len(middles))])
		axes,fig = square_tiles(len(pairs),figsize=(12,12),hspace=0.4,wspace=0.4)
		for pnum,(pairname,pair) in enumerate(zip(pairnames,pairs)):
			ax = axes[pnum]
			for sn_general in dict(replicate_mapping):
				status('rendering %s'%sn_general)
				nmols,total_areas = [],[]
				for sn in dict(replicate_mapping)[sn_general]:
					top_mn = actinlink_monolayer_indexer(sn=sn,abstractor=abstractor)
					imono = dat[sn]['data']['monolayer_indices']
					resnames = dat[sn]['data']['resnames'][np.where(imono==top_mn)]
					pair = pairs[pairnames.index(pairname)]
					group_1 = np.where(np.in1d(resnames,pair[0]))[0]
					group_2 = np.where(np.in1d(resnames,pair[1]))[0]
					nmol = float(len(group_1))
					nmols.append(nmol)
					total_area = np.product(vecs_lib[sn].mean(axis=0)[:2])
					total_areas.append(total_area)
				# ensure replicates are like
				if len(set(nmols))!=1: raise Exception
				else: 
					nmol = list(set(nmols))[0]
					sn = dict(replicate_mapping)[sn_general][0]
					factor_nonlike = float(np.in1d(resnames,pair[0]).sum())/\
						np.in1d(resnames,pair[1]).sum()
				density = nmol/np.mean(total_area)
				raw = np.concatenate([counted_all[abstractor][sn][pairname] 
					for sn in dict(replicate_mapping)[sn_general]])
				if pairname[0]==pairname[1]: raw[:,0]-=1
				counts = raw.mean(axis=0)
				ax.plot(middles,counts/areas/density*factor_nonlike,color=color_by_simulation(sn),lw=2)
			ax.axhline(1.0,lw=1,c='k')
			ax.set_xlim((0,3.0))
			ax.set_title('%s-%s'%tuple([work.vars['names']['short'].get(p,p) for p in pairname]))
			ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
			ax.set_ylabel('$g(r)$')
			ax.set_xlabel('$r\,(nm)$')
		# legend at last tile
		legendspec = []
		for sn_general,sns in replicate_mapping:
			legendspec.append(dict(name=extra_labels[sn_general],
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.))
		frame = legend.get_frame()
		frame.set_edgecolor('k')
		frame.set_facecolor('w')
		picturesave('fig.lipids_rdfs.%s'%abstractor,directory=work.plotdir,meta={},extras=[legend])

	specials = {
		'fig.lipids_rdfs.compare_cholesterol':
			dict(abstractor = 'lipid_chol_com',
			sns_special = ['pip2_10','pip2_20','pip2_30'],
			special_pairs = [
				(('DOPS','DOPS'),(['DOPS'],['DOPS'])),
				(('PI2P','PI2P'),(['PI2P'],['PI2P'])),
				(('DOPS','PI2P'),(['DOPS'],['PI2P'])),],),
		'fig.lipids_rdfs.compare_20':
			dict(abstractor = 'lipid_com',
			sns_special = ['pip2_20_no_chol','pip2_20'],
			special_pairs = [
				(('DOPS','DOPS'),(['DOPS'],['DOPS'])),
				(('PI2P','PI2P'),(['PI2P'],['PI2P'])),
				(('DOPS','PI2P'),(['DOPS'],['PI2P'])),],),}
	specials = {'fig.lipids_rdfs.compare_20':specials['fig.lipids_rdfs.compare_20']}

	# special plots show some comparisons we are interested
	for special_name,special_spec in specials.items():
		abstractor = special_spec['abstractor']
		replicate_mapping = [('pip2_20_no_chol',['mdia2bilayer_nochl2','mdia2bilayer_nochl3']),
			('pip2_10',['mdia2bilayer10','mdia2bilayer10_2']),
			('pip2_20',['mdia2bilayerphys','mdia2bilayerphys2']),
			('pip2_30',['mdia2bilayer30','mdia2bilayer30_2'])][1 if abstractor=='lipid_chol_com' else 0:]
		special_pairs = special_spec['special_pairs']
		abstractor = special_spec['abstractor']
		sns_special = special_spec['sns_special']
		axes,fig = square_tiles(len(special_pairs),figsize=(8,8),hspace=0.4,wspace=0.4)
		for pnum,(pairname,pair) in enumerate(special_pairs):
			ax = axes[pnum]
			for sn_general in sns_special:
				status('rendering %s'%sn_general)
				nmols,total_areas = [],[]
				for sn in dict(replicate_mapping)[sn_general]:
					top_mn = actinlink_monolayer_indexer(sn=sn,abstractor=abstractor)
					imono = dat[sn]['data']['monolayer_indices']
					resnames = dat[sn]['data']['resnames'][np.where(imono==top_mn)]
					pair = pairs[pairnames.index(pairname)]
					group_1 = np.where(np.in1d(resnames,pair[0]))[0]
					group_2 = np.where(np.in1d(resnames,pair[1]))[0]
					nmol = float(len(group_1))
					nmols.append(nmol)
					total_area = np.product(vecs_lib[sn].mean(axis=0)[:2])
					total_areas.append(total_area)
				# ensure replicates are like
				if len(set(nmols))!=1: raise Exception
				else: 
					nmol = list(set(nmols))[0]
					sn = dict(replicate_mapping)[sn_general][0]
					factor_nonlike = float(np.in1d(resnames,pair[0]).sum())/\
						np.in1d(resnames,pair[1]).sum()
				density = nmol/np.mean(total_area)
				raw = np.concatenate([counted_all[abstractor][sn][pairname] 
					for sn in dict(replicate_mapping)[sn_general]])
				if pairname[0]==pairname[1]: raw[:,0]-=1
				counts = raw.mean(axis=0)
				ax.plot(middles,counts/areas/density*factor_nonlike,color=color_by_simulation(sn),lw=2)
				ax.axhline(1.0,lw=1,c='k')
				ax.set_xlim((0,3.0))
				ax.set_title('%s-%s'%tuple([work.vars['names']['short'].get(p,p) for p in pairname]))
				ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
				ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
				ax.set_ylabel('$g(r)$')
				ax.set_xlabel('$r\,(nm)$')
		# legend at last tile
		legendspec = []
		for sn_general in sns_special:
			sns = dict(replicate_mapping)[sn_general]
			legendspec.append(dict(name=extra_labels[sn_general],
				patch=mpl.patches.Rectangle((0,0),1.0,1.0,fc=color_by_simulation(sns[0]))))
		patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
		legend = ax.legend(patches,labels,loc='upper left',bbox_to_anchor=(1.0,0.0,1.,1.))
		frame = legend.get_frame()
		frame.set_edgecolor('k')
		frame.set_facecolor('w')
		picturesave(special_name,directory=work.plotdir,meta={},extras=[legend])

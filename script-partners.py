#!/usr/bin/python -i

compsign,project_name = 'mesh','project-ptdins'
execfile('../../header.py')

import itertools

debug = False
get_datlist = True
routine = [
	'partner_persistence',
	][:]

if 'datlist' not in globals() and not debug:

	computer(focus,compute_mesh,headerdat,simdict)
	if get_datlist: datlist = loader(focus,headerdat)

elif debug:

	panel = focus.keys()[0]
	simname = ts = focus[panel].keys()[0]
	timestamp = focus[panel][ts]['time']
	grofile,trajfile = smx.get_slices(ts,simdict,
		timestamp=timestamp,
		wrap=calculations[compsign]['wrap'],
		groupname=calculations[compsign]['groupname'])
	self = mm = MonoMesh(simname,grofile,trajfile,skip_mesher=True,
		metadat=metadat,**calculations['mesh'])
	mm.nframes = 10
	mm.mesher()
		
#---POST
#-------------------------------------------------------------------------------------------------------------

def partner_watcher(sn,*dats,**kwargs):

	dat = dats[0]

	#---from binding_combinator
	resnames = array(dat['resnames'])
	nframes = dat['nframes']
	lipids = array(list(resnames[sort(unique(resnames,return_index=True)[1])]))

	#---replicate some of the binding_combinator and compute_bridging_norms codes for pairs
	nn = 2-1
	combos = array([''.join(j) for j in 
		itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) 
		if sum([int(k) for k in j])==nn+1])
	combonames = [tuple(v) for v in [concatenate([[lipids[ww]]*int(w) 
		for ww,w in enumerate(l)]) for l in  combos]]
	nn = 3-1
	combos3 = array([''.join(j) for j in 
		itertools.product(''.join([str(i) for i in range(nn+2)]),repeat=len(lipids)) 
		if sum([int(k) for k in j])==nn+1])
	combonames3 = [tuple(v) for v in [concatenate([[lipids[ww]]*int(w) 
		for ww,w in enumerate(l)]) for l in  combos3]]

	
	#---determine monolayer-specific residue indices
	imono = dat['imono']
	nmols = [shape(dat['mono'+str(mn)])[1] for mn in range(2)]
	resnames = array(dat['resid_names'])
	reslist = list(array(resnames)[sort(unique(resnames,return_index=True)[1])])
	rxm = [[
		array([where(where(imono==mn)[0]==i)[0][0] 
		for i in where(all((imono==mn,resnames==rn),axis=0))[0]])
		for rn in reslist] for mn in range(2)]
	rxmlook = [zeros(n) for n in nmols]
	for mn in range(2):
		for ri,r in enumerate(rxm[mn]):
			if r != []: rxmlook[mn][r] = ri
		
	#---record the link count trajectories
	combolookup = sum([array(combonames)==r for r in reslist],axis=2).T
	combolookup_str = [''.join(['%s'%s for s in i]) for i in combolookup]
	combolookup3 = sum([array(combonames3)==r for r in reslist],axis=2).T
	combolookup3_str = [''.join(['%s'%s for s in i]) for i in combolookup3]
	tally = [zeros((nframes,len(combolookup))) for mn in range(2)]
	tally_random = [zeros((nframes,len(combolookup))) for mn in range(2)]
	tally3 = [zeros((nframes,len(combolookup3))) for mn in range(2)]
	tally3_random = [zeros((nframes,len(combolookup3))) for mn in range(2)]
	for mn in range(2):
		st = time.time()
		for fr in range(nframes):
			status('[COMPUTE] counting partners mn=%d'%mn,i=fr,looplen=nframes,start=st)
			simplices = dat[i2s(fr,mn,'simplices')]
			gids = dat[i2s(fr,mn,'ghost_ids')]
			links = array(list(set([tuple(sort(j)) for j in 
				concatenate([transpose((gids[simplices][:,i],gids[simplices][:,(i+1)%3])) 
				for i in range(3)])])))
			#---uniquely identify each link type
			links_type = transpose([sum(rxmlook[mn][links]==i,axis=1) for i in range(len(reslist))])
			links_type_str = [''.join(['%s'%s for s in i]) for i in links_type]
			counts = dict(scipy.stats.itemfreq(links_type_str))
			tally[mn][fr] = array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

			#---sidestepping the probability questions and randomizing the identifier
			#---...which will effectively randomize the identities of the vertices
			#---...and then later we will take the average of these and then use it to see if observed
			#---...links are more or less likely to appear in our graph than in a random one
			rxmlook_rand = [random.permutation(r) for r in rxmlook]
			links_type_wrong = transpose([sum(rxmlook_rand[mn][links]==i,axis=1) 
				for i in range(len(reslist))])
			links_type_str = [''.join(['%s'%s for s in i]) for i in links_type_wrong]
			counts = dict(scipy.stats.itemfreq(links_type_str))
			tally_random[mn][fr] = array([int(counts[i]) if i in counts else 0 for i in combolookup_str])

			#---identify each triple type
			triples = gids[simplices]
			triple_type = transpose([sum(rxmlook[mn][triples]==i,axis=1) for i in range(len(reslist))])
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type]
			counts = dict(scipy.stats.itemfreq(triple_type_str))
			tally3[mn][fr] = array([int(counts[i]) if i in counts else 0 for i in combolookup3_str])

			#---randomize the triples
			triple_type_wrong = transpose([sum(rxmlook_rand[mn][triples]==i,axis=1) 
				for i in range(len(reslist))])
			triple_type_str = [''.join(['%s'%s for s in i]) for i in triple_type_wrong]
			counts = dict(scipy.stats.itemfreq(triple_type_str))
			tally3_random[mn][fr] = array([int(counts[i]) if i in counts else 0 for i in combolookup3_str])

	results = {}
	results ['nframes'] = array(nframes)
	results ['resnames'] = array(resnames).astype(str)
	for mn in range(2): results ['pairs%d'%mn] = array(tally[mn])
	for mn in range(2): results ['pairs_random%d'%mn] = array(tally_random[mn])
	for mn in range(2): results ['triples%d'%mn] = array(tally3[mn])
	for mn in range(2): results ['triples_random%d'%mn] = array(tally3_random[mn])
	results ['combonames'] = array(combonames).astype(str)
	results ['combonames3'] = array(combonames3).astype(str)
	return results
				
compute_post_post_basic(
	'partners',partner_watcher,
	headerdat=headerdat,focus=focus,
	compsign_original='mesh',
	upstream=['mesh'])

#---MAIN
#-------------------------------------------------------------------------------------------------------------

if 'partner_persistence' in routine:

	import matplotlib.patheffects as path_effects

	mn = 0

	for bridge_type in ['pairs','triples']:


		sns = list(unique([i for j in [comparisons[k] for k in focus] for i in j]))
		if 'postdat' not in globals():
			postdat = loader(focus,headerdat,compsign='partners',compsign_original='bridging')
		#---for each simulation store the rxmlook for later
		extradat = dict([(sn,{}) for sn in sns])
		for si,sn in enumerate(sns):
			combos = postdat[sn][{'pairs':'combonames','triples':'combonames3'}[bridge_type]]
			dat = datlist[sn]
			composition = dict([(i[0],int(i[1])) for i in scipy.stats.itemfreq(datlist[sn]['resid_names'])])
			nmols = [shape(dat['mono'+str(mn)])[1] for mn in range(2)]
			#---determine monolayer-specific residue indices
			imono = dat['imono']
			nmols = [shape(dat['mono'+str(mn)])[1] for mn in range(2)]
			resnames = array(dat['resid_names'])
			reslist = list(array(resnames)[sort(unique(resnames,return_index=True)[1])])
			rxm = [[
				array([where(where(imono==mn)[0]==i)[0][0] 
				for i in where(all((imono==mn,resnames==rn),axis=0))[0]])
				for rn in reslist] for mn in range(2)]
			rxmlook = [zeros(n) for n in nmols]
			for mn in range(2):
				for ri,r in enumerate(rxm[mn]):
					if r != []: rxmlook[mn][r] = ri
			extradat[sn]['rxmlook'] = rxmlook
			extradat[sn]['reslist'] = reslist	
		
		mn = 0

		#---compute all of the mean numbers of links
		counter = dict([(sn,{}) for sn in sns])
		for sn in sns:
			rxmlook = extradat[sn]['rxmlook']
			reslist = extradat[sn]['reslist']		
			composition = dict([(i[0],int(i[1])) for i in scipy.stats.itemfreq(datlist[sn]['resid_names'])])
			combonames = postdat[sn][{'pairs':'combonames','triples':'combonames3'}[bridge_type]]
			for code in ['','_random']:
				counts = [mean(postdat[sn]['%s%s%d'%(bridge_type,code,mn)][:,cc]) 
					for cc,c in enumerate(combonames)]
				counter[sn]['count%s'%code] = array(counts)
			#---running window average for compute the starting and ending values
			winsize = 10
			first_last = array([mean(postdat[sn]['%s%s%d'%(bridge_type,'',mn)][win,:],axis=0) 
				for win in [arange(winsize),-1*arange(1,winsize+1)]])
			first =  postdat[sn]['%s%s%d'%(bridge_type,'',mn)][0,:]
			counter[sn]['change'] = (first_last[1]-first_last[0])/counter[sn]['count_random']
			counter[sn]['first'] = first/counter[sn]['count_random']
			counter[sn]['normcount'] = nan_to_num(counter[sn]['count']/counter[sn]['count_random'])

		#---organize all combinations so that we group by link type across simulations
		allcombos = list(set([tuple(j) for j in concatenate([postdat[sn][
			{'pairs':'combonames','triples':'combonames3'}[bridge_type]] for sn in sns])]))
		all_counts,all_counts_change,all_counts_start = {},{},{}
		for cc,combo in enumerate(allcombos):
			pops,pops_change,pops_start = [],[],[]
			for sn in sns:
				combonames = [list(i) for i in postdat[sn][
					{'pairs':'combonames','triples':'combonames3'}[bridge_type]]]
				pairsums = [sum(extradat[sn]['rxmlook'][mn]==extradat[sn]['reslist'].index(c)) 
					if c in extradat[sn]['reslist'] else 0 for c in combo]
				if list(combo) in combonames and not any([i==0 for i in pairsums]):
					pops.append(counter[sn]['normcount'][combonames.index(list(combo))])
					pops_change.append(counter[sn]['change'][combonames.index(list(combo))])
					pops_start.append(counter[sn]['first'][combonames.index(list(combo))])
				else: 
					pops.append(0)
					pops_change.append(0)
					pops_start.append(0)
			if any([i>0 for i in pops]): 
				all_counts[combo] = pops
				all_counts_change[combo] = pops_change
				all_counts_start[combo] = pops_start

		relabel = dict(list(set([(metadat[sn]['ptdins_resname'],metadat[sn]['ptdins_label']) 
			for sn in sns])))

		for plot_type in ['relative','relative_arrow','change'][1:2]:

			ylims = {
				'pairs':{
					'relative':(-0.6,0.5),
					'relative_arrow':(-0.6,0.5),
					'change':(-1.,0.6)}[plot_type],
				'triples':{
					'relative':(-1.75,3),
					'relative_arrow':(-1.75,3),
					'change':(-1.,0.6)}[plot_type],
					}[bridge_type]	
			relabel = dict(list(set([(metadat[sn]['ptdins_resname'],metadat[sn]['ptdins_label']) 
				for sn in sns])))
			fig,gs,axes = generate_axes(nrows=1,ncols=1,figsize=({'pairs':16,'triples':20}[bridge_type],6))
			xlabels,xticks = [],[]
			withchar = ' and '
			barparse = 1
			ax = axes[0][0]
			for cc,combo in enumerate(all_counts):

				if plot_type in ['relative','relative_arrow']: dat = all_counts[combo]
				elif plot_type == 'change': dat = all_counts_change[combo]
				nzs = [dd for dd,d in enumerate(dat) if d!=0.]
				xvals = arange(barparse,barparse+len(nzs))
				yvals = [dat[dd]+{'relative':-1,'relative_arrow':-1,'change':0}[plot_type] for dd in nzs]
				alphas = [0.65 if metadat[sns[n]]['composition_name']=='symmetric' else 1. for n in nzs]
				hatching = ['x' if metadat[sns[n]]['composition_name']=='symmetric' else None for n in nzs]
				colors= [color_dictionary('phosphate_position',metadat[sn]) 
					for si,sn in enumerate(sns) if si in nzs]
				for ii in range(len(xvals)):
					ax.bar(xvals[ii:ii+1],yvals[ii:ii+1],width=1,linewidth=0,
						color=colors[ii],alpha=alphas[ii],hatch=hatching[ii],zorder=3)
					if plot_type == 'relative_arrow':
						y0 = all_counts_start[combo][nzs[ii]]-1
						h = all_counts_change[combo][nzs[ii]]
						ax.add_patch(mpl.patches.FancyArrow(xvals[ii]+0.5,y0,0,h,zorder=5,
							linewidth=0,color=colors[ii],head_width=0.5,
							head_length={'pairs':0.02,'triples':0.1}[bridge_type],
							width=0.25,length_includes_head=False,alpha=0.3,
							path_effects=[path_effects.withStroke(linewidth=1,foreground='w')]))
				#---accumulate xlabels
				xlabels.extend([r'\texttt{'+metadat[sn]['ion_label']+withchar+metadat[sn]['ptdins_label']+\
					(' ('+metadat[sn]['composition_name']+')').rjust(12+2)+'}'
					for sn in [sns[i] for i in nzs]])
				xticks.extend(xvals+0.5)
				tb = ax.text((barparse*2+len(nzs))/2.,
					{'pairs':{'relative':0.4,'change':0.5,'relative_arrow':0.4}[plot_type],
					'triples':{'relative':2.2,'change':0.5,'relative_arrow':2.2}[plot_type]}[bridge_type],
					'\n'.join([c if c not in relabel else relabel[c] for c in combo]),
					ha='center',va='center',rotation=0,color='k',
					fontsize={'pairs':10,'triples':8}[bridge_type])
				ax.add_patch(mpl.patches.Rectangle((barparse-0.5,ylims[0]),
					len(nzs)+1,ylims[1]-ylims[0],angle=0.0,alpha=0.1,zorder=2,
					color=['k','w'][cc%2],linewidth=0))
				barparse += len(nzs)+1
			ax.set_xticklabels(xlabels,rotation=90,fontsize=fs['tiny'],ha='center')
			ax.set_xticks(xticks)
			ax.set_xlim(0.5,barparse-0.5)
			ax.set_ylim(ylims)
			ax.axhline(0,lw=2,c='k',zorder=4)
			ax.set_title({'pairs':'average lipid pairings','triples':'average lipid triples'}[bridge_type])
			ax.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='on')
			ax.tick_params(axis='y',which='both',left='off',right='off',labelbottom='off')
			ax.set_ylabel('observed %s versus random'%bridge_type,fontsize=fs['axlabel'])
			ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(
				lambda y, position : '%.f'%(100*y)+r'$\%$'))
			plt.savefig(figspot+'fig.lipid_pairings.'+plot_type+'.'+bridge_type+'.'+\
				'_'.join([re.findall('^membrane-(v[0-9]{3})',i)[0] for i in sns])+'.png',
				dpi=300,bbox_inches='tight')
			plt.close()	


		#---ptdins layout
		layout_direction = ['rows','columns','one_column'][1]
		compnames_ordered = calculations[compsign]['comparisons']
		layout = [
			['membrane-v509','membrane-v510','membrane-v511'],
			['membrane-v530','membrane-v531','membrane-v532'],
			[None,'membrane-v533','membrane-v534'],
			[None,None,None],
			]
		sns = list(set([i for j in layout for i in j if i != None]))

		#---plot layout
		ncols,nrows = len(layout),max([len(j) for j in layout])
		(ncols,nrows) = (ncols,nrows) if layout_direction == 'columns' else (nrows,ncols)
		no_plot = [u for v in [[([ii,jj] if layout_direction == 'rows' else [jj,ii]) 
			for jj,j in enumerate(i) if not j] for ii,i in enumerate(layout)] for u in v]
		if layout_direction == 'one_column': ncols,no_plot = 1,None
		fig,gs,axes = generate_axes(nrows=nrows,ncols=ncols,figsize=(18,8),no_plot=no_plot)

		legend_sns = ['membrane-v530','membrane-v534']
		colorset = colorscale('jet',count=len(postdat[sn][
			{'pairs':'combonames','triples':'combonames3'}[bridge_type]]))
		if bridge_type == 'pairs':
			clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[j] for j in range(5)+[7,8]]
		else:
			clrs = [brewer2mpl.get_map('Set1','qualitative',9).mpl_colors[j%9] for j in range(10)+range(10)+range(10)]
		
		#---loop over row/column in the layout
		for lnum,layset in enumerate(layout):
			#---loop over tiles in the row/column
			for sn in layset:
				if sn != None:
					rxmlook = extradat[sn]['rxmlook']
					reslist = extradat[sn]['reslist']		
					if layout_direction == 'one_column': sub_index = mp
					else: sub_index = layset.index(sn)
					ax = axes[lnum][sub_index] if layout_direction in ['rows','one_column'] \
						else axes[sub_index][lnum]
					mn = 0
					nframes = postdat[sn]['nframes']
					for ii,combo in enumerate(postdat[sn][
						{'pairs':'combonames','triples':'combonames3'}[bridge_type]]):
						side_lipids = list(set(rxmlook[mn]))
						composition_mn = [sum(rxmlook[mn]==i) for i in side_lipids]
						pairsums = [sum(rxmlook[mn]==reslist.index(c)) for c in combo]
						if not any(array(pairsums)==0):
							dat = postdat[sn]['%s%d'%(bridge_type,mn)][:,ii]/\
								mean(postdat[sn]['%s_random%d'%(bridge_type,mn)][:,ii])
							ps_per_frame = float(calculations[compsign]['timeslices'][sn]['time'].
								split('-')[-1])
							ax.plot(arange(nframes)*ps_per_frame/1000,dat,
								label='-'.join([c if c not in relabel else relabel[c] for c in combo]),
								color=clrs[ii])
							counts,edges = histogram(dat)
							mids = (edges[:-1]+edges[1:])/2.
					if sn in legend_sns:
						legend = ax.legend(loc='upper left',
							fontsize={'pairs':fs['tiny'],'triples':fs['tiny']-3}[bridge_type],
							bbox_to_anchor=(1.05,0,1.,1.))
						for label in legend.get_lines(): label.set_linewidth(2)
					ax.set_title(metadat[sn]['ptdins_label']+' with '+metadat[sn]['ion_label']+'\n'+
						metadat[sn]['composition_name'],fontsize=fs['title']-6)
					ax.set_xlabel('time (ns)',fontsize=fs['axlabel']-4)
		plt.subplots_adjust(top=0.85)
		gs.update(hspace=0.8)
		plt.suptitle('lipid pairings relative to random',fontsize=20)
		plt.savefig(figspot+'fig.lipid_pairings_trajectories.'+bridge_type+'.'+\
			'_'.join([re.findall('^membrane-(v[0-9]{3})',i)[0] for i in sns])+'.png',
			dpi=300,bbox_inches='tight')
		plt.close()

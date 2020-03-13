#!/usr/bin/env python

"""
Compute lipidwise areas.
"""

import seaborn as sb

@autoload(plotrun)
def load():
	data = plotload('lipid_area_analysis')
	sns = work.sns()

if __name__=='__main__':

	# violin with custom hatching looks best
	style = ['violin','bar'][0]
	#! custom hatches requires a hack to seaborn so violinplot accepts hatches
	custom_hatches = True

	#! custom top monolayer because the higher-z method failed on v536
	top_mono = lambda sn,abstractor: 1 if sn=='membrane-v536' and abstractor=='lipid_chol_com' else 0

	abstractor_name_path = ('calculation','lipid_areas2d','upstream','lipid_mesh',
		'upstream','lipid_abstractor','selector','loop')
	abstractor_names = delve(work.plotspec.specs,*abstractor_name_path)

	# custom hatching
	sb.set_style("white")
	mpl.rcParams['hatch.linewidth'] = 1.5
	mpl.rcParams['hatch.color'] = 'k'

	# loop over abstractor styles
	for abstractor_name in abstractor_names:

		for average_method in ['lipidwise','framewise','both']:

			data.set('lipid_areas2d',{
				('upstream','lipid_mesh','upstream','lipid_abstractor','selector'):abstractor_name})
			renamer_ptdins = lambda x: 'PtdIns' if x in ['PI2P','SAPI','P35P'] else x
			resnames_by_sn = [np.unique(data.this[sn]['resnames'][
				data.this[sn]['monolayer_indices']==top_mono(sn,abstractor_name)]) for sn in sns]
			resnames_all = [tuple([renamer_ptdins(j) for j in i]) for i in resnames_by_sn]
			if len(set(resnames_all))!=1: raise Exception('inconsistent residue names')
			resnames = list(set(resnames_all))[0]

			if style in ['violin','bar']:
				bar_formats = make_bar_formats(work.sns(),work)
				colors,hatches = [[bar_formats[sn][k] for sn in sns] for k in ['c','hatch']]
				palette = dict(zip(range(len(sns)),colors))
				labels = ['%s, %s'%(work.meta[sn]['ptdins_label'],work.meta[sn]['ion_label']) for sn in sns]
			else: raise Exception

			axes,fig = square_tiles(len(resnames),figsize=(10,10),favor_rows=True)
			for rnum,resname in enumerate(resnames):
				ax = axes[rnum]
				resnames_this = [data.this[sn]['resnames'] for sn in sns]
				# mean over axis 0 is framewise, axis 1 is lipidwise, otherwise you can also use reshape(-1)
				# note that averaging framewise gives you nframes number of readings, and shows the largest effects
				if average_method in ['lipidwise','framewise']:
					dists = [data.this[sn]['areas%d'%top_mono(sn,abstractor_name)].T[
						np.where(resnames_this[snum][data.this[sn][
							'monolayer_indices']==top_mono(sn,abstractor_name)]==resnames_by_sn[snum][rnum])
						].mean(axis={'lipidwise':1,'framewise':0}[average_method]) 
						for snum,sn in enumerate(sns)]
				else:
					# similar to the above but with no mean. instead all the data are in one bin
					dists = [data.this[sn]['areas%d'%top_mono(sn,abstractor_name)].T[
						np.where(resnames_this[snum][data.this[sn][
							'monolayer_indices']==top_mono(sn,abstractor_name)]==resnames_by_sn[snum][rnum])
						].reshape(-1) for snum,sn in enumerate(sns)]
				#! auditing
				if resname=='PtdIns' and abstractor_name=='lipid_com':
					status('resname: %s, abstractor_name: %s, means: %s'%(
						resname,abstractor_name,[np.mean(i) for i in dists]),tag='audit')
				"""
				>>> [data.this[sn]['areas%d'%top_mono(sn,abstractor_name)].T[np.where(resnames_this[snum][data.this[sn]['monolayer_indices']==top_mono(sn,abstractor_name)]==resnames_by_sn[snum][rnum])].mean(axis=0).mean() for snum,sn in enumerate(sns)][0]
				0.5155107911899196
				>>> [data.this[sn]['areas%d'%top_mono(sn,abstractor_name)].T[np.where(resnames_this[snum][data.this[sn]['monolayer_indices']==top_mono(sn,abstractor_name)]==resnames_by_sn[snum][rnum])].mean(axis=1).mean() for snum,sn in enumerate(sns)][0]
				0.51551079118991949
				>>> [data.this[sn]['areas%d'%top_mono(sn,abstractor_name)].T[np.where(resnames_this[snum][data.this[sn]['monolayer_indices']==top_mono(sn,abstractor_name)]==resnames_by_sn[snum][rnum])].reshape(-1).mean() for snum,sn in enumerate(sns)][0]
				0.5155107911899196

				how does the variation sum? what is the analytical form of it (i.e. do the moments sum in quadrature)? when the distributions are not perfectly normal, how does the sum deviate from the analytical form? if we can say that there is the most or "all" of the variation in the no-segmentation i.e. reshape(-1) case, then can we say that the two axis contribute to the variation differently? perhaps this explains why the distributions are wide for lipidwise. the lipids are all stuck in different microenvironments and their means are widely distributed, however the mean for each frame is much tighter, which means that far less variation in the areas is observed when all the lipids are averaged together. more variation for lipids for each frame is really putting one lipid over multiple frames "together" and in that case, because microenvironments, there is far more variation

				means are identical but the plots seem different in a very minor way. the dots on the violin means shift up and down coherently across simulations, so I suspect it is a plotting artefact

				[AUDIT] resname: PtdIns, abstractor_name: lipid_com, means: [0.68071527315876756, 0.71185134852641474, 0.6795538179720545, 0.69621366005163732, 0.71481869700710277, 0.70947392955744559, 0.70584131677206052, 0.66441957389293693]
				[STORE] searching pictures
				[STORE] saving picture to fig.lipid_areas2d.v2.png
				[AUDIT] resname: PtdIns, abstractor_name: lipid_com, means: [0.68071527315876756, 0.71185134852641463, 0.6795538179720545, 0.69621366005163732, 0.71481869700710277, 0.70947392955744559, 0.70584131677206052, 0.66441957389293693]
				[STORE] searching pictures
				[STORE] saving picture to fig.lipid_areas2d.v3.png
				[AUDIT] resname: PtdIns, abstractor_name: lipid_com, means: [0.68071527315876756, 0.71185134852641474, 0.67955381797205461, 0.69621366005163732, 0.71481869700710277, 0.70947392955744559, 0.70584131677206052, 0.66441957389293693]
				"""

				if style=='violin':
					#! added hatches to class _ViolinPlotter in seaborn/categorical.py
					if custom_hatches:
						vp = sb.violinplot(ax=ax,data=dists,palette=palette,
							labels=labels,linewidth=0,hatches=hatches)
						ii = 0
						for i in vp.get_children():
							if isinstance(i, mpl.collections.PolyCollection):
								i.set_hatch(hatches[ii])
								ii += 1
					else:
						#! check the distributions more carefully with a standard violin plot with whiskers
						vp = sb.violinplot(ax=ax,data=dists,palette=palette,
							labels=labels)
						#! if you want to compare distributions on the same axis: ax.set_ylim((0,2.0))
				elif style=='bar':
					bp = sb.boxplot(ax=ax,data=dists,palette=palette,linewidth=1,whis=1.0)
					for snum,bar in enumerate(bp.artists): bar.set_hatch(hatches[snum])
				else: raise Exception
				ax.set_title(work.vars.get('names',{}).get('short',{}).get(resname,resname))
				ax.set_xticklabels([])
				ax.set_ylabel('area per lipid $({nm}^{2})$')

			legendspec = []
			for snum,sn in enumerate(sns):
				legendspec.append(dict(name=labels[snum],
					patch=mpl.patches.Rectangle((-0.5,-0.5),1.5,1.5,fc=colors[snum],hatch=hatches[snum])))
			patches,labels = [list(j) for j in zip(*[(i['patch'],i['name']) for i in legendspec])]
			legend = ax.legend(patches,labels,loc='upper left',fontsize=14,bbox_to_anchor=(1.05,1.0),
				borderaxespad=0,handleheight=2.0)
			frame = legend.get_frame()
			frame.set_edgecolor('k')
			frame.set_facecolor('w')
			extras = [legend]
			meta = {'lipid_abstractor_selector':abstractor_name,'averaging_method':average_method}
			picturesave('fig.lipid_areas2d',work.plotdir,backup=False,version=True,meta=meta,extras=extras,form='pdf')

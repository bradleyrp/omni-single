#!/usr/bin/python -i

compsign,project_name = 'hydration','project-ptdins'
execfile('../../header.py')

from mpl_toolkits.axes_grid.inset_locator import inset_axes

from hydration import *

routine = [
	'plot_hydration',
	'plot_hydration_image',
	'plot_hydration_scatter',
	][:1]

debug = False

if 'datlist' not in globals() and not debug:

	computer(focus,compute_hydration,headerdat,simdict,get_slices=True)
	datlist = loader(focus,headerdat)

#---ptdins layout
layout_direction = ['rows','columns','one_column'][1]
compnames_ordered = calculations[compsign]['comparisons']
layout = [
	['membrane-v509','membrane-v510','membrane-v511'],
	['membrane-v530','membrane-v531','membrane-v532'],
	[None,'membrane-v533','membrane-v534'],
	[None,'membrane-v514','membrane-v515'],
	][1:-1]
sns = [i for j in layout for i in j if i != None]
#---plot layout
ncols,nrows = len(layout),max([len(j) for j in layout])
(ncols,nrows) = (ncols,nrows) if layout_direction == 'columns' else (nrows,ncols)
no_plot = [u for v in [[([ii,jj] if layout_direction == 'rows' else [jj,ii]) 
	for jj,j in enumerate(i) if not j] for ii,i in enumerate(layout)] for u in v]
if layout_direction == 'one_column': ncols,no_plot = 1,None


if 'postdat' not in globals():

	#---precompute histograms 
	postdat = {}
	#---loop over row/column in the layout
	for lnum,layset in enumerate(layout):
		#---loop over tiles in the row/column
		for sn in layset:
			if sn != None and sn not in postdat:
				
				dat = datlist[sn]
				cuts = unique([key.split(':')[0] for key in dat if ':' in key])
				cuts = [cuts[j] for j in argsort([float(a.split('-')[0]) for a in cuts])]
				colors = colorscale(name='jet',count=len(cuts))
				nframes = dat['nframes']
				imdat,means = [],[]
				for ci,cut in enumerate(cuts):
					data = concatenate([dat[cut+':%d'%fr] for fr in range(nframes)])
					edges = bins = (arange(10)-0.5)
					counts,edges = histogram(data,bins=edges,normed=False)
					mids = (edges[1:]+edges[:-1])/2.
					imdat.append(counts)
					means.append(mean(data))
				postdat[sn] = {
					'imdat':array(imdat)/dat['nframes'],
					'mids':array(mids),
					'means':array(means),
					'cuts':array([[float(j) for j in i.split('-')] for i in cuts]),
					}

if not debug and 'plot_hydration' in routine:

	fig,gs,axes = generate_axes(nrows=nrows,ncols=ncols,figsize=(12,12),no_plot=no_plot)

	#---loop over row/column in the layout
	for lnum,layset in enumerate(layout):
		#---loop over tiles in the row/column
		for sn in layset:
			if sn != None:
				if layout_direction == 'one_column': sub_index = mp
				else: sub_index = layset.index(sn)
				ax = axes[lnum][sub_index] if layout_direction in ['rows','one_column'] \
					else axes[sub_index][lnum]

				dat = datlist[sn]
				cuts = unique([key.split(':')[0] for key in dat if ':' in key])
				cuts = [cuts[j] for j in argsort([float(a.split('-')[0]) for a in cuts])]
				#---HACK!
				cuts = [cuts[0],cuts[1],cuts[2],cuts[-2],cuts[-1]]
				colors = colorscale(name='jet',count=len(cuts))#+1 is a hack
				nframes = dat['nframes']
				for ci,cut in enumerate(cuts):
					data = concatenate([dat[cut+':%d'%fr] for fr in range(nframes)])
					edges = bins = (arange(10)-0.5)
					counts,edges = histogram(data,bins=edges,normed=False)
					mids = (edges[1:]+edges[:-1])/2.
					ax.plot(mids,counts/dat['nframes'],'-',color=colors[ci],
						mew=0,lw=2,
						label=(r'$cutoff='+cut+'$'+'\n'+r'$mean=%.2f'%mean(data)+'$' 
							if ci==0 or ci==len(cuts)-1 else None))
				ax.set_title(metadat[sn]['ptdins_label']+' with '+metadat[sn]['ion_label']+\
					' ('+metadat[sn]['composition_name']+')')
				ax.set_xlabel('waters',fontsize=fs['axlabel'])
				ax.set_ylabel('observations (per frame)',fontsize=fs['axlabel'])
				h,l = ax.get_legend_handles_labels()
				ax.legend([h[ii] for ii,i in enumerate(l) if i!='None'],
					[l[ii] for ii,i in enumerate(l) if i!='None'],
					loc='center left',fontsize=fs['tiny'],bbox_to_anchor=(0.9,0.,1.,1.))
				add_axgrid(ax)
	gs.update(wspace=0.8,hspace=0.4)
	plt.savefig(figspot+'fig.hydration_distributions2.'+'.'.join([
		re.findall('^membrane-(v[0-9]+)$',s)[0] for s in sns])+'.png',dpi=300)
	plt.close()

if not debug and 'plot_hydration_scatter' in routine:

	fig,gs,axes = generate_axes(nrows=len(sns),ncols=1,figsize=(6,10),no_plot=no_plot)
	maxcount = max([postdat[sn]['imdat'].max() for sn in sns])
	
	for si,sn in enumerate(sns):
		ax = axes[si][0]
		npts = len(postdat[sn]['means'])
		xvals = concatenate(postdat[sn]['cuts'])
		yinds = concatenate(transpose((arange(npts),arange(npts))))
		yvals = postdat[sn]['means'][yinds]
		ax.plot(xvals,yvals,lw=1.5,color='k')
		ax.set_ylim(1.5,6.6)
		ax.set_ylabel('waters',fontsize=fs['axlabel'])
		ax.set_title(metadat[sn]['ptdins_label']+' with '+metadat[sn]['ion_label']+\
			' ('+metadat[sn]['composition_name']+')')
		if si == len(sns)-1: ax.set_xlabel('ion distance to lipid $\mathrm{(\AA)}$',fontsize=fs['axlabel'])
		ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
		ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		xticks = list(postdat[sn]['cuts'][:,0])+list(postdat[sn]['cuts'][-1,-1:])
		xticks = [0,2.2,4.4,5,6,7,8,9,10,20]
		ax.set_xticks(xticks)
		ax.set_xticklabels(xticks,rotation=0)
		add_axgrid(ax)
	
	gs.update(wspace=0.8,hspace=0.6)
	plt.savefig(figspot+'fig.hydration_distributions_means.'+'.'.join([
		re.findall('^membrane-(v[0-9]+)$',s)[0] for s in sns])+'.png',dpi=300)
	plt.close()

if not debug and 'plot_hydration_image' in routine:

	fig,gs,axes = generate_axes(nrows=nrows,ncols=ncols,figsize=(12,12),no_plot=no_plot)

	maxcount = max([postdat[sn]['imdat'].max() for sn in sns])

	#---loop over row/column in the layout
	for lnum,layset in enumerate(layout):
		#---loop over tiles in the row/column
		for sn in layset:
			if sn != None:
				if layout_direction == 'one_column': sub_index = mp
				else: sub_index = layset.index(sn)
				ax = axes[lnum][sub_index] if layout_direction in ['rows','one_column'] \
					else axes[sub_index][lnum]
				im = ax.imshow(array(postdat[sn]['imdat']).T,
					interpolation='nearest',origin='lower',
					vmax=maxcount)
				ax.set_title(metadat[sn]['ptdins_label']+' with '+metadat[sn]['ion_label']+\
					' ('+metadat[sn]['composition_name']+')')
				ax.set_ylabel('waters',fontsize=fs['axlabel'])
				ax.set_xlabel('ion distance to lipid $\mathrm{(\AA)}$',fontsize=fs['axlabel'])
				ax.set_xticks(arange(len(imdat)))
				ax.set_yticks(arange(len(imdat[0])))
				cutlabels = [r'\texttt{%.1f - %.1f}'%tuple([float(i) for i in j.split('-')]) for j in cuts]
				ax.set_xticklabels(cutlabels,rotation=90)
				ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
				ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
				if sn == 'membrane-v534':
					axins = inset_axes(ax,width="5%",height="100%",loc=3,
						bbox_to_anchor=(1.05,0.,1.,1.),
						bbox_transform=ax.transAxes,
						borderpad=0)
					cbar = plt.colorbar(im,cax=axins,orientation="vertical")
					axins.set_ylabel('observations (per frame)',fontsize=fs['axlabel']-2,
						rotation=270,labelpad=20)
	gs.update(wspace=0.8,hspace=0.4)
	plt.savefig(figspot+'fig.hydration_distributions_image.'+'.'.join([
		re.findall('^membrane-(v[0-9]+)$',s)[0] for s in sns])+'.png',dpi=300)
	plt.close()
	

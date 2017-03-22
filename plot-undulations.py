#!/usr/bin/env python

from codes import undulate
from copy import deepcopy

sns = work.vars['collections']['all']
colordict = work.vars['colordict']
routine = ['spectra','height'][:1]

def undulation_panel(ax,data,meta,keys=None,art=None,title=None,lims=None):
	"""
	Plot several undulation spectra on one panel.
	"""
	for sn in (keys if type(keys)==np.ndarray or not keys else data.keys()):
		dat = data[sn]['data']
		mesh = data[sn]['data']['mesh']
		vecs = data[sn]['data']['vecs']
		surf = np.mean(data[sn]['data']['mesh'],axis=0)
		#---! need to add fitting limit to meta file here
		uspec = undulate.calculate_undulations(surf,vecs,chop_last=True,perfect=True,lims=lims,raw=False)
		x,y = uspec['x'],uspec['y']
		color = colordict[sn]
		label = meta[sn]['label']
		label += '\n'+r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
		ax.plot(x,y,'o-',lw=2,markersize=5,markeredgewidth=0,c=color,label=label)
		ax.set_title(title)
	undulate.add_undulation_labels(ax,art=art)
	undulate.add_std_legend(ax,loc='upper right',art=art)
	undulate.add_axgrid(ax,art=art)
	
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

if 'spectra' in routine:

	art = {'fs':{'legend':12}}
	for wavevector_limit in [0.5,1.0,2.0]:
		status('plotting undulations with limit %f'%wavevector_limit,tag='plot')
		axes,fig = panelplot(layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},figsize=(5,5))
		undulation_panel(axes[0],data,work.meta,keys=np.array(sns),
			art=art,title='undulation spectra',lims=(0,wavevector_limit))
		meta = deepcopy(calc)
		meta['wavevector_limit'] = wavevector_limit
		descriptor = 'qlim.%.1f'%wavevector_limit
		picturesave('fig.%s.%s'%(plotname,descriptor),work,work.plotdir,backup=False,version=True,meta=meta)

if 'height' in routine:

	import render
	from render.wavevids import print_birdseye_snapshot_render

	framecounts = [data[sn]['data']['nframes'] for sn in data]
	nframes = min(framecounts)
	data_protein,calc = plotload('protein_abstractor',work)
	protein_pts_all = np.array([data_protein[sn]['data']['points_all'] for sn in sns])
	protein_pts = [[protein_pts_all[ii][fr][...,:2] for ii,sn in enumerate(sns)] for fr in range(nframes)]
	mvecs = np.array([data[sn]['data']['vecs'][:nframes] for sn in sns]).mean(axis=1)
	avgzs = [data[sn]['data']['mesh'].mean(axis=0).mean(axis=0) for sn in sns]
	avgzs = [i-i.mean() for i in avgzs]
	extrema = np.array(avgzs).min(),np.array(avgzs).max()
	extrema = [np.abs(extrema).max()*j for j in [-1,1]]
	nprots_list = [work.meta[sn]['nprots'] for sn in sns]
	titles = [work.meta[sn]['label'] for sn in sns]

	#---settings direct from pipeline-uvideos
	ppn = 4
	frameskipper = 2
	cmap_name = 'redblue'
	cmap_names = {'seismic':'seismic','redblue':'RdBu_r'}
	cmap_mpl_name = cmap_names[cmap_name]
	fs = {'axlabel':14,'title':20,'legend':14,'ticks':14,'tiny':10,'text_note':12,}
	print_review = False
	nrows,ncols = 2,2
	figsize = 8,8
	panelspecs = dict(layout={'out':{'grid':[1,1]},'ins':[{'grid':[nrows,ncols],'hspace':0.4}]},figsize=figsize)
	handle = 'wavevid-'+'.'.join([re.findall('simulation-(v[0-9]+).+',s)[0] for s in sns])+'.%s'%cmap_name

	print_birdseye_snapshot_render(
		avgzs,[p.mean(axis=0)[...,:2] for p in protein_pts_all],mvecs,nprots_list,
		handle='OUT',outdir=work.paths['post_plot_spot'],
		pbc_expand=1.0,smooth=1.,panelspecs=panelspecs,fn='fig.average_height',
		extrema=extrema,fs=fs,titles=titles,cmap_mpl_name=cmap_mpl_name,
		cbar_label_specs=dict(rotation=0,labelpad=-20,y=1.1),
		zlabel=r'$\mathrm{\langle z \rangle\,(nm)}$')

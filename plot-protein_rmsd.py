#!/usr/bin/env python

#---settings
plotname = 'protein_rmsd'
rmsd_bin_step = 0.1

#---colors
import matplotlib.patheffects as path_effects
colormap = lambda i,n,name='jet': mpl.cm.__dict__[name](float(i)/n)

#---load the upstream data
data,calc = plotload(plotname,work)

#---prepare an axis
axes,fig = panelplot(
	layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,2]}]},
	figsize=(12,8))

#---PLOT
counter,xpos,xlim_left = 0,[],0
max_rmsd = max([max(data[sn]['data']['rmsds']) for sn in data])
rmsd_bins = np.arange(0,max_rmsd*1.1,0.1)
for snum,sn in enumerate(data):
	color = colormap(snum,len(data))
	#---trajectory
	ax = axes[0]
	ts = np.array(data[sn]['data']['timeseries'])
	rmsds = data[sn]['data']['rmsds']
	ts -= ts.min()
	#---convert to ns
	ts = ts/1000.
	name_label = ' '.join(work.meta.get(sn,{}).get('name',sn).split('_'))
	ax.plot(ts,rmsds,label=name_label,lw=2,color=color)
	ax.set_xlabel(r'time (ns)')
	#---histograms
	ax = axes[1]
	counts,bins = np.histogram(rmsds,bins=rmsd_bins,normed=True)
	ax.fill_betweenx((bins[1:]+bins[:-1])/2.,counter,counter+counts,alpha=1.0,color=color,lw=0)
	max_rmsd = max([max_rmsd,max(rmsds)])
	#---compute drift
	m,b = np.polyfit(ts[len(rmsds)/10:],rmsds[len(rmsds)/10:],1)
	drift = m*(ts[-1]-ts[len(rmsds)/10])
	drift_start = ts[len(rmsds)/10]*m+b
	#---drift arrows
	ax.arrow(counter,drift_start,0,drift,
		head_width=max(counts)*0.1,head_length=max_rmsd*0.05,fc=color,ec='w',lw=1.5,
		path_effects=[path_effects.Stroke(linewidth=4,foreground=color),
		path_effects.Normal()],zorder=3)
	xpos.append(counter)
	counter += max(counts)*1.1
	if snum == 0: xlim_left = -1*max(counts)*0.1
	ax.set_xticks(xpos)
	ax.set_xticklabels([' '.join(work.meta.get(sn,{}).get('name',sn).split('_'))
		for sn in data],rotation=45,ha='right')
	ax.set_xlim(xlim_left,counter)
for ax in axes: 
	ax.set_ylim(0,max_rmsd*1.1)
	ax.set_ylabel(r'RMSD $(\textrm{\AA})$')
picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={})

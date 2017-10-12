#!/usr/bin/env python

routine = ['leaflet_area']

#---block: load the calculation data
if 'data' not in globals(): 
	data,calcs = plotload(plotname,work)
	data_prot,calcs_prot = plotload('protein_abstractor',work)
	sns = work.sns()

#---block: plot leaflet 3D areas
if 'leaflet_area' in routine:
		
	sn = sns[0]
	fig = plt.figure(figsize=(4,8))
	ax = plt.subplot(111)
	postdat = {}
	for snum,sn in enumerate(sns):
		top_mono = work.meta[sn].get('index_top_monolayer',0)
		areas = [data[sn]['data']['areas%d'%mn].sum(axis=1).mean() for mn in range(2)]
		areas_std = [data[sn]['data']['areas%d'%mn].sum(axis=1).std() for mn in range(2)]
		postdat[sn] = dict(mean=areas,std=areas_std)
		ax.bar(snum,areas[top_mono]-areas[[1,0][top_mono]],bottom=areas[[1,0][top_mono]],
			width=1.0,zorder=2)
		for mn in range(2):
			ax.errorbar([snum],[areas[mn]],yerr=[areas_std[mn]],alpha=1.0,lw=4.0,c='k',zorder=3)
	ax.set_xticks(range(len(sns)))
	ax.set_xticklabels([work.meta[sn]['label'] for sn in sns],rotation=90)
	picturesave('fig.%s'%plotname,work.plotdir,backup=False,version=True,meta={})
	plt.close()
	print(dict([(sn,'%.3f'%(np.ptp(postdat[sn]['mean'])/np.mean(postdat[sn]['mean']))) for sn in sns]))

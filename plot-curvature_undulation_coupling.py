#!/usr/bin/env python

from codes.curvature_coupling.curvature_coupling import InvestigateCurvature

avail = ['curvature_field_review','individual_reviews']
routine = work.plots[plotname].get('specs',{}).get('routine',avail)

if 'data' not in globals():
	try: data,calc = plotload(plotname)
	except:	data,calc = None,None
	#---get all upstream curvature sweeps
	ups = work.calc_meta.unroll_loops(work.calcs[plotname],return_stubs=True)[1]
	for up in ups: up['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for up in ups:
		#---temporarily set the plots
		work.plots[plotname]['calculation'] = {plotname:up['specs']}
		#---look up the right upstream data
		dat,cal = plotload(plotname)
		tag = up['specs']['design']
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])

if 'curvature_field_review' in routine:
	figsize = (10,10)
	cmap_name = 'RdBu_r'	
	layout = {'out':{'grid':[len(datas),1]},'ins':[{'grid':[1,len(work.sns())]} for d in datas]}
	axes,fig = panelplot(layout,figsize=figsize)
	cmax = max([np.abs(d[sn]['cf']).max() for d in datas.values() for sn in work.sns()])
	for dnum,(tag,dat) in enumerate(datas.items()):
		for snum,sn in enumerate(work.sns()):
			ax = axes[dnum][snum]
			ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
				vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name])
			error = datas[tag][sn]['bundle'][sn]['fun']
			ax.set_title('%s\n%s: %.4f'%(work.meta[sn].get('label',sn),tag,error))
	picturesave('fig.curvature_optimized_meta',work.plotdir,backup=False,version=True,
		meta={'maximum C_0':cmax,'sns':work.sns(),'designs':datas.keys()})

if 'individual_reviews' in routine:
	tag = datas.keys()[0]
	tag = 'v8'
	sn = work.sns()[0]
	figsize = (10,10)
	cmap_name = 'RdBu_r'	
	nthings = 2
	axes,fig = square_tiles(nthings,figsize)
	#---global values for consistency between reviews
	#---plot the mean curvature field
	ax = axes[0]
	cmax = np.abs(datas[tag][sn]['cf']).max()
	ax.imshow(datas[tag][sn]['cf'].T,origin='lower',interpolation='nearest',
		vmax=cmax,vmin=-1*cmax,cmap=mpl.cm.__dict__[cmap_name])
	#---plot the points for the gaussians
	ax = axes[1]
	pts = np.concatenate(datas[tag][sn]['drop_gaussians_points'].transpose(1,0,2)[::10])
	ax.scatter(*pts.T,c='k')
	#---! set axis limits here
	plt.show()

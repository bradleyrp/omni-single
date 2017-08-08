#!/usr/bin/env python

"""
INSPECT THE CURVATURE COUPLING RESULTS
"""

#---block: get the curvature investigation class
#---! this is pretty much redundant!
import codes.curvature_coupling.InvestigateCurvature
codes.curvature_coupling.InvestigateCurvature.plotload = plotload
codes.curvature_coupling.InvestigateCurvature.plotname = plotname
codes.curvature_coupling.InvestigateCurvature.work = work
#---alternate loader
from import_regular_membrane import plotloader_for_dextran
if 'ic' not in globals(): 
	ic = codes.curvature_coupling.InvestigateCurvature.InvestigateCurvature(plotloader=plotloader_for_dextran)

#---block: plot an error landscape
if False:
	landscapes = ic.generate_standard_manual_landscapes(isotropy=1.0)
	axes,fig = square_tiles(len(landscapes),figsize=(10,10))
	for snum,sn in enumerate(work.sns()):
		ax = axes[snum]
		ax.imshow(landscapes[sn].T,interpolation='nearest',origin='lower',cmap=mpl.cm.__dict__['gist_stern'])
	plt.show()

#---block: sweep isotropy
if False:
	theta = 90.0
	isotropy_sweep = ic.spec.get('isotropy',[1.0])
	postdat = {}
	for isotropy in isotropy_sweep:
		postdat[isotropy] = ic.generate_standard_manual_landscapes(isotropy=isotropy)
	axes,fig = panelplot(figsize=(7,10),layout={
		'out':{'grid':[len(work.sns()),1]},'ins':[{'grid':[1,len(isotropy_sweep)]} for sn in work.sns()]})
	#---determine a uniform range
	zmin,zmax = [[f([f(postdat[i][sn].reshape(-1)) for i in isotropy_sweep for sn in work.sns()])] 
		for f in [min,max]]
	for snum,sn in enumerate(work.sns()):
		for inum,isotropy in enumerate(isotropy_sweep):
			ax = axes[snum][inum]
			ax.imshow(postdat[isotropy][sn].T,
				interpolation='nearest',origin='lower',cmap=mpl.cm.__dict__['gist_stern'],vmin=zmin,vmax=zmax)
	plt.show()

#---block: check the undulation spectra
if True:
	art = {'fs':{'legend':10}}
	labels = dict([(i,i) for i in work.sns()])
	data,calc = plotloader_for_dextran('curvature_dextran')
	from codes.undulate_plot import undulation_panel
	axes,fig = panelplot(layout={'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]},figsize=(10,10))
	ax = axes[0]
	uspecs = undulation_panel(ax,data['undulations'],art=art,labels=labels)
	picturesave('fig.undulation_dev',work.plotdir,backup=False,version=True,meta={})
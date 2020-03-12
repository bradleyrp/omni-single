#!/usr/bin/env python

"""
Compute lipidwise areas.
"""

import scipy
import scipy.stats
import seaborn as sb

@autoload(plotrun)
def load():
	data,calc = plotload('lipid_area_converge')
	sns = work.sns()

if __name__=='__main__':
	sns = work.sns()
	meta = {}
	legends = []
	dpi = 300
	img_form = 'pdf'
	axes,fig = square_tiles(1,
		figsize=(1,1),hspace=0.4,wspace=0.4)
	ax = axes[0]
	drifts = {}
	drifts_pct = {}
	for sn in sns:
		series = np.product(data[sn]['data']['vecs'][:,:2],axis=1)
		start,end = calc['extras'][sn]['start'],calc['extras'][sn]['end']
		times = np.linspace(start,end,len(series))
		ax.plot(series)
		regress = scipy.stats.linregress(x=times,y=series)
		slope, intercept, r_value, p_value, std_err = regress
		drifts[sn] = (intercept+slope*end)-(intercept+slope*start)
		drifts_pct[sn] = drifts[sn]/series.mean()
	picturesave('fig.lipid_areas.convergence.all',directory=work.plotdir,
		version=True,meta=meta,extras=legends,dpi=dpi,form=img_form)
	print('[STATUS] maximum drift is %0.3f'%max(drifts_pct.values()))
	long_ones = ['membrane-v565','membrane-v563']
	max_drift_short = max([j for i,j in drifts_pct.items() if i not in long_ones])
	print('[STATUS] maximum drift excluding 500ns simulations is %0.3f'%max_drift_short)

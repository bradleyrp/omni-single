#!/usr/bin/env python

"""
Extract curvature undulation coupling data for the dextran project for further analysis.
"""

if 'data' not in globals():
	whittle_calc = {'curvature_undulation_coupling_dex':{
		'design':'v3','upstream':{
		'import_readymade_meso_v1_nanogel':{},'import_readymade_meso_v1_membrane':{}}}}
	data,calc = plotload(plotname,whittle_calc=whittle_calc)
	sns = work.sns()

#---plot the ratio of the energies
fig = plt.figure()
ax = plt.subplot(111)
for snum,sn in enumerate(sns):
	dat = data[sn]['data']
	ax.scatter(dat['qs'],dat['ratios'],c='kb'[snum],s=1)
ax.set_xscale('log')
ax.set_yscale('log')
picturesave('fig.TEST',work.plotdir)

#---save the data to two-column format
np.savetxt(os.path.join(work.plotdir,'test.dat'),np.transpose((dat['qs'],dat['ratios'])))

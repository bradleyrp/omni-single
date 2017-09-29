#!/usr/bin/env python

"""
Extract curvature undulation coupling data for the dextran project for further analysis.
"""

#---choose what to plot
routine = [
	'entropy_via_curvature_undulation',
	'entropy_via_undulations',
	'spectra_comparison',
	][1:2]
#---high cutoff for undulations since it lacks this setting
high_cutoff_undulation = 0.1

from base.tools import gopher
from hypothesis import sweeper
from codes.undulate import calculate_undulations
from codes.undulate_plot import undulation_panel
import copy

#---share writes
os.umask(0o002)

#---fetch the entropy function
if 'entropy_function' not in globals():
	entropy_loader = {
		'module':'codes.dextran_entropy_calculation',
		'function':'undulations_to_entropy'}
	entropy_function = gopher(entropy_loader)

if 'entropy_via_curvature_undulation' in routine:

	#---load all possible upstream curvature undulation coupling data and get the entropy
	if 'sweep' not in globals():
		#---get all of the data at once. try plotload with whittle_calc if it's too large
		sweep = work.collect_upstream_calculations_over_loop('curvature_undulation_coupling_dex')

	entropy_survey = {}
	#---compute entropies for each simulation, coupling hypothesis
	for tag in sweep['datas']:
		for sn in sweep['datas'][tag]:
			dat = sweep['datas'][tag][sn]
			high_cutoff = sweep['calcs'][tag][sn]['calcs']['specs']['fitting']['high_cutoff']
			kappa,sigma,constant = dat['x'][:3]
			result = entropy_function(dat['qs'],dat['ratios'],
				high_cutoff=high_cutoff,kappa=kappa,sigma=sigma)
			entropy_survey[(tag,sn)] = result

if 'entropy_via_undulations' in routine:

	#---load all possible upstream curvature undulation coupling data and get the entropy
	if 'data' not in globals():
		#---use plotload to get membrane shapes
		data,calc = plotload('import_readymade_meso_v1_membrane')
		sns = work.sns()

	#---inject undulation fits to the incoming membrane data for upcoming entropy calculation
	for sn in data:
		dat = data[sn]['data']
		surf = dat['mesh'].mean(axis=0)
		vecs = dat['vecs']
		#---choose a specific fitting method
		#---! talk this over with Ravi, others to make sure correct theory
		fit_style = 'band,perfect,curvefit-crossover'
		residual_form = 'log'
		midplane_method = 'average'
		fit_tension = True
		lims = (0.0,high_cutoff_undulation)
		#---compute undulations
		result = calculate_undulations(surf,vecs,chop_last=False,residual_form=residual_form,
			midplane_method=midplane_method,fit_style=fit_style,fit_tension=True,
			lims=lims,perfect=True,raw=False)
		#---save this for later
		dat['undulation_postprocessing'] = result

	#---use height-height undulations to calculate the entropy
	entropy_survey_undulations = {}
	for sn in data:
		dat = data[sn]['data']
		post = dat['undulation_postprocessing']
		result = entropy_function(post['q_binned'],post['energy_binned'],high_cutoff=high_cutoff_undulation,
			kappa=post['kappa'],sigma=post['sigma'])
		entropy_survey_undulations[sn] = result

if 'spectra_comparison' in routine:

	"""
	notes on the fitting method: 
		1. basic (fits the exponent, no difference in residual forms)
			fits are good but kappa is impossible to interpret with flexible exponents
			exponents for flat,average,average_normal are 4.294,3.906,3.742
				per usual average is best (not sure why average_normal is not best)
		2. simple
			works great but is obviously a bit crude
		3. curvefit 
			note that we removed the manual log transformation and now just fit directly in log space
				including sigma, with no algebra and no tricks with the crossover
			also note that the linear fit gives the predictable response 
				whereby it always overshoots and downweights errors lower in the scale
				hence we will not consider the linear part anymore
			you get different results depending on the limits
		4. fit is an advanced scipy.optimize method that, ironically, is getting stuck a lot
			so we are tabling it for now
	remaining questions: 
		1. should we enforce a strict cutoff and fit different parts of the spectra
	"""

	from codes.undulate_plot import add_undulation_labels,add_axgrid,add_std_legend

	data,calc = plotload('import_readymade_meso_v1_membrane',work)
	data_average_normal,_ = plotload('undulations_average_normal',work)

	def hqhq(q_raw,kappa,sigma,area,exponent=4.0):
		return 1.0/(area/2.0*(kappa*q_raw**(exponent)+sigma*q_raw**2))

	sns = work.sns()
	sn = sns[0]
	art = {'fs':{'legend':8}}

	lims = (0.0,high_cutoff_undulation)
	plotspecs = sweeper(**{
		'fit_style':['band,perfect,curvefit','band,perfect,curvefit-crossover'][:1],
		'midplane_method':['flat','average','average_normal'],
		'lims':[(0.0,high_cutoff_undulation),(0.04,high_cutoff_undulation)],
		'residual_form':['log','linear'][:1],
		'fit_tension':[False,True][:]})

	axes,fig = square_tiles(len(plotspecs),figsize=18,favor_rows=True,wspace=0.4,hspace=0.1)
	for pnum,plotspec in enumerate(plotspecs):
		ax = axes[pnum]

		fit_style = plotspec['fit_style']
		residual_form = plotspec['residual_form']
		midplane_method = plotspec['midplane_method']
		fit_tension = plotspec['fit_tension']
		lims = plotspec['lims']

		if midplane_method=='average_normal': 
			custom_heights = data_average_normal[sn]['data']['average_normal_heights']
		else: custom_heights = None
		dat = data[sn]['data']
		mesh = dat['mesh']
		vecs = dat['vecs']
		surf = np.mean(dat['mesh'],axis=0)

		uspec = calculate_undulations(surf,vecs,
			fit_style=fit_style,custom_heights=custom_heights,lims=lims,
			midplane_method=midplane_method,residual_form=residual_form,fit_tension=fit_tension)

		label = 'structure: %s'%re.sub('_',' ',midplane_method)
		label += '\n residuals: %s'%residual_form
		label += '\n method: \n%s'%fit_style
		label += '\n'+r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
		if uspec['sigma']!=0.0:
			label += '\n'+r'$\mathrm{\sigma='+('%.3f'%uspec['sigma'])+'\:{k}_{B} T {nm}^{-2}}$'
		colors = ['b','k','r']

		q_binned,energy_binned = uspec['q_binned'],uspec['energy_binned']
		ax.plot(q_binned,energy_binned,'.',lw=0,markersize=10,markeredgewidth=0,
			c=colors[0],label=None,alpha=0.2)
		q_fit,energy_fit = np.transpose(uspec['points'])
		ax.plot(q_fit,energy_fit,'.',lw=0,markersize=4,markeredgewidth=0,
			c=colors[1],label=label,alpha=1.,zorder=4)
		#---alternate exponent
		if 'linear_fit_in_log' in uspec:
			exponent = uspec['linear_fit_in_log']['c0']
			status('alternate exponent %.3f'%exponent,tag='note')
			ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
				#---distinctive green color for alternate exponents
				exponent=-1.0*exponent,area=uspec['area']),lw=3,zorder=5,c='g')
		elif 'crossover' in uspec:
			ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=0.0,
				area=uspec['area']),lw=3,zorder=5,c='g')
			ax.plot(q_fit,hqhq(q_fit,kappa=0,sigma=uspec['sigma'],
				area=uspec['area']),lw=3,zorder=5,c='g')
			ax.axvline(uspec['crossover'],c='k',lw=1.0)
		#---standard exponent
		else:
			ax.plot(q_fit,hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],
				area=uspec['area']),lw=1,zorder=5,c=colors[2])
			#---inset axis shows the residuals
			#---! note that there is a big whereby the log tick marks still show up sometimes
			axins = inset_axes(ax,width="30%",height="30%",loc=3)
			diffs = (energy_fit-
				hqhq(q_fit,kappa=uspec['kappa'],sigma=uspec['sigma'],area=uspec['area']))
			axins.set_title('residuals',fontsize=6)
			if residual_form=='log': 
				axins.plot(10**diffs,'.-',lw=0.5,ms=1,color='k')
				axins.set_xscale('log')
				axins.set_yscale('log')
				axins.axhline(1.0,c='k',lw=0.5,alpha=0.5)
			else: 
				axins.axhline(0.0,c='k',lw=0.5,alpha=0.5)
				axins.plot(diffs,'.-',lw=0.5,ms=1,color='k')
			axins.set_xticks([])
			axins.set_yticks([])
			axins.set_xticklabels([])
			axins.set_yticklabels([])
			axins.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
			axins.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
		add_undulation_labels(ax,art=art)
		add_std_legend(ax,loc='upper right',art=art)
		add_axgrid(ax,art=art)

	picturesave('fig.undulation_survey',work.plotdir,backup=False,version=True,meta={})

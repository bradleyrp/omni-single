#!/usr/bin/env python

"""
Extract curvature undulation coupling data for the dextran project for further analysis.
"""

#---! RYAN MOVE THIS
def collect_upstream_calculations_over_loop(calcname):
	"""
	Some plotting and analysis benefits from checking all calculations in an upstream loop (which is 
	contrary to the original design of )
	!!! This script is a candidate for inclusion in omnicalc.py
	"""
	global plotname,work
	#---! move to a separate function
	plotspecs = work.plots.get(plotname,work.calcs.get(plotname,{})).get('specs',{})
	if not calcname: calcname = plotspecs.get('calcname',plotname)
	#---load the canonical upstream data that would be the focus of a plot in standard omnicalc
	#---! load the upstream data according to the plot. note that this may fail in a loop hence needs DEV!
	try: data,calc = plotload(plotname)
	except:
		data,calc = None,None
		status('failed to load a single upstream calculation however this plot script has requested '
			'all of them so we will continue with a warning. if you have downstream problems consider '
			'adding a specific entry to plots to specify which item in an upstream loop you want',
			tag='warning')
	#---in case there is no plot entry in the metadata we copy it
	if plotname not in work.plots: work.plots[plotname] = copy.deepcopy(work.calcs[calcname])
	#---load other upstream data
	#---get all upstream curvature sweeps
	upstreams,upstreams_stubs = work.calc_meta.unroll_loops(work.calcs[calcname],return_stubs=True)
	datas,calcs = {},{}
	#---loop over upstream calculations and load each specifically, using plotload with whittle_calc
	for unum,upstream in enumerate(upstreams_stubs):
		#---use the whittle option to select a particular calculation
		dat,cal = plotload(calcname,whittle_calc={calcname:upstream['specs']})
		tag = upstreams_stubs[unum]['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
	#---singluar means the typical "focus" of the upstream calculation, plural is everything else
	return dict(datas=datas,calcs=calcs,data=data,calc=calc)

#---settings
routine = ['entropy_via_curvature_undulation','entropy_via_undulations'][-1:]
high_cutoff_undulation = 0.2

from base.tools import gopher
from codes.undulate import calculate_undulations

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
		sweep = collect_upstream_calculations_over_loop('curvature_undulation_coupling_dex')

	entropy_survey = {}
	#---compute entropies for each simulation, coupling hypothesis
	for tag in sweep['datas']:
		for sn in sweep['datas'][tag]:
			dat = sweep['datas'][tag][sn]['data']
			#---! get the high cutoff!
			#---! need to pipe kappa,sigma in through kwargs below before continuing!
			result = entropy_function(dat['qs'],dat['ratios'],high_cutoff=high_cutoff)
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
		#---! check if we really want perfect collapser
		result = calculate_undulations(surf,vecs,chop_last=False,
			lims=(0,high_cutoff_undulation),perfect=True,raw=False)
		#---save this for later
		dat['undulation_postprocessing'] = result

	#---use height-height undulations to calculate the entropy
	entropy_survey_undulations = {}
	for sn in data:
		dat = data[sn]['data']
		post = dat['undulation_postprocessing']
		result = entropy_function(post['x'],post['y'],high_cutoff=high_cutoff_undulation,
			kappa=post['kappa'],sigma=post['sigma'])
		entropy_survey_undulations[sn] = result

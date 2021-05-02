#!/usr/bin/env python

"""
Plot all of the upstream loops for the curvature undulation coupling analysis.

"""

from codes.curvature_coupling.curvature_coupling_plots import individual_reviews_plotter

#---autoplot settings
plotrun.routine = ['individual_reviews']
global seepspace

# moved from plot-curvature_undulation_coupling_pixel.py on 2021.04.17
def collect_upstream_calculations_over_loop():
	"""
	Some plotting and analysis benefits from checking all calculations in an upstream loop (which is 
	contrary to the original design of )
	!!! This script is a candidate for inclusion in omnicalc.py
	"""
	global plotname,work
	#---! move to a separate function
	plotspecs = work.metadata.plots.get(plotname,work.metadata.calculations.get(plotname,{})).get('specs',{})
	calcname = plotspecs.get('calcname',plotname)
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
	if plotname not in work.metadata.plots: work.metadata.plots[plotname] = copy.deepcopy(work.metadata.calculations[calcname])
	#---load other upstream data
	#---get all upstream curvature sweeps
	upstreams,upstreams_stubs = work.calcs.unroll_loops(work.metadata.calculations[calcname],return_stubs=True)
	#for u in upstreams_stubs: u['specs'].pop('upstream',None)
	datas,calcs = {},{}
	for unum,upstream in enumerate(upstreams_stubs):
		#---use the whittle option to select a particular calculation
		dat,cal = plotload(calcname,whittle_calc={calcname:upstream['specs']})
		tag = upstreams_stubs[unum]['specs']['design']
		if type(tag)==dict: tag = 'v%d'%unum
		datas[tag] = dict([(sn,dat[calcname][sn]['data']) for sn in work.sns()])
		calcs[tag] = dict([(sn,cal) for sn in work.sns()])
	#---singluar means the typical "focus" of the upstream calculation, plural is everything else
	return dict(datas=datas,calcs=calcs,data=data,calc=calc)

@autoload(plotrun)
def loader():
	"""Load data."""
	#---only load once
	if 'data' not in globals():
		#---begin loading sequence
		if plotname not in work.metadata.plots: raise Exception('add %s to the plots metadata'%plotname)
		plotspecs = work.metadata.plots[plotname].get('specs',{})
		calcname = plotspecs.get('calcname',plotname)
		#---new method for getting all upstream calculations in the loop
		# moved to the function above
		combodat = collect_upstream_calculations_over_loop()
		#---use the whittle option to select a particular calculation
		data,datas,calcs = combodat['data'],combodat['datas'],combodat['calcs']
		protein_abstractor_name = plotspecs.get('protein_abstractor_name','protein_abstractor')
		undulations_name = plotspecs.get('undulations_name','undulations')
		#---check for alternate loaders
		alt_loader = plotspecs.get('loader',None)
		if alt_loader:
			from base.tools import gopher
			postdat = gopher(alt_loader)(data=data)
		#---compile the necessary (default) data into a dictionary
		else:
			postdat = dict([(sn,dict(
				vecs=data[undulations_name][sn]['data']['vecs'].mean(axis=0),
				points_protein=data[protein_abstractor_name][sn]['data']['points_all']))
				for sn in work.sns()])
		#---end loading sequence
		#---export to globals
		global seepspace
		seepspace = 'data,datas,postdat,undulations_name,calcs,protein_abstractor_name'.split(',')
		for key in seepspace: globals()[key] = locals()[key]

@autoplot(plotrun)
def individual_reviews():
	"""
	Plot a few versions of the main figure.
	This is the standard plot. See other functions below for spinoffs.
	"""
	plotspec = {
		'coupling_review':{
			'viewnames':['average_height','average_height_pbc','neighborhood_static',
				'neighborhood_dynamic','average_field','example_field','example_field_pbc',
				'spectrum','spectrum_zoom']},}
	global seepspace
	seep = dict([(key,globals()[key]) for key in seepspace])
	for out_fn,details in plotspec.items(): 
		individual_reviews_plotter(out_fn=out_fn,seep=seep,**details)

@autoplot(plotrun)
def individual_reviews_ocean():
	"""
	Custom version of individual_reviews for the ocean project.
	Not included in the default plots.
	"""
	plotspec = {
		'coupling_review.simple':{
			'viewnames':['average_height','example_field_no_neighborhood'],'figsize':(6,6),'horizontal':True},
		'coupling_review.center_debug':{
			'viewnames':['neighborhood_static','neighborhood_dynamic',
				'average_height','average_height_pbc','average_field','example_field','example_field_pbc',
				'average_height_center','average_field_center','average_field_center_smear',
				'curvature_field_center'],
				'figsize':(16,16)},
		#---! this fails on the helix-0 tests because they have multiple proteins
		'coupling_review.simple_centered':{
			'viewnames':['average_height_center','curvature_field_center'],
			'figsize':(8,8),'horizontal':True,'wspace':0.7},}
	#---turn some off when developing
	for i in ['coupling_review.center_debug']: plotspec.pop(i)
	global seepspace
	seep = dict([(key,globals()[key]) for key in seepspace])
	for out_fn,details in plotspec.items(): 
		individual_reviews_plotter(out_fn=out_fn,seep=seep,**details)

@autoplot(plotrun)
def compare_curvature_estimates():
	"""Find the best curvature estimates for each method."""
	sns = work.sns()
	comp = dict([(sn,(datas.keys()[i],datas.values()[i][sn]['bundle'][sn]['fun'])) 
		for i in [np.argmin([datas[tag][sn]['bundle'][sn]['fun'] for tag in datas])] for sn in sns])
	print(comp)

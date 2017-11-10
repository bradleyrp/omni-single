#!/usr/bin/env python

"""
Plot bilayer areas for the ptdins project.
"""

"""
A brief primer on the new-style plot scripts.
Run with `make plot undulations` for interactive mode (plots everything).
Or use `make plot undulations <args>` to name the plot functions you want.
Check name main for prototyping code at the bottom.
Modify plot_super.routine (a list) to name functions to run, otherwise everything.
To plot everything on replot after setting plot_super.routine, reset it to None.
"""

from codes.undulate_plot import undulation_panel

@autoload(plot_super)
def load():
	"""Load the data."""
	#---condition for loading the data
	if 'data' not in globals(): 
		#---load once
		global data,calc
		data,calc = plotload(plotname,work)

@autoplot(plot_super)
def plot_undulation_spectra():
	"""
	Plot the undulation spectra.
	"""
	plotspecs = work.plots[plotname].get('specs',{})
	wavevector_limits = plotspecs.get('wavevector_limits',[0.5,1.0,2.0])
	#---settings
	art = {'fs':{'legend':12}}
	layout = {'out':{'grid':[1,1]},'ins':[{'grid':[1,1]}]}
	figsize=(5,5)
	#---sweep over plots
	plots_this = sweeper(**{'layout':layout,'wavevector_limit':wavevector_limits})
	#---note that the default midplane method for undulation spectra is "average"
	#---...however you can also use "flat". note that curvature undulation coupling uses default "flat"
	midplane_method = plotspecs.get('midplane_method','average')
	for pnum,plotspec in enumerate(plots_this):
		status('plotting layout',i=pnum,tag='plot')
		wavevector_limit = plotspec['wavevector_limit']
		axes,fig = panelplot(layout,figsize=figsize)
		for aa,ax in enumerate(axes):
			#---panel settings
			title = 'undulation spectra'
			keys = None
		#---metadata and file tag
		meta = deepcopy(calc)
		meta.update(**plotspec)
		#---! ERROR WHY IS META NOT USED?
		tag = 'qlim.%.1f'%wavevector_limit
		picturesave('figT',work.plotdir,backup=False,version=True,meta=plotspec)

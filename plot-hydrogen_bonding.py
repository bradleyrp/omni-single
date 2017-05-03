#!/usr/bin/env python

"""
PLOT hydrogen bonding
"""

#---block: what to plot
routine = ['spectra','height'][:]
sns = work.sns()

#---block: load the calculation data
if 'data' not in globals(): 
	data,calc = plotload(plotname,work)

#---all simulations in one bar plot
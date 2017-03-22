#!/usr/bin/python

"""
Common figures specifications.
"""

def figplacer(sn,figplace):

	"""
	Return the outer and inner axes indices for a particular simulation given a figplace.
	"""

	rown = next(row for row,i in enumerate(figplace) if sn in i)
	coln = next(i.index(sn) for row,i in enumerate(figplace) if sn in i)
	return rown,coln
	
def replot(): 

	"""
	Rerun the script with a single command.
	This is very useful for tuning plots, or in cases where the plot apparatus is used to write new 
	calculations, for rapidly developing new calculations.
	"""

	import sys
	print("[STATUS] replotting via %s"%sys.argv[0])
	execfile(sys.argv[0])

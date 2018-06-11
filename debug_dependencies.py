#!/usr/bin/env python

"""
"""

import MDAnalysis

if 'data' not in globals():
	plotname = 'lipid_abstractor'
	data,calc = plotload('lipid_abstractor')
	sns = work.sns()
	# example trajectory
	sn = sns[0]
	struct,traj = [os.path.join(work.postdir,calc['extras'][sn]['slice_path']+'.%s'%i) for i in ['gro','xtc']]
	uni = MDAnalysis.Universe(struct,traj)
	sel = uni.select_atoms('all')
	#! note that I was expecting an older MDAnalysis on dark to have different residue behavior
	#! ... however it is more likely that I was not making use of the uniqueness features of sel.residues
	#! ... in lipid_abstractor recently, which is why it was failing on some actinlink calculations with the
	#! ... headless RDF
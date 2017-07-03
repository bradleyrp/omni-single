#!/usr/bin/env python

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

vmd_render_good_side_head_angle = """

draw delete all
set me [atomselect top "all"]
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
$me moveby [vecscale -1.0 [expr [$me_c11 get {x y z}]]]
display resetview
mouse stoprotation
rotate x to -90
if { $reverse_view == 1 } { rotate y by 90 } else { rotate y by -90 }
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
set vec_ring [vecsub [expr [$me_c14 get {x y z}]] [expr [$me_c11 get {x y z}]]]
set vec_up {0.0 0.0 1.0}
set vec_look [veccross $vec_ring $vec_up]
set relevant_angle [expr [tcl::mathfunc::acos [vecdot [vecnorm $vec_ring] $vec_up]]/3.1415926*180.0]
set tmat [transvecinv $vec_look]
$me move $tmat
set me_c11 [atomselect top "resname PIP2RESNAME and resid $my_index and name C11"]
set me_c14 [atomselect top "resname PIP2RESNAME and resid $my_index and name C14"]
# REMOVED draw cylinder [expr [$me_c11 get {x y z}]] [expr [$me_c14 get {x y z}]] radius 0.1
"""

vmd_render_good_side_bond = """

# project a vector to a certain plane
# Usage: projection $v1 x ==> project vector v1 to y-z plane
proc projection { vec axis } {
if {$axis == "x"} {
set vec "0 [lindex $vec 1] [lindex $vec 2]"
} elseif {$axis == "y"} {
set vec "[lindex $vec 0] 0 [lindex $vec 2]"
} else {
set vec "[lindex $vec 0] [lindex $vec 1] 0"
}
return $vec
}

#---move to partner_2_bead, which will be the axis
set me [atomselect top "all"]
$me moveby [vecscale -1.0 [expr [$partner_2_bead get {x y z}]]]

#---reset the view
display resetview
mouse stoprotation
rotate x to -90

#---project the vector onto the plane and perform the rotation
set vec_main [vecsub [expr [$partner_2_bead get {x y z}]] [expr [$partner_1_bead get {x y z}]]]
set vec_main_flat [projection $vec_main z]
set vec_up {0.0 0.0 1.0}
set vec_look [veccross $vec_main_flat $vec_up]
set tmat [transvecinv $vec_look]
$me move $tmat

#---tie the whole thing together because we just rotated into the wrong direction
rotate y by 90
"""

def exemplar_lipid(theta,phi,x,y,nframes,nlipids,rank=0):
	"""Select a single lipid and frame from a set of statistics."""
	resid_rel,frame = np.unravel_index(np.argsort(
		np.sum((np.array((x,y)).T-np.array([theta,phi]))**2,axis=1))[rank],(nlipids,nframes))
	return frame,resid_rel

def snapshot_namer(**kwargs):
	"""Name the snapshots in the plot folder."""
	canon_order = ['tag','sn','pos']
	mark_out,mark_in = '.','_'
	suffix = mark_out.join(
		'%s%s%s'%(str(key),mark_in,str(val))
		for key,val in [(k,kwargs[k]) for k in canon_order])
	return 'fig.head_angle.snap.'+suffix

def angle_image_format(ax,sn=None,fsbase=20):
	ax.tick_params(axis='y',which='both',left='off',right='off',labelleft='on')
	ax.tick_params(axis='x',which='both',top='off',bottom='off',labelbottom='on')
	plt.setp(ax.xaxis.get_majorticklabels(),rotation=90,fontsize=fsbase-2)
	plt.setp(ax.yaxis.get_majorticklabels(),rotation=0,fontsize=fsbase-2)
	ax.axhline(0,lw=1.5,c='w',zorder=4)
	ax.axvline(0,lw=1.5,c='w',zorder=4)
	if sn:
		ax_side = ax.twinx()
		ax_side.set_yticks([])
		ax_side.set_ylabel(work.meta[sn]['composition_name'],fontsize=fsbase,rotation=270,labelpad=25)


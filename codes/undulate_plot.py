#!/usr/bin/env python

import numpy as np
from . undulate import calculate_undulations

def undulation_panel(ax,data,keys=None,art=None,title=None,lims=None,
	colors=None,labels=None,show_fit=True,midplane_method=None,custom_heights=None):
	"""
	Plot several undulation spectra on one panel.
	"""
	#---function for the q4-scaline
	func_q4 = lambda q,kappa: kappa*q**-4
	uspecs = {}
	sns = keys if keys else data.keys()
	for sn in sns:
		dat = data[sn]['data']
		mesh = data[sn]['data']['mesh']
		vecs = data[sn]['data']['vecs']
		surf = np.mean(data[sn]['data']['mesh'],axis=0)
		#---kernel of this plot: calculate the spectra here
		uspec = calculate_undulations(surf,vecs,chop_last=True,custom_heights=custom_heights,
			perfect=True,lims=lims,raw=False,midplane_method=midplane_method)
		uspecs[sn] = uspec
		x,y = uspec['q_binned'],uspec['energy_binned']
		label = labels[sn] if labels else sn
		label += '\n'+r'$\mathrm{\kappa='+('%.1f'%uspec['kappa'])+'\:k_BT}$'
		#---colors should be a dict over the keys
		color = colors[sn] if colors else None
		ax.plot(x,y,'o-',lw=2,markersize=5,markeredgewidth=0,c=color,label=label)
		ax.set_title(title)
		if show_fit: ax.plot(x,func_q4(x,kappa=uspec['kappa'])/np.product(vecs.mean(axis=0)[:2])/100.*2.,zorder=5,c='k',lw=1)
	add_undulation_labels(ax,art=art)
	add_std_legend(ax,loc='upper right',art=art)
	add_axgrid(ax,art=art)
	return uspecs

def add_undulation_labels(ax,art=None):
	"""
	Label an axes for undulation spectra.
	"""
	ax.set_ylabel(r'$\left\langle h_{q}h_{-q}\right\rangle \left(\mathrm{nm}^{2}\right)$')
	ax.set_xlabel(r'${\left|\mathbf{q}\right|}(\mathrm{{nm}^{-1}})$')
	ax.set_xscale('log')
	ax.set_yscale('log')

def add_axgrid(ax,art=None):
	"""
	Standard plot function for adding gridlines.
	"""
	ax.grid(True,linestyle='-',zorder=0,alpha=0.35)
	ax.set_axisbelow(True)

def add_std_legend(ax,loc='upper right',lw=2.0,title=None,art=None):
	"""
	Add a legend.
	"""
	h,l = ax.get_legend_handles_labels()
	leginds = range(len(h))
	legend = ax.legend([h[i] for i in leginds],[l[i] for i in leginds],
		loc=loc,fontsize=art['fs']['legend'],ncol=1,title=title)
	for legobj in legend.legendHandles: legobj.set_linewidth(lw)

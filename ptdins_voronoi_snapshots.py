#!/usr/bin/env python

#!!! clear plotter-membrane.py when you are done and also script-aamd-apl.py

import scipy
import scipy.spatial
import random
from codes.voronoi_mesh_snapshot import plot_snapshots

@autoload(plotrun)
def load():
	data = plotload(plotname)
	sns = work.sns()

@autoplot(plotrun)
def plot(): 
	plot_snapshots(sn=work.sns()[0],data=data)

if False:
	#! ridiculous hack via: https://stackoverflow.com/questions/11578760
	import types
	from matplotlib.backend_bases import GraphicsContextBase, RendererBase
	import matplotlib.pyplot as plt
	from matplotlib.collections import LineCollection
	class GC(GraphicsContextBase):
	    def __init__(self):
	        super(GC,self).__init__()
	        self._capstyle = 'round'
	def custom_new_gc(self): return GC()
	RendererBase.new_gc = types.MethodType(custom_new_gc,RendererBase)

if __name__=='__main__': pass
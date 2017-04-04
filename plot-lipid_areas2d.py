#!/usr/bin/env python

routine = ['summary']

#---load everything
if 'data' not in globals(): 
	data,calcs = plotload(plotname,work)
	ns = next_script = sys.argv[0]

print(data['mk001']['data']['areas0'].mean(axis=0))
print('VORONOI AREA EXAMPLE ^^^^')
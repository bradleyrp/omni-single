#!/usr/bin/env python

import codes.curvature_coupling.InvestigateCurvature
codes.curvature_coupling.InvestigateCurvature.plotload = plotload
codes.curvature_coupling.InvestigateCurvature.plotname = plotname
codes.curvature_coupling.InvestigateCurvature.work = work
#---alternate loader
from import_regular_membrane import plotloader_for_dextran
ic = codes.curvature_coupling.InvestigateCurvature.InvestigateCurvature(plotloader=plotloader_for_dextran)

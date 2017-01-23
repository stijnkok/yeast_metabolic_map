#!/usr/bin/python

import sys
import os
import cbmpy as cbm
import vonda


# get model (Yeast 7 consensus model)
print 'loading Yeast 7 model...'
Y7 = cbm.CBRead.readSBML3FBC('models/Y7.xml')
Y7.name = 'Y7'
cbm.doFBA(Y7)
cbm.doFBAMinSum(Y7)
# make vonda mod
vmod = vonda.PVisualizer(
    'Yeast_7.svg',
    SVG_dir=os.path.join(os.getcwd(), 'interactive_maps'), 
    reaction_suffix='r_',
    species_suffix='s_')
Y7_fluxes = Y7.getReactionValues()
vmod.doMapReactions(
    Y7_fluxes, 
    filename_out = "Yeast7_aerobic0", 
    IsAbsoluteValues=False, 
    valuesRange=['outside',-1e-09, 1e-09], 
    minValue=0.0001, 
    maxValue=10)

# To get a flux through the phosphofructokinase and TCA cycle,
# make cytoplasmic isocitrate dehydrogenase irreversible:
Y7.setReactionLowerBound('r_0659', 0.0)
cbm.doFBA(Y7)
cbm.doFBAMinSum(Y7)
Y7_fluxes = Y7.getReactionValues()
vmod.doMapReactions(
    Y7_fluxes, 
    filename_out = "Yeast7_aerobic1", 
    IsAbsoluteValues=False, 
    valuesRange=['outside',-1e-09, 1e-09], 
    minValue=0.0001, 
    maxValue=10)
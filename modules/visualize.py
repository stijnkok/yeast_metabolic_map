import os
import cbmpy as cbm
import re
from numpy import log
import webbrowser

class Vmod:

	def __init__(self, svg, model, r_prefix='r_'):
		
		self.svg = svg
		self.model = model
		self.r_prefix = r_prefix

		with open(svg, 'rb') as f:
			svg_lines=f.readlines()
		self.svg_lines = svg_lines



	def mapFBA(self, D_fluxes=None, D_bounds=None, absminval=1e-15, out_file='FBA_result.svg'):

		svg_lines = [line for line in self.svg_lines]
		r_prefix = self.r_prefix
		
		if self.model and not D_fluxes:
			D_fluxes=self.model.getReactionValues()
		
		if self.model and not D_bounds:
			D_bounds = {}
			for r in self.model.reactions:
				D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())

		rids = set(re.findall('id="({}\w+)"'.format(r_prefix), ''.join(svg_lines)))
		maxval = max([abs(D_fluxes[rid]) for rid in rids])
		minval = min([abs(D_fluxes[rid]) for rid in rids if abs(D_fluxes[rid]) != 0])
		minval = max(minval, absminval)


		for i in range(len(svg_lines)):
			line = svg_lines[i]
			m = re.search('#({}\w+)(\S*)\s*({{.*}})'.format(r_prefix), line)
			if m:
				rid = m.group(1)
				rclass = m.group(2)
				rstyle = m.group(3)
				
				val = D_fluxes[rid]
				if D_bounds:
					if D_bounds[rid][0] == 0 and D_bounds[rid][1] == 0:
						if '.reversible' in rclass or '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.inactive' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')
					elif D_bounds[rid][0] == 0 or D_bounds[rid][1] == 0:
						if '.reversible' in rclass or '.inactive' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')
					else:
						if '.inactive' in rclass or '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.reversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')

				if '.substrate' in rclass and val <0:
					svg_lines[i] = line.replace('#substrate', '#product')
				elif '.product' in rclass and val <0:
					svg_lines[i] = line.replace('#product', '#substrate')
				elif '.fluxvalue_tooltip' in rclass:
					svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif '.fluxvalue' in rclass:
					if abs(val) < minval:
						svg_lines[i] = line.replace('fill-opacity:1', 'fill-opacity:0')
					else:
						svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif not rclass:
					if abs(val) < minval:
						svg_lines[i] = line.replace(rstyle, '{stroke:#cccccc; stroke-width:1.0; stroke-dasharray:1.5}')
					else:
						hue = 240 - (log(abs(val))-log(minval))/(log(maxval)-log(minval))*240
						svg_lines[i] = line.replace(rstyle, '{{stroke:hsl({}, 100%, 50%); stroke-width:2.0}}'.format(hue))

			m=re.search('ReactionValue:(\w*)(\d+):({}\w+)'.format(r_prefix), line)
			if m:
				abs_str = m.group(1)
				d = int(m.group(2))
				rid = m.group(3)	
				val = D_fluxes[rid]
				if abs_str=='abs':
					val=abs(val)
				if abs(val) < minval:
					svg_lines[i] = line.replace(m.group(), '0')
				else:
					val_str = '{{:0.{}e}}'.format(d-1).format(val)
					svg_lines[i] = line.replace(m.group(), val_str)
		

		with open(out_file, 'wb') as f:
			f.write(''.join(svg_lines))

		webbrowser.open_new_tab(os.path.join(os.getcwd(), out_file))

	def mapFVA(self, fva_result, D_bounds=None, absminval=1e-15, minspan=1e-15, out_file='FVA_result.svg'):
		
		svg_lines = [line for line in self.svg_lines]
		r_prefix = self.r_prefix

		D_fva = {}
		for i in range(len(fva_result[1])):
			D_fva[fva_result[1][i]] = fva_result[0][i]

		if self.model and not D_bounds:
			D_bounds = {}
			for r in self.model.reactions:
				D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())

		rids = set(re.findall('id="({}\w+)"'.format(r_prefix), ''.join(svg_lines)))
		maxspan = max([D_fva[rid][4] for rid in rids])
		minspan_ = min([D_fva[rid][4] for rid in rids if D_fva[rid][4] != 0])
		minspan = max(minspan, minspan_)
		minval = min([abs(D_fva[rid][0]) for rid in rids if abs(D_fva[rid][0]) != 0])
		minval = max(minval, absminval)

		for i in range(len(svg_lines)):
			line = svg_lines[i]
			m = re.search('#({}\w+)(\S*)\s*({{.*}})'.format(r_prefix), line)
			if m:
				rid = m.group(1)
				rclass = m.group(2)
				rstyle = m.group(3)
				
				val = D_fva[rid][0]
				span = D_fva[rid][4]

				if D_bounds:
					if D_bounds[rid][0] == 0 and D_bounds[rid][1] == 0:
						if '.reversible' in rclass or '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.inactive' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')
					elif D_bounds[rid][0] == 0 or D_bounds[rid][1] == 0:
						if '.reversible' in rclass or '.inactive' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')
					else:
						if '.inactive' in rclass or '.irreversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:1', 'stroke-opacity:0')
						if '.reversible' in rclass:
							svg_lines[i] = line.replace('stroke-opacity:0', 'stroke-opacity:1')
				if '.substrate' in line and val <0:
					svg_lines[i] = line.replace('#substrate', '#product')
				elif '.product' in line and val <0:
					svg_lines[i] = line.replace('#product', '#substrate')
				elif '.fluxvalue_tooltip' in line:
					svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif '.fluxvalue' in line:
					svg_lines[i] = line.replace('fill-opacity:1', 'fill-opacity:0')
				elif '.FVAspan_tooltip' in line:
					svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif '.FVAspan' in line:
					if span < minspan:
						svg_lines[i] = line.replace('fill-opacity:1', 'fill-opacity:0')
					else:
						svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif '.FVAmin' in line:
					if D_bounds:
						if val == D_bounds[rid][0] == D_fva[rid][2]:
							svg_lines[i] = line.replace(rstyle, '{fill:#ff0000; fill-opacity:1}')
					svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif '.FVAmax' in line:
					if D_bounds:
						if val == D_bounds[rid][1] == D_fva[rid][3]:
							svg_lines[i] = line.replace(rstyle, '{fill:#ff0000; fill-opacity:1}')
					svg_lines[i] = line.replace('fill-opacity:0', 'fill-opacity:1')
				elif not rclass:
					if span < minspan and abs(val) < minval:
						svg_lines[i] = line.replace(rstyle, '{stroke:#cccccc; stroke-width:1.0; stroke-dasharray:1.5}')
					elif span < minspan:
						svg_lines[i] = line.replace(rstyle, '{stroke:#cccccc; stroke-width:2.0}')
					else:
						hue = 240 - (log(span)-log(minspan))/(log(maxspan)-log(minspan))*240
						svg_lines[i] = line.replace(rstyle, '{{stroke:hsl({}, 100%, 50%); stroke-width:2.0}}'.format(hue))

			m=re.search('ReactionValue:(\w*)(\d+):({}\w+)'.format(r_prefix), line)
			if m:
				abs_str = m.group(1)
				d = int(m.group(2))
				rid = m.group(3)	
				val = D_fva[rid][0]
				if abs_str=='abs':
					val=abs(val)
				if abs(val) < minval:
					svg_lines[i] = line.replace(m.group(), '0')
				else:
					val_str = '{{:0.{}e}}'.format(d-1).format(val)
					svg_lines[i] = line.replace(m.group(), val_str)

			m=re.search('ReactionSpan:(\d+):({}\w+)'.format(r_prefix), line)
			if m:
				d = int(m.group(1))
				rid = m.group(2)	
				span = D_fva[rid][4]
				if span < minspan:
					svg_lines[i] = line.replace(m.group(), '0')
				else:
					val_str = '{{:0.{}e}}'.format(d-1).format(span)
					svg_lines[i] = line.replace(m.group(), val_str)

			m=re.search('ReactionMinValue:(\d+):({}\w+)'.format(r_prefix), line)
			if m:
				d = int(m.group(1))
				rid = m.group(2)	
				val = D_fva[rid][2]
				val_str = '{{:0.{}e}}'.format(d-1).format(val)
				svg_lines[i] = line.replace(m.group(), val_str)

			m=re.search('ReactionMaxValue:(\d+):({}\w+)'.format(r_prefix), line)
			if m:
				d = int(m.group(1))
				rid = m.group(2)
				val = D_fva[rid][3]
				val_str = '{{:0.{}e}}'.format(d-1).format(val)
				svg_lines[i] = line.replace(m.group(), val_str)


		with open(out_file, 'wb') as f:
			f.write(''.join(svg_lines))

		webbrowser.open_new_tab(os.path.join(os.getcwd(), out_file))

def main():

	model = cbm.CBRead.readSBML3FBC('models/Y7.xml')
	vmod = Vmod('metabolic_maps/Yeast_7.svg', model)

	cbm.doFBA(model)
	cbm.doFBAMinSum(model)
	vmod.mapFBA(absminval=1e-5, out_file = 'metabolic_maps/Y7_fba_default.svg')

	model.setReactionLowerBound('r_0659', 0.0)
	cbm.doFBA(model)
	cbm.doFBAMinSum(model)
	vmod.mapFBA(absminval=1e-5, out_file = 'metabolic_maps/Y7_fba_irreversible_IDH.svg')

	FVA = cbm.doFVA(model)
	vmod.mapFVA(FVA, absminval=1e-5, minspan=1e-5, out_file = 'metabolic_maps/Y7_fva.svg')

	gene_ids = ['YER081W', 'YIL074C', 'YKL141W']
	for gene_id in gene_ids:
		model.setGeneInactive(gene_id)
	model.updateNetwork()
	cbm.doFBA(model)
	cbm.doFBAMinSum(model)	
	vmod.mapFBA(absminval=1e-5, out_file = 'metabolic_maps/Y7_fba_Otero_2013_succinate.svg')


	# delft stratagy

	model = cbm.CBRead.readSBML3FBC('models/Y7.xml')
	vmod = Vmod('metabolic_maps/Yeast_7.svg', model)
	model.createObjectiveFunction('r_2056')
	cbm.doFBA(model)
	cbm.doFBAMinSum(model)
	vmod.mapFBA(absminval=1e-5, out_file = 'metabolic_maps/Y7_fba_succinate_obj_default.svg')

	vmod = Vmod('metabolic_maps/Yeast_7delft.svg', model)
	model.setReactionLowerBound('r_0659', 0)
	model.setReactionBounds('r_0959', 0, 0) # knock out pyruvate decarboxylase
	model.createReaction('r_FRDc', name='cytosolic fumarate reductase')
	model.createReactionReagent('r_FRDc', 's_0725', -1) # fumarate
	model.createReactionReagent('r_FRDc', 's_1203', -1)  # NADH
	model.createReactionReagent('r_FRDc', 's_1458', 1)  # succinate
	model.createReactionReagent('r_FRDc', 's_1198', 1) # NAD

	model.createReaction('r_PDHc', name='bacterial PDH')
	model.createReactionReagent('r_PDHc', 's_1399', -1) # pyruvate
	model.createReactionReagent('r_PDHc', 's_1198', -1) # NAD
	model.createReactionReagent('r_PDHc', 's_0529', -1) # CoA
	model.createReactionReagent('r_PDHc', 's_0373', 1)  # Acetyl-CoA
	model.createReactionReagent('r_PDHc', 's_1203', 1) # NADH
	model.createReactionReagent('r_PDHc', 's_0456', 1) # CO2
	model.setReactionLowerBound('r_PDHc', 0)

	vmod.mapFBA(absminval=1e-5, out_file = 'metabolic_maps/Y7_fba_succinate_obj_delft.svg')

	FVA = cbm.doFVA(model)
	vmod.mapFVA(FVA, absminval=1e-5, minspan=1e-5, out_file = 'metabolic_maps/Y7_fva_succinate_obj_delft.svg.svg')

if __name__ == '__main__':
	main()	


	# # D_fluxes=model.getReactionValues()

	# # with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# # 	svg_lines=f.readlines()

	# # svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', absminval=1e-5)

	# # with open('metabolic_maps/Y7fba_aerobic.svg', 'wb') as f:
	# # 	f.write(''.join(svg_lines))

	# # webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_aerobic.svg'))

	# # with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# # 	svg_lines=f.readlines()

	# # FVA = cbm.doFVA(model)
	# # D_bounds = {}
	# # for r in model.reactions:
	# # 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())
	# # svg_lines = mapFVA(svg_lines, FVA, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5, minspan=1e-5)

	# # with open('metabolic_maps/Y7fva_aerobic.svg', 'wb') as f:
	# # 	f.write(''.join(svg_lines))

	# # webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fva_aerobic.svg'))

	# # succinate production

	# gene_ids = ['YER081W', 'YIL074C', 'YKL141W']
	# for gene_id in gene_ids:
	# 	model.setGeneInactive(gene_id)
	# model.updateNetwork()
	# model.setReactionLowerBound('r_0659', 0.0)
	# cbm.doFBA(model)
	# cbm.doFBAMinSum(model)
	# D_fluxes=model.getReactionValues()
	# D_bounds = {}
	# for r in model.reactions:
	# 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())

	# with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# 	svg_lines=f.readlines()	
	
	# svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5)

	# with open('metabolic_maps/Y7fba_succinate_stratagy_1.svg', 'wb') as f:
	# 	f.write(''.join(svg_lines))	

	# webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_succinate_stratagy_1.svg'))

	# model = cbm.CBRead.readSBML3FBC('models/Y7.xml')
	# model.setReactionLowerBound('r_0659', 0.0)
	# # model.setReactionBounds('r_0959', 0, 0) # knock out pyruvate decarboxylase
	# model.createReaction('r_FRDc', name='cytosolic fumarate reductase')
	# model.createReactionReagent('r_FRDc', 's_0725', -1) # fumarate
	# model.createReactionReagent('r_FRDc', 's_1203', -1)  # NADH
	# model.createReactionReagent('r_FRDc', 's_1458', 1)  # succinate
	# model.createReactionReagent('r_FRDc', 's_1198', 1) # NAD
	# model.setReactionBounds('r_FRDc', 0, 0)

	# model.createReaction('r_PDHc', name='bacterial PDH')
	# model.createReactionReagent('r_PDHc', 's_1399', -1) # pyruvate
	# model.createReactionReagent('r_PDHc', 's_1198', -1) # NAD
	# model.createReactionReagent('r_PDHc', 's_0529', -1) # CoA
	# model.createReactionReagent('r_PDHc', 's_0373', 1)  # Acetyl-CoA
	# model.createReactionReagent('r_PDHc', 's_1203', 1) # NADH
	# model.createReactionReagent('r_PDHc', 's_0456', 1) # CO2
	# model.setReactionBounds('r_PDHc', 0, 0)

	# cbm.doFBA(model)
	# cbm.doFBAMinSum(model)
	# D_fluxes=model.getReactionValues()
	# D_bounds = {}
	# for r in model.reactions:
	# 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())
	
	# with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# 	svg_lines=f.readlines()	
	
	# svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5)

	# with open('metabolic_maps/Y7fba_succinate_stratagy_2.svg', 'wb') as f:
	# 	f.write(''.join(svg_lines))	

	# webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_succinate_stratagy_2.svg'))

	# gene_ids = ['YER081W', 'YIL074C', 'YKL141W']
	# for gene_id in gene_ids:
	# 	model.setGeneInactive(gene_id)
	# model.updateNetwork()
	# model.setReactionLowerBound('r_0659', 0.0)
	# cbm.doFBA(model)
	# cbm.doFBAMinSum(model)
	# D_fluxes=model.getReactionValues()
	# D_bounds = {}
	# for r in model.reactions:
	# 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())

	# with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# 	svg_lines=f.readlines()	
	
	# svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5)

	# with open('metabolic_maps/Y7fba_succinate_stratagy_1and2.svg', 'wb') as f:
	# 	f.write(''.join(svg_lines))	

	# webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_succinate_stratagy_1and2.svg'))

	# # succinate maximum yield

	# # model = cbm.CBRead.readSBML3FBC('models/Y7.xml')
	# # model.createObjectiveFunction('r_2056')
	# # cbm.doFBA(model)
	# cbm.doFBAMinSum(model)
	# D_fluxes=model.getReactionValues()
	# D_bounds = {}
	# for r in model.reactions:
	# 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())
	
	# with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# 	svg_lines=f.readlines()	
	
	# svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5)

	# with open('metabolic_maps/Y7fba_succinate_max.svg', 'wb') as f:
	# 	f.write(''.join(svg_lines))	

	# webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_succinate_max.svg'))	

	# model.setReactionLowerBound('r_0659', 0.0)
	# model.setReactionUpperBound('r_0959', 0.0)

	# model.createReaction('R_NEW_FRD_c', name='cytosolic fumarate reductase')
	# model.createReactionReagent('R_NEW_FRD_c', 's_0725', -1) # fumarate
	# model.createReactionReagent('R_NEW_FRD_c', 's_1203', -1) # NADH
	# model.createReactionReagent('R_NEW_FRD_c', 's_1458', 1)  # succinate
	# model.createReactionReagent('R_NEW_FRD_c', 's_1198', 1)  # NAD

	# model.createReaction('R_NEW_PDH_c', name='bacterial PDH')
	# model.createReactionReagent('R_NEW_PDH_c', 's_1399', -1) # pyruvate
	# model.createReactionReagent('R_NEW_PDH_c', 's_1198', -1) # NAD
	# model.createReactionReagent('R_NEW_PDH_c', 's_0529', -1) # CoA
	# model.createReactionReagent('R_NEW_PDH_c', 's_0373', 1)  # Acetyl-CoA
	# model.createReactionReagent('R_NEW_PDH_c', 's_1203', 1) # NADH
	# model.createReactionReagent('R_NEW_PDH_c', 's_0456', 1) # CO2
	# model.setReactionLowerBound('R_NEW_PDH_c', 0.0)

	# cbm.doFBA(model)
	# cbm.doFBAMinSum(model)
	# D_fluxes=model.getReactionValues()
	# D_bounds = {}
	# for r in model.reactions:
	# 	D_bounds[r.id] = (r.getLowerBound(), r.getUpperBound())
	
	# with open('metabolic_maps/Yeast_7.svg', 'rb') as f:
	# 	svg_lines=f.readlines()	
	
	# svg_lines = mapFBA(svg_lines, D_fluxes, r_prefix='r_', D_bounds=D_bounds, absminval=1e-5)

	# with open('metabolic_maps/Y7fba_succinate_max_delft.svg', 'wb') as f:
	# 	f.write(''.join(svg_lines))	

	# webbrowser.open_new_tab(os.path.join(os.getcwd(), 'metabolic_maps','Y7fba_succinate_max_delft.svg'))	


	
    #r_2148 {stroke:#cccccc; stroke-width:1.0; stroke-dasharray:1.5}
    #r_2148.substrate {marker-end:url(#substrate)}
    #r_2148.product   {marker-end:url(#product)}
    #r_2148.fluxvalue {fill-opacity:0}
    #r_2148.FVAspan {fill-opacity:0}
    #r_2148.FVAmin {fill-opacity:0}
    #r_2148.FVAmax {fill-opacity:0}
    #r_2148.reversible {stroke-opacity:1}
    #r_2148.irreversible {stroke-opacity:0}
    #r_2148.inactive {stroke-opacity:0}

# def mapFBA(svg_doc, D_fluxes, r_prefix=None):
# 	if r_prefix:
# 		rids = set(re.findall('id="({}\w+)"'.format(r_prefix), svg_doc))
# 	else:
# 		rids = D_fluxes.keys()
# 	maxval = max([abs(D_fluxes[rid]) for rid in rids])
# 	# minval = min([abs(D_fluxes[rid]) for rid in rids if abs[D_fluxes[rid]] != 0])
# 	minval = 1e-5

# 	pstroke = re.compile('stroke:.*?;')
# 	pstroke_width = re.compile('stroke-width:.*?;')
# 	pstroke_dasharray = re.compile('stroke-dasharray:.*?;')

# 	for rid in rids:
# 		val = D_fluxes[rid]
# 		m = re.search('#{}\s*{{.*?}}'.format(rid), svg_doc, re.DOTALL)
# 		if m:
# 			if abs(val) < minval:
# 				new_stroke = 'stroke:#cccccc;'
# 				new_stroke_width = 'stroke-width:1.0;'
# 				new_stroke_dasharray = 'stroke-dasharray:1.5'
# 			else:
# 				hue = 240 - (log(abs(val))-log(minval))/(log(maxval)-log(minval))*240
# 				new_stroke = 'stroke:hsl({}, 100%, 50%);'.format(hue)
# 				new_stroke_width = 'stroke-width:2.0;'
# 				new_stroke_dasharray = ''

# 			newstyle = m.group()
# 			stroke = pstroke.search(newstyle)
# 			if stroke:
# 				newstyle = newstyle.replace(stroke.group(), new_stroke)
# 			else:
# 				newstyle = newstyle.replace('{', '{'+new_stroke)

# 			stroke_width = pstroke_width.search(newstyle)
# 			if stroke_width:
# 				newstyle = newstyle.replace(stroke_width.group(), new_stroke_width)
# 			else:
# 				newstyle = newstyle.replace('{', '{'+new_stroke_width)
			
# 			stroke_dasharray = pstroke_dasharray.search(newstyle)
# 			if stroke_dasharray:
# 				newstyle = newstyle.replace(stroke_dasharray.group(), new_stroke_dasharray)
# 			else:
# 				newstyle = newstyle.replace('{', '{'+new_stroke_dasharray)

# 			svg_doc = svg_doc.replace(m.group(), newstyle)

# 		m = re.search('#{}[.]fluxvalue_tooltip\s*{{.*?}}'.format(rid), svg_doc, re.DOTALL)
# 		if m:
# 			svg_doc = svg_doc.replace(m.group(), m.group().replace('fill-opacity:0', 'fill-opacity:1'))

# 		m = re.search('#{}[.]fluxvalue\s*{{.*?}}'.format(rid), svg_doc, re.DOTALL)
# 		if m:
# 			if abs(val)<minval:
# 				svg_doc = svg_doc.replace(m.group(), m.group().replace('fill-opacity:1', 'fill-opacity:0'))
# 			else:
# 				svg_doc = svg_doc.replace(m.group(), m.group().replace('fill-opacity:0', 'fill-opacity:1'))

# 		if val <0:
# 			m = re.search('#{}[.]substrate\s{{.*?}}'.format(rid), svg_doc, re.DOTALL)
# 			if m:
# 				svg_doc = svg_doc.replace(m.group(), m.group().replace('#substrate','#product'))
# 			m = re.search('#{}[.]product\s*{{.*?}}'.format(rid), svg_doc, re.DOTALL)
# 			if m:
# 				svg_doc = svg_doc.replace(m.group(), m.group().replace('#product','#substrate'))

# 		for m in set(re.findall('ReactionValue:((\w*)(\d+)):{}'.format(rid), svg_doc)):
# 			d = int(m[2])
# 			if m[1]=='abs':
# 				val_str = '{{:0.{}e}}'.format(d-1).format(abs(val))
# 			else:
# 				val_str = '{{:0.{}e}}'.format(d-1).format(val)
# 			if abs(val) < minval:
# 				svg_doc = svg_doc.replace('ReactionValue:'+m[0]+':'+rid, '0')
# 			else:
# 				svg_doc = svg_doc.replace('ReactionValue:'+m[0]+':'+rid, val_str)
	
# 	return svg_doc
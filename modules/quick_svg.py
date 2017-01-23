import sys
import os
import networkx as nx

from pysvg.core import TextContent
from pysvg.structure import svg, g
from pysvg.text import text
from pysvg.builders import StyleBuilder
from pysvg.shape import circle, path
from pysvg.parser import parse


graph_thf = nx.read_graphml('Y7thf.graphml')

rids = [nid for nid in graph_thf.nodes() if nid.startswith('r')]

scale = 1.2
for nid in graph_thf.nodes():
	n = graph_thf.node[nid]
	n['x']=n['x']*scale
	n['y']=n['y']*scale

circles = []
paths   = []
texts   = []
sids = []

for rid in rids:
	n = graph_thf.node[rid]
	svgcircle = circle(
		cx = n['x'],
		cy =-n['y'],
		r  = 5)
	svgcircle.set_id(rid)
	svgcircle.set_style('fill:#999999')
	circles.append(svgcircle)
	svgtext = text(
		content = '00e+00',
		x = n['x']+10,
		y = -n['y'])
	svgtext.set_id('rval_'+rid)
	svgtext.set_style('fill:black;font-size:8px;font-family:Raleway')
	texts.append(svgtext)

	for mid in graph_thf.predecessors(rid):
		m = graph_thf.node[mid]
		if mid not in sids:
			# substrates: labels
			svgtext = text(
				content = m['sLabel'],
				x = m['x'],
				y = -m['y'])
			svgtext.set_id(mid)
			svgtext.set_style('fill:black;font-size:10px;font-family:Raleway')
			texts.append(svgtext)
			sids.append(mid)

		# substrates: arrows
		d_x = m['x']-n['x']
		d_y = m['y']-n['y']
		svgpath = path('m {},{} {},{}'.format(n['x'], -n['y'], d_x, -d_y))
		svgpath.set_id('path_'+rid+m['label'][:6])
		svgpath.set_style('fill:none;stroke:#cccccc;stroke-width:2;stroke-linecap:round;marker-end:url(#substrate)')
		paths.append(svgpath)

	for mid in graph_thf.successors(rid):
		m = graph_thf.node[mid]
		if mid not in sids:
			# substrates: labels
			svgtext = text(
				content = m['sLabel'],
				x = m['x'],
				y = -m['y'])
			svgtext.set_id(mid)
			svgtext.set_style('fill:black;font-size:10px;font-family:Raleway')
			texts.append(svgtext)
			sids.append(mid)

		# substrates: arrows
		d_x = m['x']-n['x']
		d_y = m['y']-n['y']
		svgpath = path('m {},{} {},{}'.format(n['x'], -n['y'], d_x, -d_y))
		svgpath.set_id('path_'+rid+m['label'][:6])
		svgpath.set_style('fill:none;stroke:#cccccc;stroke-width:2;stroke-linecap:round;marker-end:url(#product)')
		paths.append(svgpath)

defs = TextContent("""
<defs>
  <marker id="product" markerWidth="6" markerHeight="6" refX="0.75" refY="2.5" orient="auto" markerUnits="strokeWidth">
	<path d="M1,1 Q3,2.5 1,4 L4,2.5 z " stroke-linejoin="round" stroke-linecap="round" stroke="#000" stroke-width="0.5" fill="#000"/>
  </marker>
  <marker id="substrate" markerWidth="6" markerHeight="6" refX="0.75" refY="2.5" orient="auto" markerUnits="strokeWidth">
	<circle cx="2.6" cy="2.5" r="1.1" fill="#000"/>
  </marker>
</defs>
""")

svgdoc= svg()
for e in [defs] + paths + circles + texts:
	svgdoc.addElement(e)
svgdoc.save('Y7_thf.svg')



#######
""""

need to add 
r_0506
r_0507
r_0508


"""

cofactors = ['s_0799', 's_1005', 's_0460', 's_1488', 's_0307', 's_0421', 's_1200', 's_1205']

import cbmpy as cbm
Y7 = cbm.CBRead.readSBML3FBC('models/Y7.xml')

substrate = """<text
  style="font-size:10px;font-family:Raleway"
  id="{s_id}"
  y="0"
  x="0">{s_name}</text>
<path
  d="m 0,25 0,-25"
  id="path_{r_id}{s_id}"
  inkscape:connector-curvature="0"
  style="fill:none;stroke:#cccccc;stroke-width:2;stroke-linecap:round;marker-end:url(#substrate)"
  sodipodi:nodetypes="cc" />
"""

r_circle = """  <circle
     cy="25"
     cx="0"
     r="5"
     id="{r_id}"
     style="fill:#999999" />
  <text
     style="font-size:8px;font-family:Raleway;fill:#000000"
     id="rval_{r_id}"
     y="25"
     x="0">-00e+00</text>
"""
product = """<text
  style="font-size:10px;font-family:Raleway"
  id="{s_id}"
  y="50"
  x="0">{s_name}</text>
<path
  d="m 0,25 0,25"
  id="path_{r_id}{s_id}"
  inkscape:connector-curvature="0"
  style="fill:none;stroke:#cccccc;stroke-width:2;stroke-linecap:round;marker-end:url(#product)"
  sodipodi:nodetypes="cc" />
"""
# cof paths
d_t = ['m 0,25 c 0,-25 10,-25 20,-25', 'm 0,25 c 0,-25 -10,-25 -20,-25']
d_b = ['m 0,25 c 0,25 10,25 20,25', 'm 0,25 c 0,25 -10,25 -20,25']
#cof coords
xy_t = [(-25, 0), (25, 0)]
xy_b = [(-25, 50), (25, 50)]

cof = """<text
  style="font-size:10px;font-family:Raleway"
  id="{s_id}{r_id}"
  y="{y}"
  x="{x}">{s_name}</text>
<path
  d="{d}"
  id="path_{r_id}{s_id}{r_id}"
  inkscape:connector-curvature="0"
  style="fill:none;stroke:#cccccc;stroke-width:2;stroke-linecap:round;marker-end:url(#{role})"
  sodipodi:nodetypes="cc" />
"""

svg_rxns=[]
for rid in ['r_0506','r_0507','r_0508']:
	r = Y7.getReaction(rid)
	sid = [s for s in r.getSubstrateIds() if s not in cofactors][0]
	s_name = Y7.getSpecies(sid).name
	cf_sids = [s for s in r.getSubstrateIds() if s in cofactors]
	cf_s_names = [Y7.getSpecies(cf_id).name for cf_id in cf_sids]

	pid = [p for p in r.getProductIds() if p not in cofactors][0]
	p_name = Y7.getSpecies(pid).name
	cf_pids = [p for p in r.getProductIds() if p in cofactors]
	cf_p_names = [Y7.getSpecies(cf_id).name for cf_id in cf_pids]

	svg_rxn = substrate.format(s_id=sid, s_name=s_name, r_id=rid)
	svg_rxn+= product.format(s_id=pid, s_name=p_name, r_id=rid)
	for i in range(len(cf_sids)):
		svg_rxn+= cof.format(s_id=cf_sids[i], r_id=rid, x= xy_t[i][0], y= xy_t[i][1], s_name=cf_s_names[i], d=d_t[i], role='substrate')
	for i in range(len(cf_pids)):
		svg_rxn+= cof.format(s_id=cf_pids[i], r_id=rid, x= xy_b[i][0], y= xy_b[i][1], s_name=cf_p_names[i], d=d_b[i], role='product')
	svg_rxn+= r_circle.format(r_id=rid)
	svg_rxn = svg_rxn.replace(' [mitochondrion]', '')
	svg_rxns.append(svg_rxn)

# r_0990 'sedoheptulose 1,7-bisphosphate [cytoplasm] = dihydroxyacetone phosphate [cytoplasm] + D-erythrose 4-phosphate [cytoplasm] '
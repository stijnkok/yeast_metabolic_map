import json
import argparse
import itertools
from PIL import ImageFont
from svg_paths import get_paths
from pysvg.structure import svg, g
from pysvg.text import text
from pysvg.shape import path, circle
from pysvg.core import TextContent


def get_label_width(label, font, font_size):
	font_size = float(font_size)
	width, _ = font.getsize(label)
	width = width*(font_size/font.size)
	return width

def get_svgdata(d, font, font_size, scale, padding, padding_labels, normalize, overlap, cofactors = None, cap_labels = True, scale_labels= False, defdir='v', reverse_cof = []):
	""" Get all information to make an svg-file with the metabolic map. Output in a dictionary.
	d 				Dictionary with information extracted from json-output from Nicholas (dictionary)
						>> See read_json_data function in json_to_svg.py
	font 			Font for the labels (ImageFont truetype object)
	font_size 		Font size (int/float)
	scale 			Scale factor(s) for all coordinates. If only one number is given scaling is done in 
					both x and y direction with the given scale factor. If a list of numbers is given, the
					first number is scale factor in x direction and second number in y direction (list of int/float)
	padding 		Extra white space added to the borders of the whole svg document (list of int/float)
	padding_labels	White space around label text (list of int/float)
	normalize		Translates all coordinates to positive coordinates, so everything will be visible when
					the svg is opened in a browser (bool)
	overlap 		Allow for overlapping labels and lines (bool)
	cofactors		Dictionary with cofactors per reaction (dictionary)
						>> See get_cofactors_from_sbml function in parse_cofactors.py
	cap_labels		Cap labels that are very long and might cause overlap (bool)
	scale_labels	Scale labels that are very long (bool)
	defdir			Default direction of arrows from reaction nodes ('v'/'vertical' or 'h'/'horizontal')
	reverse_cof		List of reaction nodes that should have the cofactors placed in reverse direction compared to default

	OUTPUT dictionary keys: 'rxn_nodes' (reaction nodes), 'paths' (svg paths), 'labels', 'font_size', 'font_family'
	"""

	if len(scale)==1:
		scale.append(scale[0])
	if len(padding) == 1:
		padding.append(padding[0])
	if padding_labels:
		if len(padding_labels) == 1:
			padding_labels.append(padding_labels[0])
	else:
		padding_labels = [font_size, font_size]

	role = d.pop('edge_type')

	# adjust coordinates
	min_x = min([p[0] for p in d['pos'].values()])
	max_y = max([p[1] for p in d['pos'].values()])

	for n, p in d['pos'].items():
		if normalize:
			p = (p[0]-min_x, p[1]-max_y) # normalize
		p = (p[0]*scale[0], p[1]*scale[1]) # scale
		p = (p[0] + padding[0], p[1] - padding[1]) # add padding
		p = (p[0], -p[1]) # convert to svg coordinates
		d['pos'][n] = p

	for e in d['extra_nodes']:
		for i in range(len(d['extra_nodes'][e])):
			p = d['extra_nodes'][e][i]
			if normalize:
				p = (p[0]-min_x, p[1]-max_y) # normalize
			p = (p[0]*scale[0], p[1]*scale[1]) # scale
			p = (p[0] + padding[0], p[1] - padding[1]) # add padding
			p = (p[0], -p[1]) # convert to svg coordinates
			d['extra_nodes'][e][i] = p

	# get max width of labels
	snpos = [d['pos'][n] for n in d['nodes'] if d['node_type'][n]=='species']
	snpos.sort(key = lambda p: p[1])
	dx = []
	for _, group in itertools.groupby(snpos, key = lambda p: p[1]):
		dx+= [abs(c[0]-c[1]) for c in itertools.combinations([p[0] for p in group], 2)]
	max_w = min(dx) + 2*padding_labels[0]
	
	# get labels, label size
	if cofactors:
		for r in cofactors:
			for s in cofactors[r]:
				d['label'][s] = cofactors[r][s]['label']
	d['label_size'] = {}
	link_text = {}
	fs = {}
	if cap_labels and get_label_width('...', font, font_size) + 2*padding_labels[0] > max_w:
		print 'Please scale up or use smaller font size to prevent overlapping labels'
		for n in d['label']:
			link_text[n] = d['label'][n]
			d['label'][n] = '...'
			w = get_label_width(d['label'][n], font, font_size) + 2*padding_labels[0]
			h = font_size + 2*padding_labels[1]
			d['label_size'][n] = (w, h)
			fs[n] = font_size

	elif cap_labels:
		for n in d['label']:
			link_text[n]= d['label'][n]
			lab = d['label'][n].replace(' [cytoplasm]', '') # yeast concensus models adjustment
			w = get_label_width(lab, font, font_size) + 2*padding_labels[0]
			# cap too long labels
			if w > max_w:
				for i in range(len(lab),0,-1):
					newlab = lab[:i]+'...'
					neww = get_label_width(newlab, font, font_size) + 2*padding_labels[0]
					if neww <= max_w:
						lab = newlab
						break
				else:
					lab = '...'
			d['label'][n] = lab
			w = get_label_width(d['label'][n], font, font_size) + 2*padding_labels[0]
			h = font_size + 2*padding_labels[1]
			d['label_size'][n] = (w, h)
			fs[n] = font_size
	
	elif scale_labels:
		for n in d['label']:
			link_text[n]= d['label'][n]
			lab = d['label'][n].replace(' [cytoplasm]', '')
			d['label'][n] = lab
			w = get_label_width(lab, font, font_size)
			# scale too long labels
			if w > max_w:
				fs[n] = font_size * max_w/w
				w = max_w
			h = font_size + 2*padding_labels[1]
			d['label_size'][n] = (w, h)

	else:
		for n in d['label']:
			link_text[n]= d['label'][n]
			lab = d['label'][n].replace(' [cytoplasm]', '')
			d['label'][n] = lab
			w = get_label_width(lab, font, font_size)
			h = font_size + 2*padding_labels[1]
			d['label_size'][n] = (w, h)

	# get paths
	paths_svg, d['pos'] = get_paths(
		min_path_length = font_size, 
		max_bend = 0.5*max_w, 
		direction_default = defdir, 
		prevent_overlap = not(overlap), 
		cofactors= cofactors, 
		reverse_cofactor_direction = reverse_cof,
		**d) 

	labels = {}
	for n in d['nodes']:
		if d['node_type'][n] == 'species':
			labels[n] = {}
			labels[n]['x'] = d['pos'][n][0] - 0.5*(d['label_size'][n][0]-2*padding_labels[0])
			labels[n]['y'] = d['pos'][n][1] + 0.25 * font_size
			labels[n]['text'] = d['label'][n]
			# labels[n]['text'] = link_text[n].replace(' [cytoplasm]', '')
			labels[n]['link_text'] = link_text[n]
			labels[n]['font_size'] = fs[n]
	if cofactors:
		for r in cofactors:
			for s in cofactors[r]:
				n = str(s)+str(r)
				labels[n] = {}
				labels[n]['x'] = d['pos'][n][0] - 0.5*(d['label_size'][s][0]-2*padding_labels[0])
				labels[n]['y'] = d['pos'][n][1] + 0.25 * font_size
				labels[n]['text'] = d['label'][s]
				labels[n]['text'] = link_text[s].replace(' [cytoplasm]', '')
				labels[n]['link_text'] = link_text[s]
				labels[n]['font_size'] = fs[s]

	pstyle = 'fill:none;stroke:#cccccc;stroke-width:{};fill-opacity:0;stroke-linecap:round;marker-end:url(#{})'
	paths = {}
	for e in paths_svg:
		paths[e]= {'d': paths_svg[e], 'style':pstyle.format(0.2*font_size, role[e])}

	rn = {}
	for n in d['nodes']:
		if d['node_type'][n] == 'reaction':
			rn[n] = {'x':d['pos'][n][0], 'y':d['pos'][n][1]}
	
	# determine height and width of document
	max_x = max([p[0] for p in d['pos'].values()])
	max_y = max([p[1] for p in d['pos'].values()])
	
	data = {'rxn_nodes': rn, 
			'paths': paths, 
			'labels': labels, 
			'font_size':font_size, 
			'font_family': font.getname()[0],
			'height': max_y + font_size*2,
			'width': max_x + max_w}
	return data


def get_svgdoc(paths, rxn_nodes, labels, font_size, font_family, height, width, easy_edit = True, vonda_compatible = False):
	""" Assemble svg document (easy_edit version, easily editable in inkscape)"""

	# create svg document
	doc = svg()
	doc.set_height(height)
	doc.set_width(width)
	# add marker defenitions
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
	doc.addElement(defs)

	# add paths
	for p in paths:
		pth = path(paths[p]['d'])
		if vonda_compatible:
			pth.set_id(str(p[0]))
		elif easy_edit:
			pth.set_id('path_{}{}'.format(*p))
		else:
			pth.set_id(str(p))
		pth.set_style(paths[p]['style'])
		doc.addElement(pth)

	# add labels
	for n in labels:
		x = labels[n]['x']
		y = labels[n]['y']
		txt = text(labels[n]['text'], x, y)
		txt.set_id(str(n))
		txt.set_font_family(font_family)
		txt.set_font_size(font_size)
		doc.addElement(txt)

	# add reaction nodes
	for c in rxn_nodes:
		x = rxn_nodes[c]['x']
		y = rxn_nodes[c]['y']
		crc = circle(x, y , 5)
		crc.set_id(c)
		crc.set_style('fill:#999999')
		doc.addElement(crc)

		# add reaction value placeholders
		if easy_edit:
			txt = text('-0.0e+00', x + 10, y + 4)
			txt.set_id('rval_' + c)
			txt.set_font_family(font_family)
			txt.set_font_size(8)
			doc.addElement(txt)

	return doc

def main(args):
	file_name = ' '.join(args.json_file)
	with open(file_name) as json_data:
		svg_data = json.load(json_data)
	svgdoc = get_svgdoc(**svg_data)
	svgdoc.save(args.svg_name)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('json_file', metavar= 'file_name.json', nargs= '+')
	parser.add_argument('--svg_name', '-o', default = 'temp.svg')
	args = parser.parse_args()
	main(args)
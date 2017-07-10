import sys
import argparse
import networkx as nx
from PIL import ImageFont
from svg_assembly import get_svgdata, get_svgdoc
from readers import read_graph, get_cofactors_from_sbml

def compatible_graph(graph):
	compatible = True
	unknown_node_types = set()
	try:
		for n in graph.nodes():
			graph.node[n]['x']
			graph.node[n]['y']
			if not graph.node[n]['node_type'] in ['reaction', 'species', 'ctrl']:
				unknown_node_types.add(graph.node[n]['node_type'])
	except KeyError as e:
		print 'KeyError:', e
		if e.args[0]=='node_type':
			print "please add 'node_type' attribute to the nodes in graph (set to 'reaction', 'species' or 'ctrl')"
		compatible = False

	if len(unknown_node_types)>0:
		for nt in unknown_node_types:
			print "WARNING: unrecognized nodetype {}".format(nt)

	return compatible

def main(args):
	
	file_name = ' '.join(args.graph_file)
	# get file with layout data
	ext = file_name.split('.')[-1]
	if ext == 'graphml':
		graph = nx.read_graphml(file_name)
	elif ext == 'gml':
		graph = nx.read_gml(file_name)
	elif ext == 'gexf':
		graph = nx.read_gexf(file_name)
	elif ext == 'g6':
		graph = nx.read_graph6(file_name)
	elif ext == 's6':
		graph = nx.read_sparse6(file_name)
	elif ext == 'gpickle' or ext == 'p':
		graph = nx.read_gpickle(file_name)
	elif ext == 'yaml':
		graph = nx.read_yaml(file_name)
	else:
		print "Graph file format not supported. Supported fileformats: graphml, gexf, gml, g6, s6, gpickle, yaml"

	for n in graph.nodes():
		# from nicholas' tulip output
		if not 'x' in graph.node[n].keys():
			if 'graphics' in graph.node[n].keys():
				g = graph.node[n].pop('graphics')
				graph.node[n]['x'] = g['x']
				graph.node[n]['y'] = g['y']
				if 'h' in g.keys() and not 'node_type' in graph.node[n].keys():
					if g['h'] == 1:
						graph.node[n]['node_type'] = 'reaction'
					elif g['h'] == 2.5:
						graph.node[n]['node_type'] = 'species'
					else:
						graph.node[n]['node_type'] = 'ignore'
		if not 'label' in graph.node[n].keys():
			graph.node[n]['label'] = n

	if compatible_graph(graph):
		d = read_graph(graph)

		font = ImageFont.truetype(args.font_file, 1000)

		if args.add_cofactors_from_sbml:
			sbml_file = ' '.join(args.add_cofactors_from_sbml)
			cofactors = get_cofactors_from_sbml(d, sbml_file)
			for r in cofactors:
				for s in cofactors[r]:
					d['edge_type'][(r,s)]=cofactors[r][s]['role']
		else:
			cofactors = None
		
		svg_data = get_svgdata(
			d= d,
			font=font, 
			font_size= args.font_size, 
			scale=args.scale,
			padding=args.padding,
			padding_labels= args.padding_labels,
			normalize=args.normalize,
			overlap=args.overlap,
			defdir = args.r_direction,
			cofactors = cofactors,
			reverse_cof = args.reverse_cof)
		
		doc = get_svgdoc(**svg_data)
		doc.save(args.svg_name)
		print 'output svg saved in', args.svg_name


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('graph_file', metavar = 'file_name.graphml', nargs = '+')
	parser.add_argument('--add_cofactors_from_sbml', '-cof', metavar = 'model.xml', nargs = '+', help = "Add the omitted cofactors from this sbml (3fbc) model to the graph.")
	parser.add_argument('--svg_name', '-o', default = 'temp.svg', help = "The name/path of the output svg.")
	parser.add_argument('--scale', '-s', type = float, nargs='+', default = [1.0, 1.0], metavar= '1.0', help = "Scale up the graph with this factor. Example: -s 10.0 (10 in both x- and y-direction) Example: -s 20 10 (20 in x-direction, 10 in y-direction")
	parser.add_argument('--padding', type = float, nargs='+', default = [0.0, 0.0], metavar= '0', help = "Extra space (pixels) added to the edges of the svg, e.g. so that all labels are visible in a browser.")
	parser.add_argument('--padding_labels', type = float, nargs = '+', default = [10.0, 10.0], metavar= '10', help = "Space (pixels) around the text of the labels. Can also accept two terms, for x and y-direction.")
	parser.add_argument('--font_file', default = 'C:\Users\User\Documents\Raleway\Raleway-Regular.ttf', metavar = 'C:\Users\User\Documents\Raleway\Raleway-Regular.ttf', help= "The font to be used. Raleway can be downloaded from https://github.com/google/fonts/blob/master/ofl/raleway/Raleway-Regular.ttf")
	parser.add_argument('--font_size', default = 10.0, metavar = '10.0', help = "Font size of the labels.")
	parser.add_argument('--reverse_cof', nargs = '+', default = [], metavar = 'R_0001 R_0020 R_0033', help = "List of reaction ids for which the cofactors must be placed in the opposite way to the default (default is substrates top, products bottom).")
	parser.add_argument('--normalize', dest = 'normalize', action = 'store_true', help = "Translate all coordinates to positive coordinates, so the svg can be viewed in a browser.") 
	parser.add_argument('--overlap', dest = 'overlap', action = 'store_true', help = "Allow for overlap of the reaction arrow paths and metabolite labels.")
	parser.add_argument('--auto_direction', dest = 'r_direction', action = 'store_false', help = "Determine the direction (horizontal/vertical) of the arrow through reaction nodes automatically. (default is all reactions vertical)")
	parser.set_defaults(normalize = False)
	parser.set_defaults(overlap = False)
	parser.set_defaults(r_direction = 'vertical')
	args = parser.parse_args()
	main(args)
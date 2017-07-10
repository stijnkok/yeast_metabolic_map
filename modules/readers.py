import json
import sys
import re
import networkx as nx
import cbmpy as cbm


def read_json_data(data):
	"""
	unpack & repack json dictionary from Nicholas's output
	"""

	graph = data['graph']

	# dictionary of labels
	label = graph['properties']['name']['nodesValues']

	# dictionary of pathways
	pathway = graph['properties']['pathway']['nodesValues']

	# dictionary of node types
	n_type = {}
	for node, node_type in graph['properties']['type']['nodesValues'].items():
		n_type[node] = {'1': 'species', '2':'reaction'}[node_type]

	# list of nodes
	nodes = n_type.keys()

	# nodes that are for layout purposes only (don't represent reaction or species)
	if 'edgesValues' in graph['properties']['viewLayout'].keys():
		helper_nodes = graph['properties']['viewLayout']['edgesValues']
	else:
		helper_nodes = {}
	for e_id in helper_nodes:
		n = helper_nodes[e_id]
		l = []
		for m in n.split('), ('):
			l.append( [float(k) for k in m.strip('(').strip(')').split(',')[:2]])
		helper_nodes[e_id] = l
	
	# list of edges
	edges = []
	
	# dictionary with edge types
	e_type = {}
	extra_nodes = {}
	for i in range(len(graph['edges'])):
		e = graph['edges'][i]
		# convert integers to strings, so i can use them directly as a key later.
		n0 = str(e[0])
		n1 = str(e[1])

		if n_type[n0]=='reaction':
			edges.append((n0, n1))
			if str(i) in helper_nodes.keys():
				extra_nodes[(n0, n1)] = helper_nodes[str(i)]
			else:
				extra_nodes[(n0, n1)] = []
			e_type[(n0, n1)] = 'product'
		else:
			# make reaction node always the source node.
			edges.append((n1, n0))
			if str(i) in helper_nodes.keys():
				extra_nodes[(n1, n0)] = helper_nodes[str(i)][::-1]
			else:
				extra_nodes[(n1, n0)] = []
			e_type[(n1, n0)] = 'substrate'

	# dictionary of (x,y) positions
	pos ={}
	for key, val in graph['properties']['viewLayout']['nodesValues'].items():
		# convert strings with coordinats to tuples with floats.
		val = val.strip('(').strip(')').split(',')
		pos[key] = (float(val[0]), float(val[1]))

	# change node_ids to cmod ids
	cmod_id = {}
	copies = {}
	for n, lab in graph['properties']['viewLabel']['nodesValues'].items():
		if lab in copies:
			copies[lab].append(n)
			lab += '_copy_'+str(len(copies[lab])-1)
		else:
			copies[lab] = [n]
		cmod_id[n] = lab

	# return a dictionary with all the extracted information
	d = {'edges': [], 'nodes': [], 'node_type': {}, 'edge_type':{}, 'extra_nodes':{}, 'pos':{}, 'label': {}, 'pathway': {}}
	for n in nodes:
		new_n = cmod_id[n]
		d['nodes'].append(new_n)
		d['node_type'][new_n] = n_type[n]
		d['pos'][new_n] = pos[n]
		d['label'][new_n] = label[n]

	for n in nodes:
		new_n = cmod_id[n]
		try:
			d['pathway'][new_n] = pathway[n]
		except KeyError:
			d['pathway'][new_n] = "Undefined"
	for e in edges:
		new_e = (cmod_id[e[0]], cmod_id[e[1]])
		d['edges'].append(new_e)
		d['edge_type'][new_e] = e_type[e]
		d['extra_nodes'][new_e] = extra_nodes[e]
	return d


def read_graph(graph):
	label = {}
	node_type = {}
	edge_type = {}
	extra_nodes = {}
	pos = {}
	pathway = {}
	edges = []
	nodes = [n for n in graph.nodes() if graph.node[n]['node_type'] != 'ctrl']
	helper_nodes = [n for n in graph.nodes() if graph.node[n]['node_type'] == 'ctrl']

	# get real edges, ignore helper nodes
	for e in graph.edges():
		hn = []
		if graph.node[e[0]]['node_type'] == 'reaction':
			n0 = e[0] # reaction
			n1 = e[1] # species (or helper node)
			if n1 in helper_nodes:
				while graph.node[n1]['node_type'] == 'ctrl':
					hn.append(n1)
					n1 = graph.successors(n1)[0]
			edge_type[(n0, n1)] =  'product'
			edges.append((n0, n1))
			extra_nodes[(n0, n1)] = []
			for n in hn:
				extra_nodes[(n0, n1)].append((graph.node[n]['x'], graph.node[n]['y']))
		
		elif graph.node[e[0]]['node_type'] == 'species':
			n0 = e[1] # reaction (or helper node)
			n1 = e[0] # species
			if n0 in helper_nodes:
				while graph.node[n0]['node_type'] == 'ctrl':
					hn.append(n0)
					n0 = graph.successors(n0)[0]
			edge_type[(n0, n1)] =  'substrate'
			edges.append((n0, n1))
			extra_nodes[(n0, n1)] = []
			for n in hn[::-1]:
				extra_nodes[(n0, n1)].append((graph.node[n]['x'], graph.node[n]['y']))


	for n in nodes:
		label[n] = graph.node[n]['label']
		node_type[n] = graph.node[n]['node_type']
		pathway[n] = graph.node[n]['pathway']
		pos[n]= (graph.node[n]['x'], graph.node[n]['y'])

	d =  {'edges': edges, 'nodes': nodes, 'node_type': node_type, 'edge_type':edge_type, 'extra_nodes':extra_nodes, 'pos':pos, 'label': label, 'pathway': pathway}
	return d


def get_cofactors_from_sbml(d, sbml_file):
	alt_lab = {'carbon dioxide [cytoplasm]': 'CO2',
			'phosphate [cytoplasm]': 'Pi',
			'diphosphate [cytoplasm]': 'PPi',
			'ammonium [cytoplasm]': 'NH4+'}

	model = cbm.CBRead.readSBML3FBC(sbml_file)

	reagents = {}
	for e in d['edges']: reagents[e[0]] = set()
	for e in d['edges']: reagents[e[0]].add(re.sub('_copy_[0-9]+', '', e[1]))

	cofactors = {}
	for r in reagents: cofactors[r] = {}
	for r in reagents:
		R = model.getReaction(re.sub('_copy_[0-9]+', '', r))
		cof_ids = set(R.getSpeciesIds()).difference(reagents[r])
		for s in cof_ids:
			role = R.getReagentWithSpeciesRef(s).role
			lab = model.getSpecies(s).name
			if lab in alt_lab:
				lab = alt_lab[lab]
			cofactors[r][s] = {'role': role, 'label': lab}

	return cofactors
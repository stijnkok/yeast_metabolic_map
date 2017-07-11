import json
import sys
import argparse
import networkx as nx
from readers import read_json_data

def infodict_to_graph(d, scale = [1,1], padding = [0,0], normalize = False):
	"""
	Input  - Dictionary with layout info.
	Output - networkx.DiGraph object.
	"""

	## adjust coordinates ##
	if len(scale)==1:
		scale.append(scale[0])
	if len(padding) == 1:
		padding.append(padding[0])

	# get minimum x and maximum y coordinate for normalization
	min_x = min([p[0] for p in d['pos'].values()])
	max_y = max([p[1] for p in d['pos'].values()])

	for n, p in d['pos'].items():
		if normalize:
			p = (p[0]-min_x, p[1]-max_y)           # normalize
		p = (p[0]*scale[0], p[1]*scale[1])         # scale
		p = (p[0] + padding[0], p[1] - padding[1]) # add padding
		d['pos'][n] = p

	for e in d['extra_nodes']:
		for i in range(len(d['extra_nodes'][e])):
			p = d['extra_nodes'][e][i]
			if normalize:
				p = (p[0]-min_x, p[1]-max_y)           # normalize
			p = (p[0]*scale[0], p[1]*scale[1])         # scale
			p = (p[0] + padding[0], p[1] - padding[1]) # add padding
			d['extra_nodes'][e][i] = p

	## make graph ##
	graph = nx.DiGraph()

	# add edges and nodes
	for e in d['edges']:
		if d['edge_type'][e] == 'product':
			# include helper/control nodes
			hn_ids = ['hn'+str(i)+str(e) for i in range(len(d['extra_nodes'][e]))]
			for i in range(len(hn_ids)):
				d['pos'][hn_ids[i]] = d['extra_nodes'][e][i]
				d['label'][hn_ids[i]] = 'ctrl'
				d['node_type'][hn_ids[i]] = 'ctrl'
				graph.add_node(hn_ids[i])
			node_id_lst = [e[0]] + hn_ids + [e[1]]
			# add edges to graph (this also adds the nodes)
			for i in range(len(node_id_lst)-1):
				graph.add_edge(node_id_lst[i], node_id_lst[i+1])
		else:
			# include helper/control nodes
			hn_ids = ['hn'+str(i)+str(e) for i in range(len(d['extra_nodes'][e]))]
			for i in range(len(hn_ids)):
				d['pos'][hn_ids[i]] = d['extra_nodes'][e][i]
				d['label'][hn_ids[i]] = 'ctrl'
				d['node_type'][hn_ids[i]] = 'ctrl'
			node_id_lst = [e[1]] + hn_ids[::-1] + [e[0]]
			# add edges to graph (this also adds the nodes)
			for i in range(len(node_id_lst)-1):
				graph.add_edge(node_id_lst[i], node_id_lst[i+1])

	# add node attributes
	for n in graph.nodes():
		graph.node[n]['x'] = d['pos'][n][0]
		graph.node[n]['y'] = d['pos'][n][1]
		graph.node[n]['label'] = d['label'][n]
		graph.node[n]['node_type'] = d['node_type'][n]
		graph.node[n]['pathway'] = d['pathway'][n]
		if d['node_type'][n] == 'ctrl':
			graph.node[n]['size'] = 1.0
		elif d['node_type'][n] == 'reaction':
			graph.node[n]['size'] = 1.5
		else:
			graph.node[n]['size'] = 2.5
	
	return graph

def main(args):
	file_name = ' '.join(args.json_file)
	with open(file_name) as json_data:
		data = json.load(json_data)
	d = read_json_data(data)
	graph = infodict_to_graph(d, args.scale, args.padding, args.normalize)
	nx.write_graphml(graph, args.graph_name)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('json_file', metavar= 'file_name.json', nargs= '+')
	parser.add_argument('--graph_name', '-o', default = 'temp.graphml')
	parser.add_argument('--scale', '-s', type = float, nargs='+', default = [20.0, 20.0], metavar= '20.0', help = "Scale up the graph with this factor. Example: -s 10.0 (10 in both x- and y-direction) Example: -s 20 10 (20 in x-direction, 10 in y-direction")
	parser.add_argument('--padding', type = float, nargs='+', default = [20.0, 20.0], metavar= '20', help = "Extra space (pixels) added to the edges of the graph, e.g. so that all labels are visible in a browser when converted to svg.")
	parser.add_argument('--normalize', dest = 'normalize', action='store_true', help = "Translate all coordinates to positive coordinates, so if converted to svg, it can be viewed in a browser.") 
	args = parser.parse_args()
	main(args)
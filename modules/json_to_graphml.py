import json
import sys
import argparse
import networkx as nx
from readers import read_json_data

def infodict_to_graph(d):
	graph = nx.DiGraph()

	for e in d['edges']:
		if d['edge_type'][e] == 'product':
			hn_ids = ['hn'+str(i)+str(e) for i in range(len(d['extra_nodes'][e]))]
			for i in range(len(hn_ids)):
				d['pos'][hn_ids[i]] = d['extra_nodes'][e][i]
				d['label'][hn_ids[i]] = 'ctrl'
				d['node_type'][hn_ids[i]] = 'ctrl'
				graph.add_node(hn_ids[i])
			node_id_lst = [e[0]] + hn_ids + [e[1]]
			for i in range(len(node_id_lst)-1):
				graph.add_edge(node_id_lst[i], node_id_lst[i+1])
		else:
			hn_ids = ['hn'+str(i)+str(e) for i in range(len(d['extra_nodes'][e]))]
			for i in range(len(hn_ids)):
				d['pos'][hn_ids[i]] = d['extra_nodes'][e][i]
				d['label'][hn_ids[i]] = 'ctrl'
				d['node_type'][hn_ids[i]] = 'ctrl'
			node_id_lst = [e[1]] + hn_ids[::-1] + [e[0]]
			for i in range(len(node_id_lst)-1):
				graph.add_edge(node_id_lst[i], node_id_lst[i+1])
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
	graph = infodict_to_graph(d)
	nx.write_graphml(graph, args.graph_name)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('json_file', metavar= 'file_name.json', nargs= '+')
	parser.add_argument('--graph_name', '-o', default = 'temp.graphml')
	args = parser.parse_args()
	main(args)
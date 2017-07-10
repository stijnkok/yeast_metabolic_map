import sys
import argparse
import networkx as nx

def snap_to_grid(graph, grid, offset):
	""" Align x,y position of nodes in the graph with a grid """
	G = nx.DiGraph()
	G.add_nodes_from(graph.nodes())
	G.add_edges_from(graph.edges())
	for n in graph.nodes():
		x = graph.node[n]['x']
		y = graph.node[n]['y']
		x = offset[0] + round((x-offset[0])/grid[0])*grid[0]
		y = offset[1] + round((y-offset[1])/grid[1])*grid[1]
		G.node[n]['x'] = x
		G.node[n]['y'] = y
		G.node[n]['label'] = graph.node[n]['label']
		G.node[n]['node_type'] = graph.node[n]['node_type']
	return G

def main(args):
	file_name = ' '.join(args.graphml_file)
	graph = nx.read_graphml(file_name)
	graph = snap_to_grid(graph, args.grid, args.offset)
	gx, gy = args.grid
	gx, gy = int(gx), int(gy)
	nx.write_graphml(graph, file_name.replace('.graphml', '')+'_grid_{}_{}.graphml'.format(gx, gy))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('graphml_file', metavar= 'file_name.graphml', nargs='+')
	parser.add_argument('--grid', '-g', type = float, nargs='+', default = [1.0, 1.0])
	parser.add_argument('--offset', type = float, nargs= '+', default = [0.0, 0.0])
	args = parser.parse_args()
	main(args)
import json
import sys
import math
import numpy as np
import networkx as nx
from itertools import groupby, combinations, product
from PIL import ImageFont
from svg.path import Path, Line, Arc, CubicBezier

from pysvg.structure import svg, g
from pysvg.text import text
from pysvg.shape import path
from pysvg.core import TextContent


def get_cofactors_layout(rpos, cofactors, direction, label, label_size, dist = (10, 15)):
	path = {'v': 'm {x},{y} v {v} c 0,{dy} {dx1},{dy} {dx2},{dy}', 
	        'h': 'm {x},{y} h {h} c {dx},0 {dx},{dy1} {dx},{dy2}'}
	gr = []
	subs = [m for m in cofactors.keys() if cofactors[m]['role']=='substrate']
	prds = [m for m in cofactors.keys() if cofactors[m]['role']=='product']
	
	while len(subs)>0 and len(prds)>0:
		mtch = []
		for s, p in product(subs, prds):
			lab1 = label[s]
			lab2 = label[p]
			mtch.append([sum([lab1.count(ch) for ch in lab2]), s, p])
		mtch.sort()
		best_mtch = mtch[-1][1:]
		gr.append(best_mtch)
		subs.remove(best_mtch[0])
		prds.remove(best_mtch[1])
	
	gr += [[m, None] for m in subs]
	gr += [[None, m] for m in prds]

	if direction == 'v':
		dx = -dist[0]
		dy = -dist[1]
		add_dy = dy

		for i in range(len(gr)):
			g = gr[i]
			dx = dx * -1
			for j in range(2):
				if g[j]:
					x = rpos[0] + dx + 0.5*label_size[g[j]][0]*dx/abs(dx)
					y = rpos[1] + dy
					cofactors[g[j]]['pos'] = (x, y)
					cofactors[g[j]]['path'] = path['v'].format(
						x = rpos[0], 
						y = rpos[1], 
						v = dy-abs(add_dy)*dy/abs(dy),
						dy = abs(add_dy)*dy/abs(dy), 
						dx1 = dx/2., 
						dx2 = dx)
				dy = dy * -1
			if i%2:
				dy+= add_dy
	else:
		if dist[1] <0:
			dist = (dist[0], abs(dist[1]))
			for i in range(len(gr)):
				gr[i] = gr[i][::-1]

		dx_path = {0:[-dist[1], dist[1]], 1:[-dist[1], dist[1]]}
		dx_pos = {0:[0,0], 1:[0,0]}
		dy = -dist[0]

		for i in range(len(gr)):
			g = gr[i]
			dy = dy * -1
			
			if g[0]: # left (substrate)
				x = rpos[0] + dx_pos[i%2][0] - 0.5*label_size[g[0]][0]
				y = rpos[1] + dy + 0.5*label_size[g[0]][1]*dy/abs(dy)
				cofactors[g[0]]['pos']  = (x, y)
				cofactors[g[0]]['path'] = path['h'].format(
					x = rpos[0],
					y = rpos[1],
					h = dx_pos[i%2][0],
					dx = dx_path[i%2][0], 
					dy1 = dy/2., 
					dy2 = dy)

				dx_pos[i%2][0] += -label_size[g[0]][0]

			if g[1]:# right (product)
				x = rpos[0] + dx_pos[i%2][1] + 0.5*label_size[g[1]][0]
				y = rpos[1] + dy + 0.5*label_size[g[1]][1]*dy/abs(dy)
				cofactors[g[1]]['pos']  = (x, y)
				cofactors[g[1]]['path'] = path['h'].format(
					x = rpos[0],
					y = rpos[1],
					h = dx_pos[i%2][1],
					dx = dx_path[i%2][1], 
					dy1 = dy/2., 
					dy2 = dy)

				dx_pos[i%2][1] += label_size[g[1]][0]
			
	return cofactors

def eucdist(a, b):
	""" Distance between points a & b (x,y tuples)"""
	return math.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)

def least_circular_node(pos, nodes):
	""" Find the node that is not on a circle. Nodes that are all evenly placed 
	on a circle will all have the exact same set of distances to each other. This
	function compares each nodes' set of distances to other nodes to the distances-set 
	of other nodes, and finds the node that is the most different. 
	"""
	
	# get all distances between nodes
	dist= {}

	for n in nodes:
		dist[n] = []

	for p in product(nodes, nodes):
		dist[p[0]].append(round(eucdist(pos[p[0]], pos[p[1]]), 0))
		dist[p[0]].sort()

	common_dist=[]

	for n0 in nodes:
		
		ndist = len(dist[n0]) # number of distances
		shared_dist = []      # list of the fraction of distances shared with each node
		for n1 in nodes:
			# add fraction of distances shared between n0 and n1
			nshared = sum([dist[n0][i] in dist[n1] for i in range(ndist)])
			shared_dist.append(nshared/float(ndist))
		common_dist.append((sum(shared_dist)/float(len(nodes)), n0)) # mean of the fractions

	common_dist.sort()
	return common_dist[0][1], sum([cd[0] for cd in common_dist])/len(common_dist)

def get_nodes_on_circle(pos, nodes):
	""" Find nodes on circle by excluding outlier nodes """
	n, circularity = least_circular_node(pos, nodes)
	while circularity < 0.95:
		nodes.remove(n)
		n, circularity = least_circular_node(pos, nodes)
	return nodes

def intersect_circ_rect(cx, cy, r, x, y, w, h):
	"""
	Get intersection of a circle with center cy,cx and radius r
	with a rectangle with center x,y and width w and height h.
	"""
	# get lines that define the rectangle
	x1, x2, y1, y2 = x - 0.5*w, x + 0.5*w, y + 0.5*h, y - 0.5*h
	# find intersections
	# using r^2 = x^2 + y^2
	intersections = []
	for x_ in [x1, x2]:
		# check if x-line intersects with circle
		if (x_ - cx)**2 < r**2:
			# get y-coordinates of the intersections
			for sgn in [-1,1]:
				y_ =  sgn*math.sqrt(r**2 - (x_ - cx)**2) + cy
				# check if y-coordinate is between y1 and y2
				if y_ >= min(y1, y2) and y_ <= max(y1, y2):
					intersections.append((x_, y_))
	for y_ in [y1, y2]:
		# check if y-line intersects with circle
		if (y_ - cy)**2 < r**2:
			# get x-coordinates of the intersections
			for sgn in [-1,1]:
				x_ =  sgn*math.sqrt(r**2 - (y_ - cy)**2) + cx
				# check if x-coordinate is between x1 and x2
				if x_ >= min(x1, x2) and x_ <= max(x1, x2):
					intersections.append((x_, y_))
	return intersections 

def arc_path(x0, y0, r, swp, x1, y1):
	""" return circle segment svg path """
	return 'M{},{} A{},{} 0 0 {} {},{}'.format(x0,y0, r,r, swp, x1,y1)

def path_node_direction(r_pos, s_pos_lst, default='v'):
	""" 
	Determine direction of path start. All given species coordinates 
	get 1 vote on the direction being horizontal or vertical. If any
	species has the same x or y coordinate, it becomes the only vote.
	If undecided, direction is default (vertical).
	"""
	h = 0
	v = 0
	for x, y in s_pos_lst:
		if r_pos[0]==x:
			h=0
			v=1
			break
		elif r_pos[1]==y:
			h=1
			v=0
			break
		elif abs(r_pos[0] - x) > abs(r_pos[1] - y):
			h+=1
		else:
			v+=1

	if h==v:
		start_curve = default
	elif h>v:
		start_curve = 'h'
	else:
		start_curve = 'v'
	
	return start_curve

def path_direction_end(dx, dy, start_curve, label_size, min_bend = 15, min_length = 10):
	"""
	Determine direction of path end based on position of node relative to reaction (dx & dy),
	the direction of the path at the start and the size of the label.
	"""
	width, height = label_size

	if start_curve == 'v':
		if abs(dx) > min_bend + 0.5*width:
			end_curve = 'h'
		else:
			assert abs(dy) > height + min_length, 'species label and reaction are overlapping or too close'
			end_curve = 'v'
	
	elif start_curve == 'h':
		if abs(dy) > min_bend + 0.5*height:
			end_curve = 'v'
		else:
			assert abs(dx) > min_length + 0.5*width, 'species label and reaction are overlapping or too close'
			end_curve = 'h'
	
	return end_curve

def path_position_end(pos, dx, dy, end_curve, label_size, padding = (0,0)):
	"""
	Get adjusted path end position based on the label position (mid point),
	label size and direction of path end.
	"""
	width, height = label_size
	if end_curve == 'v':
		x = pos[0]
		y = pos[1] - 0.5*height*dy/abs(dy)
	elif end_curve == 'h':
		x = pos[0] - 0.5*width*dx/abs(dx)
		y = pos[1]
	
	return [x, y]

def adjust_duplicate_path_ends(s, p, spos, p_nodes):
	"""
	Adjust duplicate path end points.
	"""

	# group path end nodes by position
	gr = (([],[]),([],[]))
	for r in p:
		if p[r][0] == spos[0]:
			if p[r][1] < spos[1]:
				# adjust x top
				gr[0][0].append(r)
			else:
				# adjust x bottom
				gr[0][1].append(r)
		
		elif p[r][1] == spos[1]:
			if p[r][0] < spos[0]:
				# adjust y left
				gr[1][0].append(r)
			else:
				# adjust y right
				gr[1][1].append(r)
	
	# adjust positions if there are duplicate path end points
	for i in range(2):
		for j in range(2):
			for r in gr[i][j]:
				# find path ends that should be in middle, so 
				# that path stays a straight line instead of 
				# becoming s-shaped.
				if p_nodes[(r,s)][-2][i] == spos[i]:
					mid = True
					mid_r = r
					break
			else:
				mid = False

			if mid:
				gr[i][j].remove(mid_r) # don't move this position
				p[mid_r][i] = spos[i]

			# sort r_nodes in group based on r_node x or y positions
			gr[i][j].sort(key = lambda r: p_nodes[(r,s)][-2][i])

			if mid and not len(gr[i][j])%2:
				# even:
				# place half the path ends on one side of mid,
				# place the other half on the other side.
				gr1 = gr[i][j][:len(gr[i][j])/2][::-1]
				gr2 = gr[i][j][len(gr[i][j])/2:]
				for r in gr1:
					p[r][i] = spos[i] -10 - 10*gr1.index(r)
				for r in gr2:
					p[r][i] = spos[i] +10 + 10*gr2.index(r)

			elif mid and len(gr[i][j])%2:
				# odd
				for r in gr[i][j]:
					# adjust x or y position
					p[r][i] = spos[i] + (len(gr[i][j])*-5+5) + 10*gr[i][j].index(r)
				# decide which side of mid has the majority of path ends
				# based on reaction position.
				num_side1 = len([r for r in gr[i][j] if p_nodes[(r,s)][-2][i] < spos[i]])
				num_side2 = len([r for r in gr[i][j] if p_nodes[(r,s)][-2][i] > spos[i]])
				if num_side1 > num_side2:
					gr1 = [r for r in gr[i][j] if p[r][i] <= spos[i]][::-1]
					gr2 = [r for r in gr[i][j] if p[r][i] >  spos[i]]
				else:
					gr1 = [r for r in gr[i][j] if p[r][i] <  spos[i]][::-1]
					gr2 = [r for r in gr[i][j] if p[r][i] >= spos[i]]
				# place path ends on each side
				for r in gr1:
					p[r][i] = spos[i] -10 - 10*gr1.index(r)
				for r in gr2:
					p[r][i] = spos[i] +10 + 10*gr2.index(r)				
			else:
				for r in gr[i][j]:
					# adjust x or y position
					p[r][i] = spos[i] + (len(gr[i][j])*-5+5) + 10*gr[i][j].index(r)
			
	return p

def overlapping(points, label_pos, label_size):
	"""
	Check if any point in list of points overlaps with a label.
	"""

	l_x = label_pos[0]-0.5*label_size[0]
	l_y = label_pos[1]-0.5*label_size[1]

	for point in points:
		if point.real > l_x and point.real < l_x + label_size[0] and point.imag > l_y and point.imag < l_y + label_size[1]:
			overlap = True
			break
	else:
		overlap = False
	return overlap

def get_path_segments(start, end, start_direction, end_direction, max_bend= 40, adjust= None):
	"""get an svg.Path() object from 'start' to 'end'. 
	directionality is 'v' (vertical) or 'h' (horizontal).
	resulting path type is something like:
	R------------S (l-shaped)

	R-----
		  \
		   \
			-----S (s-shaped)
	or 
	R-------
			\
			 |
			 |
			 S (j-shaped)
	depending on given path directionality of reaction and species. """

	start_direction = start_direction[0].lower()
	end_direction = end_direction[0].lower()
	shape = start_direction + end_direction
	dx = end.real - start.real
	dy = end.imag - start.imag

	if start.real==end.real or start.imag==end.imag:
		segs = [Line(start, end)]

	elif shape == 'vh':
		# j-shaped
		if adjust == None:
			adjust = [0,0]
		segs = path_segments_vh(start, end, max_bend, adjust)
	
	elif shape == 'hv':
		# j-shaped
		if adjust == None:
			adjust = [0,0]
		segs = path_segments_hv(start, end, max_bend, adjust)

	elif shape == 'vv':
		# s-shaped
		if adjust == None:
			adjust = [0, 0]
			midpoint = [start.real + 0.5*dx, start.imag + 0.5*dy]
		elif len(adjust) == 2:
			midpoint = [start.real + 0.5*dx, start.imag + 0.5*dy + adjust[0]*dy/abs(dy)]
		else:
			midpoint = [start.real + 0.5*dx, start.imag + (adjust[2] + adjust[0])*dy/abs(dy)]
		# same as vh -- hv
		seg1 = path_segments_vh(
			start,
			complex(*midpoint),
			max_bend, 
			adj = [adjust[0], 0])
		seg2 = path_segments_hv(
			complex(*midpoint),
			end,
			max_bend,
			adj = [0, adjust[1]])
		segs = seg1+seg2

	elif shape == 'hh':
		# s-shaped
		if adjust == None:
			adjust = [0, 0]
			midpoint = [start.real + 0.5*dx, start.imag + 0.5*dy]
		elif len(adjust) == 2:
			midpoint = [start.real + 0.5*dx + adjust[0]*dx/abs(dx), start.imag + 0.5*dy]
		else:
			midpoint = [start.real + (adjust[2] + adjust[0])*dx/abs(dx), start.imag + 0.5*dy]
		# same as hv -- vh
		seg1 = path_segments_hv(
			start,
			complex(*midpoint),
			max_bend,
			adj = [adjust[0], 0])
		seg2 = path_segments_vh(
			complex(*midpoint),
			end,
			max_bend,
			adj = [0, adjust[1]])
		segs = seg1+seg2

	return segs

def path_segments_vh(start, end, max_bend, adj):
	# j-shaped path 'vertical' to 'horizontal'
	dx = end.real - start.real
	dy = end.imag - start.imag
	dv = adj[0] * dy/abs(dy)
	dh = adj[1] * dx/abs(dx)

	v = Line(start, complex(start.real, start.imag+dv))
	h = Line(complex(end.real-dh, end.imag), end)


	p1 = v.end + complex(0, 0.75*(dy-dv))
	p2 = v.end + complex(0.25*(dx-dh), dy-dv)
	p3 = v.end + complex(dx-dh, dy-dv)
	
	if abs(dx) - abs(dh) > max_bend:
		dh = dx - max_bend*dx/abs(dx)
		h = Line(complex(end.real-dh, end.imag), end)
		p2 = complex(v.end.real + 0.25*max_bend*dx/abs(dx), p2.imag)
		p3 = complex(v.end.real + max_bend*dx/abs(dx), p3.imag)
	if abs(dy) - abs(dv) > max_bend:
		dv = dy - max_bend*dy/abs(dy)
		v = Line(start, complex(start.real, start.imag+dv))
		p1 = complex(p1.real, v.end.imag + 0.75*max_bend*dy/abs(dy))
		p2 = complex(p2.real, v.end.imag + max_bend*dy/abs(dy))
		p3 = complex(p3.real, v.end.imag + max_bend*dy/abs(dy))

	c = CubicBezier(complex(start.real, start.imag+dv), p1, p2, p3)
	
	segs = [v,c,h]
	if v.start == v.end:
		segs.remove(v)
	if h.start == h.end:
		segs.remove(h)
	
	return segs

def path_segments_hv(start, end, max_bend, adj):
	# j-shaped path 'horizontal' to 'vertical' 
	dx = end.real - start.real
	dy = end.imag - start.imag
	dv = adj[0] * dy/abs(dy)
	dh = adj[1] * dx/abs(dx)

	h = Line(start, complex(start.real+dh, start.imag))
	v = Line(complex(end.real, end.imag-dv), end)

	p1 = h.end + complex(0.75*(dx-dh), 0)
	p2 = h.end + complex(dx-dh, 0.25*(dy-dv))
	p3 = h.end + complex(dx-dh, dy-dv)

	if abs(dx) > max_bend:
		dh = dx - max_bend*dx/abs(dx)
		h = Line(start, complex(start.real+dh, start.imag))
		p1 = complex(h.end.real + 0.75*max_bend*dx/abs(dx), p1.imag)
		p2 = complex(h.end.real + max_bend*dx/abs(dx), p2.imag)
		p3 = complex(h.end.real + max_bend*dx/abs(dx), p3.imag)
	if abs(dy) > max_bend:
		dv = dy - max_bend*dy/abs(dy)
		v = Line(complex(end.real, end.imag-dv), end)
		p2 = complex(p2.real, h.end.imag + 0.25*max_bend*dy/abs(dy))
		p3 = complex(p3.real, h.end.imag + max_bend*dy/abs(dy))

	c = CubicBezier(complex(start.real+dh, start.imag), p1, p2, p3)

	segs = [h,c,v]
	if v.start == v.end:
		segs.remove(v)
	if h.start == h.end:
		segs.remove(h)

	return segs
	
def get_paths(edges, nodes, node_type, extra_nodes, pos, label, label_size, pathway, cofactors = None, min_path_length = 10, max_bend = 40, prevent_overlap = True, direction_default = 'v', reverse_cofactor_direction = []):
	"""
	edges: list of edges
	nodes: list of nodes
	node_type: dictionary of node types ('reaction' or 'species')
	extra_nodes: dictionary of additional layout nodes for edges (list of [x,y] positions)
	pos: dictionary with node positions ([x,y] tuples)
	label: dictionary with node labels
	label_size: dictionary with node label size ([width, height] tuples)
	"""

	if direction_default:
		direction_default = direction_default[0].lower()
	s_nodes = set()
	r_nodes = []
	for n in nodes:
		if node_type[n] == 'reaction':
			r_nodes.append(n)
		else:
			s_nodes.add(n)

	# get circular paths, exclude edges with circular paths
	arc_paths = {} # svg arcs

	cyclic_pathways = set([pw for pw in pathway.values() if 'cycle' in pw.lower()])
	for pw in cyclic_pathways:
		circle_nodes = [n for n in nodes if pathway[n] == pw]

		# get nodes that are placed on a circle
		circle_nodes = get_nodes_on_circle(pos, circle_nodes)
		
		# get center coordinates of circle
		cx = sum([pos[n][0] for n in circle_nodes])/len(circle_nodes)
		cy = sum([pos[n][1] for n in circle_nodes])/len(circle_nodes)
		# potential egdes
		# find closest node pairs for each source node
		sourcenodes = [n for n in circle_nodes if node_type[n]=='reaction']
		potential_edges = []
		for sn in sourcenodes:
			# get distances for all nodes on circle
			dist = [eucdist((pos[sn][0], pos[sn][1]), (pos[n][0], pos[n][1])) for n in circle_nodes]
			# sort nodes based on distance
			sorted_n = [n for (dst,n) in sorted(zip(dist, circle_nodes))] 
			# get edges with closest nodes
			potential_edges += [(sn, tn) for tn in sorted_n[1:]][:2]

		circle_edges = [e for e in edges if e in potential_edges]
		for e in circle_edges:
			x0, y0 = pos[e[0]][0], pos[e[0]][1]
			x, y = pos[e[1]][0], pos[e[1]][1]
			w, h = label_size[e[1]]
			r = math.sqrt((x - cx)**2 + (y - cy)**2)
			# get the 2 intersections of the circle with the rectangular label 'bounding box'
			intersections = intersect_circ_rect(cx, cy, r, x, y, w, h)
			# get the intersection closest to reaction node, this will be the path end
			dist = [eucdist((x0, y0), intrsc) for intrsc in intersections]
			intrsc = intersections[dist.index(min(dist))]
			# determine sweep-flag for the svg arc
			angle0 = np.angle(complex((x0-cx),(y0-cy)))
			angle1 = np.angle(complex((intrsc[0]-cx),(intrsc[1]-cy)))
			if angle1-angle0 < 0:
				swp = 0 # anti-clockwise
			else:
				swp = 1 # clockwise
			arc_paths[e] = arc_path(x0, y0, r, swp, *intrsc)
			edges.remove(e)

	v_nodes={}
	for r in r_nodes: v_nodes[r] = [] # voter nodes for reaction direction

	p_nodes={} # path nodes coordinates (by edge)
	p_nodes_flat = [] # list of path nodes coordinates
	direction = {}
	for e in edges:
		p_nodes[e]=[pos[e[0]]] + extra_nodes[e] + [pos[e[1]]]
		p_nodes_flat += p_nodes[e]
		v_nodes[e[0]].append(p_nodes[e][1])
		if 'cycle' in pathway[e[0]].lower() or not direction_default:
			direction[e]=['h']*len(p_nodes[e])
		else:	
			direction[e]=[direction_default]*len(p_nodes[e])

	cof_adj = {}.fromkeys(edges, 0)
	if cofactors:
		for e in edges:
			if len(cofactors[e[0]])>0:
				cof_adj[e] = 1.5*min_path_length
	else:
		cofactors = {}.fromkeys(r_nodes, {})

	path_end = {}
	for e in edges: path_end[e[1]] = {}


	for e in edges:
		# determine curve direction of path at reaction node
		if 'cycle' in pathway[e[0]].lower() or not direction_default:
			direction[e][0] = path_node_direction(pos[e[0]], v_nodes[e[0]])

		# get cofactor positions and paths
		if e[0] in reverse_cofactor_direction:
			dst = (min_path_length, -1.5*min_path_length)
		else:
			dst = (min_path_length, 1.5*min_path_length)

		cofactors[e[0]] = get_cofactors_layout(pos[e[0]], cofactors[e[0]], direction[e][0], label, label_size, dist=dst)
		for s in cofactors[e[0]]:
			s_nodes.add(str(s)+str(e[0]))
			label_size[str(s)+str(e[0])] = label_size[s]
			pos[str(s)+str(e[0])] = cofactors[e[0]][s]['pos']
		# determine curve direction of path at helper nodes
		for i in range(1, len(p_nodes[e])-1):
			n = p_nodes[e][i]
			direction[e][i] = path_node_direction(n, [p_nodes[e][i-1], p_nodes[e][i+1]])
		
		# determine curve direction of path at species node
		dx = pos[e[1]][0] - p_nodes[e][-2][0]
		dy = pos[e[1]][1] - p_nodes[e][-2][1]
		direction[e][-1] = path_direction_end(dx, dy, direction[e][-1], label_size[e[1]], min_length = min_path_length)
		
		if not 'cycle' in pathway[e[0]].lower() and direction_default == 'v' and direction[e][-1] == 'h' and (p_nodes[e][-2][0], pos[e[1]][1]) in p_nodes_flat:
			direction[e][-1] = 'v'

		# determine endpoint of path at species node
		p_nodes[e][-1] = path_position_end(pos[e[1]], dx, dy, direction[e][-1], label_size[e[1]])

		# collect path ends for duplicates adjustment
		path_end[e[1]][e[0]] = p_nodes[e][-1]

	# adjust horizontal path ends (no more than 1 on each side of species label)
	max_hz =1
	for s in path_end:
		left = []
		right = []

		for r in path_end[s]:
			if path_end[s][r][1] == pos[s][1]:
				if path_end[s][r][0] < pos[s][0]:
					left.append(r)
				else:
					right.append(r)
		for rns in [left, right]:
			if len(rns) > max_hz:
				rns.sort(key = lambda r: abs(pos[s][1] - pos[r][1])) # sort on y-distance to reaction node
				for r in rns[1:]:
					e = (r,s)
					dx = pos[s][0] - p_nodes[e][-2][0]
					dy = pos[s][1] - p_nodes[e][-2][1]
					direction[e][-1] = 'v'
					p_nodes[e][-1] = path_position_end(pos[e[1]], dx, dy, direction[e][-1], label_size[e[1]])
					path_end[s][r] = p_nodes[e][-1]

	# adjust duplicate path ends
	for s in path_end:
		path_end[s] = adjust_duplicate_path_ends(s, path_end[s], pos[s], p_nodes)
	for e in edges:
		p_nodes[e][-1] = path_end[e[1]][e[0]]


	# get path segments
	path_segs = {}
	for e in edges:
		path_segs[e]=[]
		for i in range(len(p_nodes[e])-1):
			# get list of path segments
			segs = get_path_segments(
				start = complex(*p_nodes[e][i]), 
				end = complex(*p_nodes[e][i+1]), 
				start_direction = direction[e][i], 
				end_direction = direction[e][i+1],
				max_bend = max_bend,
				adjust = [cof_adj[e],0])
			path_segs[e].append(segs)

	
	if not prevent_overlap:
		# return paths if no overlap prevention
		svg_paths = arc_paths
		for e in edges:
			p = []
			for segs in path_segs[e]:
				p+=segs
			svg_paths[e] = Path(*p).d()
		for r in cofactors:
			for s in cofactors[r]:
				svg_paths[(r,s)] = cofactors[r][s]['path']

		return svg_paths, pos
	
	## adjust path/label overlap ## 
	
	changed_direction = []
	# change path segment basic shape if overlapping with labels
	print 'checking path/label overlap...'
	for e in edges:
		for i in range(len(path_segs[e])):
			seg = path_segs[e][i]

			p = Path(*seg)
			# get 100 (x + yj) points on the path
			points = [p.point(t) for t in [k/100. for k in range(101)]]
			for n in s_nodes.difference({e[1]}):
				# check if any points are inside species labels
				if overlapping(points, pos[n], label_size[n]):
					# change j shape to s shape
					print 'changing path {} shape to prevent {} label overlap...'.format(e, n)
					direction[e][i+1] = direction[e][i]
					new_path = True
					break
			else:
				new_path = False
		
		changed_direction.append(new_path)
		
		# determine new path ends
		dx = pos[e[1]][0] - p_nodes[e][-2][0]
		dy = pos[e[1]][1] - p_nodes[e][-2][1]
		p_nodes[e][-1] = path_position_end(pos[e[1]], dx, dy, direction[e][-1], label_size[e[1]])

		# collect path ends for duplicates adjustment
		path_end[e[1]][e[0]] = p_nodes[e][-1]

	# if not any(changed_direction):
	# 	# return paths if no overlap found
	# 	svg_paths = arc_paths
	# 	for e in edges:
	# 		p = []
	# 		for segs in path_segs[e]:
	# 			p+=segs
	# 		svg_paths[e] = Path(*p).d()
	# 	for r in cofactors:
	# 		for s in cofactors[r]:
	# 			svg_paths[(r,s)] = cofactors[r][s]['path']

	# 	return svg_paths, pos


	# adjust horizontal path ends 
	max_hz =1 # (no more than 1 path end on each side of species label)
	for s in path_end:
		left = []
		right = []
		for r in path_end[s]:
			if path_end[s][r][1] == pos[s][1]:
				if path_end[s][r][0] < pos[s][0]:
					left.append(r)
				else:
					right.append(r)
		for rns in [left, right]:
			if len(rns) > max_hz:
				rns.sort(key = lambda r: abs(pos[s][1] - pos[r][1])) # sort on y-distance to reaction node
				for r in rns[1:]:
					e = (r,s)
					dx = pos[s][0] - p_nodes[e][-2][0]
					dy = pos[s][1] - p_nodes[e][-2][1]
					direction[e][-1] = 'v'
					p_nodes[e][-1] = path_position_end(pos[e[1]], dx, dy, direction[e][-1], label_size[e[1]])
					path_end[s][r] = p_nodes[e][-1]	
	
	# adjust duplicate path ends
	for s in path_end:
		path_end[s] = adjust_duplicate_path_ends(s, path_end[s], pos[s], p_nodes)
	for e in edges:
		p_nodes[e][-1] = path_end[e[1]][e[0]]

	# get new path segments
	vv_shaped=[]
	hh_shaped=[]
	path_segs = {}
	for e in edges:
		path_segs[e]=[]
		for i in range(len(p_nodes[e])-1):
			# get list of path segments
			adj = [cof_adj[e],0]
			start = complex(*p_nodes[e][i])
			end = complex(*p_nodes[e][i+1])
			segs = get_path_segments(
				start,
				end,
				start_direction = direction[e][i], 
				end_direction = direction[e][i+1],
				max_bend = max_bend,
				adjust = adj)
			path_segs[e].append(segs)

			# collect information for finding and adjusting path/path overlap
			if direction[e][i] == direction[e][i+1] == 'v':
				dx = end.real - start.real
				dy = end.imag - start.imag
				y = start.imag + 0.5*dy + adj[0]*dy/abs(dy)
				info = {'y':y, 'dx':dx, 'dy':dy, 'start':start, 'end':end, 'e,i':(e,i)}
				vv_shaped.append(info)
			if direction[e][i] == direction[e][i+1] == 'h':
				dx = end.real - start.real
				dy = end.imag - start.imag
				x = start.real + 0.5*dx + adj[0]*dx/abs(dx)
				info = {'x':x, 'dx':dx, 'dy':dy, 'start':start, 'end':end, 'e,i':(e,i)}
				hh_shaped.append(info)

	print 'adjusting path/path overlap'
	
	vv_shaped.sort(key = lambda k:k['y']) # sort by midpoint y-coordinate
	for y, info in groupby(vv_shaped, key = lambda k:k['y']):
		# determine overlap for all path segments with same midpoint y-coordinate
		segs_to_adjust = []
		for c in combinations(list(info), 2):
			min_x1 = min(c[0]['start'].real, c[0]['end'].real)
			max_x1 = max(c[0]['start'].real, c[0]['end'].real)
			min_x2 = min(c[1]['start'].real, c[1]['end'].real)
			max_x2 = max(c[1]['start'].real, c[1]['end'].real)
			if max(min_x1, min_x2) < min(max_x1, max_x2):
				# paths have overlap!
				if len(segs_to_adjust) == 0:
					segs_to_adjust.append(list(c))
				else:
					for lst in segs_to_adjust:
						if c[0] in lst and not c[1] in lst:
							lst.append(c[1])
							break
						elif c[1] in lst and not c[0] in lst:	
							lst.append(c[0])
							break
					else: # if neither c[0] or c[1] in any lst:
						segs_to_adjust.append(list(c))
		for lst in segs_to_adjust:
			lst.sort(key = lambda k:k['start'].imag)
			lst.sort(key = lambda k:k['start'].real*-k['dx']/abs(k['dx']))
			l = len(lst)
			mp_adj = [(l*-5 +5) + i*10 for i in range(l)] # e.g. [-10, 0, 10] for l = 3
			for n in range(l):
				dy = lst[n]['dy']
				e, i = lst[n]['e,i']
				start = lst[n]['start']
				end = lst[n]['end']
				adj = [cof_adj[e], 0, 0.5*abs(dy) -cof_adj[e] + mp_adj[n]]
				segs = get_path_segments(start,	end, 'v', 'v', max_bend, adj)
				path_segs[e][i] = segs

	hh_shaped.sort(key = lambda k:k['x']) # sort by midpoint y-coordinate
	for x, info in groupby(hh_shaped, key = lambda k:k['x']):
		# determine overlap for all path segments with same midpoint y-coordinate
		segs_to_adjust = []
		for c in combinations(list(info), 2):
			min_y1 = min(c[0]['start'].imag, c[0]['end'].imag)
			max_y1 = max(c[0]['start'].imag, c[0]['end'].imag)
			min_y2 = min(c[1]['start'].imag, c[1]['end'].imag)
			max_y2 = max(c[1]['start'].imag, c[1]['end'].imag)
			if max(min_y1, min_y2) < min(max_y1, max_y2):
				# paths have overlap!
				if len(segs_to_adjust) == 0:
					segs_to_adjust.append(list(c))
				else:
					for lst in segs_to_adjust:
						if c[0] in lst and not c[1] in lst:
							lst.append(c[1])
							break
						elif c[1] in lst and not c[0] in lst:
							lst.append(c[0])
							break
					else: # if neither c[0] or c[1] in any lst:
						segs_to_adjust.append(list(c))
		for lst in segs_to_adjust:
			lst.sort(key = lambda k:k['start'].real)
			l = len(lst)
			mp_adj = [(l*-5 +5) + i*10 for i in range(l)] # e.g. [-10, 0, 10] when l = 3
			for n in range(l):
				dx = lst[n]['dx']
				e, i = lst[n]['e,i']
				start = lst[n]['start']
				end = lst[n]['end']
				adj = [cof_adj[e], 0, 0.5*abs(dy) -cof_adj[e] + mp_adj[n]]
				segs = get_path_segments(start,	end, 'h', 'h', max_bend, adj)
				path_segs[e][i] = segs

	# check if still overlap
	# if overlap, try to adjust path midpoint until there is no overlap
	print 'checking path/label overlap...'

	for e in edges:
		for i in range(len(path_segs[e])):
			seg = path_segs[e][i]

			p = Path(*seg)
			# get 100 (x + yj) points on the path
			points = [p.point(t) for t in [k/100. for k in range(101)]]
			overlap =[]
			for n in s_nodes.difference({e[1]}):
				# check if any points are inside species labels
				overlap.append(overlapping(points, pos[n], label_size[n]))
			if any(overlap):
				start = complex(*p_nodes[e][i])
				end = complex(*p_nodes[e][i+1])
				if direction[e][i] == 'v' and start.real==end.real or direction[e][i] == 'h' and start.imag == end.imag:
					print 'could not find non-overlapping path for {}'.format(e)
				elif direction[e][i] == direction[e][i+1]:
					print 'adjusting path {} to prevent path/label overlap...'.format(e)
					adjustments = []
					if direction[e][i] == 'h':
						len_seg = 0.5*abs(end.real - start.real)
					else:
						len_seg = 0.5*abs(end.imag - start.imag)

					len_seg = len_seg - cof_adj[e]

					for k in range(1, int(2*len_seg/30)+1):
						adjustments.append(30*k)

					for dn in range(2,6):
						# move middle path segments toward reaction node 
						for nm in range(1, dn):
							adj = nm*(len_seg/dn)
							if adj not in adjustments:
								adjustments.append(adj)
					for dn in range(2,6):
						# move middle path segments toward species node 
						for nm in range(1, dn):
							adj = len_seg + (dn-nm)*(len_seg/dn)
							if adj not in adjustments:
								adjustments.append(adj)

					# adjust, calculate new path and check for overlap
					adjusted_paths = []
					num_overlap = []
					for adj in adjustments:
						print '.',
						sys.stdout.flush()
						new_segs = get_path_segments(start, end, direction[e][i], direction[e][i+1], max_bend, adjust=[0,0,adj])
						p = Path(*new_segs)
						points = [p.point(t) for t in [k/100. for k in range(101)]]
						overlap = []
						for n in s_nodes.difference({e[1]}):
							overlap.append(overlapping(points, pos[n], label_size[n]))
						num_overlap.append(sum(overlap))
						adjusted_paths.append(new_segs)
						if not any(overlap):
							print 'ok'
							path_segs[e][i] = new_segs
							break
					else:
						print 'could not find non-overlapping path'
						path_segs[e][i] = adjusted_paths[num_overlap.index(min(num_overlap))]

	svg_paths = arc_paths
	for e in edges:
		p = []
		for segs in path_segs[e]:
			p+=segs
		svg_paths[e] = Path(*p).d()
	for r in cofactors:
		for s in cofactors[r]:
			svg_paths[(r,s)] = cofactors[r][s]['path']

	return svg_paths, pos
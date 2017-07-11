import cbmpy as cbm
import json
import sys
import re
import argparse
from PIL import ImageFont
from readers import read_json_data, get_cofactors_from_sbml
from svg_assembly import get_svgdata, get_svgdoc

def main(args):

	# get file with layout data
	file_name = ' '.join(args.json_file)
	with open(file_name) as json_data:
		data = json.load(json_data)
	d = read_json_data(data) # layout data
	
	# get font
	font = ImageFont.truetype(args.font_file, 1000)

	# add cofactors
	if args.add_cofactors_from_sbml:
		sbml_file = ' '.join(args.add_cofactors_from_sbml)
		cofactors = get_cofactors_from_sbml(d, sbml_file)
		for r in cofactors:
			for s in cofactors[r]:
				d['edge_type'][(r,s)]=cofactors[r][s]['role']
	else:
		cofactors = None

	# change labels to ids	
	if args.ids_as_label:
		for n in d['label']:
			d['label'][n] = re.sub('_copy_[0-9]+', '', n)
		for r in cofactors:
			for s in cofactors[r]:
				cofactors[r][s]['label'] = s
	
	# get the data to assemble the svg file (editable version)
	svg_data = get_svgdata(
		d = d,
		font = font, 
		font_size = args.font_size, 
		scale = args.scale,
		padding = args.padding,
		padding_labels= args.padding_labels,
		normalize = args.normalize,
		overlap = args.overlap,
		cofactors = cofactors,
		defdir = args.r_direction,
		reverse_cof = args.reverse_cof)
	
	if args.output_json:
		# save svg data in json-format
		with open(args.output_json, 'wb') as f:
			json.dump(svg_data, f)
	else:
		# assemble svg file and save (editable version)
		doc = get_svgdoc(**svg_data)
		doc.save(args.svg_name)
		print 'output svg saved in', args.svg_name

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('json_file', metavar= 'file_name.json', nargs= '+', help = "A json file containg the output from nicholas.")
	parser.add_argument('--svg_name', '-o', default = 'temp.svg', metavar = 'temp.svg', help = "The name/path of the output svg.")
	parser.add_argument('--output_json', '-oj', default = '', metavar = 'svgdata.json', help = "Don't save svg-file, instead save data (info on label coordinates, paths etc.) for creating the svg file in a json file.")
	parser.add_argument('--add_cofactors_from_sbml', '-cof', metavar='model.xml', nargs='+', help = "Add the omitted cofactors from this sbml (3fbc) model to the graph.")
	parser.add_argument('--scale', '-s', type = float, nargs='+', default = [20.0, 20.0], metavar= '20.0', help = "Scale up the graph with this factor. Example: -s 10.0 (10 in both x- and y-direction) Example: -s 20 10 (20 in x-direction, 10 in y-direction")
	parser.add_argument('--padding', type = float, nargs='+', default = [20.0, 20.0], metavar= '20', help = "Extra space (pixels) added to the edges of the svg, e.g. so that all labels are visible in a browser.")
	parser.add_argument('--padding_labels', nargs='+', metavar= '10', help = "Space (pixels) around the text of the labels. Can also accept two terms, for x and y-direction.") 
	parser.add_argument('--font_file', default = 'fonts/Raleway/Raleway-Regular.ttf', metavar = 'C:\Users\User\Documents\Raleway\Raleway-Regular.ttf', help= "The font to be used. Raleway can be downloaded from https://github.com/google/fonts/blob/master/ofl/raleway/Raleway-Regular.ttf")
	parser.add_argument('--font_size', default = 10.0, metavar = '10.0', help = "Font size of the labels.")
	parser.add_argument('--reverse_cof', nargs='+', default = [], metavar = 'R_0001 R_0020 R_0033', help = "List of reaction ids for which the cofactors must be placed in the opposite way to the default (default is substrates top, products bottom).")
	parser.add_argument('--auto_direction', dest = 'r_direction', action='store_false', help = "Determine the direction (horizontal/vertical) of the arrow through reaction nodes automatically. (default is all reactions vertical)")
	parser.add_argument('--ids_as_label', dest = 'ids_as_label', action='store_true', help = "Use metabolite IDs instead of metabolite names as labels.")
	parser.add_argument('--normalize', dest = 'normalize', action='store_true', help = "Translate all coordinates to positive coordinates, so the svg can be viewed in a browser.") 
	parser.add_argument('--overlap', dest='overlap', action='store_true', help = "Allow for overlap of the reaction arrow paths and metabolite labels.")
	parser.set_defaults(ids_as_label = False)
	parser.set_defaults(r_direction = 'vertical')
	parser.set_defaults(normalize = False)
	parser.set_defaults(overlap = False)
	args = parser.parse_args()
	parser.print_help()
	main(args)
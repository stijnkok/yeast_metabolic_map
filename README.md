About
-----
json ---> svg (editable) ---> svg (extra features)
With this module one can take the json-file with the metabolic network layout info and make a pretty svg. See folder 'metabolic_maps' for examples. The svg can, if needed, be edited in Inkscape. The (edited) svg can then be turned into an svg-file with data on hover when viewed in a browser.

json ---> graphml (editable) ---> svg (editable) ---> svg (extra features)
If the graph in the json-file is not to your satisfaction to start with, one can also first turn the layout info from the json-file to a graphml-file, edit the graph, and turn the graphml-file into an svg.

Use in command line:
-------------------
json_to_svg.py 		: Create an svg-file from the json output from Nicholas' module. This file is editable in an svg-editor such as Inkscape. Use json_to_svg.py --help for information on available flags.
layout_final.py 	: Create an svg-file from the output svg-file of json_to_svg.py and/or graph_to_svg.py with annotations on hover that can be used with visualize.py. Use layout_final.py --help for information on available flags.

json_to_graphml.py 	: Create a graph-file (graphml) from the json output from Nicholas' module. This graph can be edited in a graph editor such as Gephi.
graph_to_svg.py 	: Create an svg from a graph file. This file is editable in an svg-editor such as Inkscape.
snap_to_grid.py 	: Change coordinates of nodes in graph file such that they are on a grid.

Functions and classes:
---------------------
readers.py 		: Contains functions for reading the json, xml and graphml files
svg_paths.py		: Contains functions to make svg-paths (used in json_to_svg.py)
svg_assembly.py		: Contains functions to assemble an svg-file (used in json_to_svg.py)

visualize.py 		: Contains classes for mapping FBA results from cbmpy to an svg graphical map of the metabolic network (as created by json_to_svg.py)

Example
-------

Make an svg from "yeast_5_01_model_xml_nucleotide reactions.json"

1) To create the svg-file with metabolic map, use the command:

python modules\json_to_svg.py "json_files\yeast_5_01_model_xml_nucleotide reactions.json" --svg_name editable_svg_files\Y5_nucleotides.svg --scale 20 15 --padding 20 --normalize --add_cofactors_from_sbml models\Y5.xml

A file with the metabolic map called 'Y5_nucleotides.svg' should now be saved in the editable_svg_files directory. This file can be edited in e.g. inkscape.

2) To create a file with annotations, clickable labels and other information on hover that can be used with the visualize.py module, use the command:

python modules\layout_final.py editable_svg_files\Y5_nucleotides.svg models\Y7.xml r_ s_ --svg_name Y5_nucleotides.svg

A file called 'temp.svg' should now be saved in the metabolic_maps folder and opened in a new tab in your browser.
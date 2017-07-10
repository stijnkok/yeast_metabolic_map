#!/usr/bin/python

"""
From the Y7_easy_edit.svg file (which is easier
to edit in inkscape), create a visualisation compatible graphical
metabolic map with clickable reactions and metabolites, 
with tooltips displaying information.
"""

from pysvg.structure import svg, g
from pysvg.text import text
from pysvg.builders import StyleBuilder
from pysvg.shape import circle, rect, path
from pysvg.parser import parse
from pysvg.core import TextContent
from pysvg.linking import a
import argparse
import sys
import os
import re
import cbmpy as cbm
import json
import webbrowser
from numpy import pi, sin, cos
from PIL import ImageFont

def parse_annotations(model):
    cbm.doFBA(model)
    cbm.doFBAMinSum(model)
    # Regex for database ids. Works for http://identifiers.org/*IDENTIFIER* type links
    chebipattern = re.compile('CHEBI:(.+)')
    keggcpattern = re.compile('kegg.compound.(.+)')
    keggrpattern = re.compile('kegg.reaction.(.+)')
    uniprpattern = re.compile('uniprot.(.+)')
    pubmdpattern = re.compile('pubmed.(.+)')

    annotations = {}
    # species annotations
    for s in model.species:
        annotations[s.id] = {}
        # get database references from miriam annotations
        miriam = s.getMIRIAMannotations()
        DBrefs = {}
        if miriam:
            if len(miriam['is']) > 0:
                lnk = miriam['is'][0]
                if lnk.startswith('urn'):
                    lnk = lnk.split(':')
                    lnk = 'http://identifiers.org/{}/{}'.format(lnk[2], lnk[3])
                annotations[s.id]['link'] = lnk
            else:
                annotations[s.id]['link'] = None
            for lnk in miriam['is'] + miriam['isDescribedBy']:
                if lnk.startswith('urn:miriam:'):
                    lnk = lnk.split(':')
                elif lnk.startswith('http://identifiers.org/'):
                    lnk = lnk.split('/')
                if lnk[-2] in DBrefs:
                    DBrefs[lnk[-2]].append(lnk[-1])
                else:
                    DBrefs[lnk[-2]] = [lnk[-1]]
        else:
            annotations[s.id]['link'] = None

        annotations[s.id]['DBrefs'] = DBrefs
        # add formula
        annotations[s.id]['formula'] = s.getChemFormula()
        # add stoichiometry
        stoichiometry = []
        for rid in s.isReagentOf():
            coef = model.getReaction(rid).getReagentWithSpeciesRef(s.id).coefficient
            stoichiometry.append((rid, coef))
        stoichiometry.sort(key=lambda k: abs(model.getReaction(k[0]).value), reverse=True)
        annotations[s.id]['stoichiometry'] = stoichiometry

    #reaction annotations
    for r in model.reactions:
        annotations[r.id] = {}
        # get database references from miriam annotations
        miriam = r.getMIRIAMannotations()
        DBrefs = {}
        if miriam:
            if len(miriam['is']) > 0:
                lnk = miriam['is'][0]
                if lnk.startswith('urn:miriam:'):
                    lnk = lnk.split(':')
                    lnk = 'http://identifiers.org/{}/{}'.format(lnk[2], lnk[3])
                annotations[r.id]['link'] = lnk
            elif len(miriam['isDescribedBy']) > 0:
                lnk = miriam['isDescribedBy'][0]
                if lnk.startswith('urn:miriam:'):
                    lnk = lnk.split(':')
                    lnk = 'http://identifiers.org/{}/{}'.format(lnk[2], lnk[3])
                annotations[r.id]['link'] = lnk
            else:
                annotations[r.id]['link'] = None

            for lnk in miriam['is'] + miriam['isDescribedBy']:
                if lnk.startswith('urn:miriam:'):
                    lnk = lnk.split(':')
                elif lnk.startswith('http://identifiers.org/'):
                    lnk = lnk.split('/')
                if lnk[-2] in DBrefs:
                    DBrefs[lnk[-2]].append(lnk[-1])
                else:
                    DBrefs[lnk[-2]] = [lnk[-1]]
        else:
            annotations[r.id]['link'] = None

        annotations[r.id]['DBrefs'] = DBrefs
        
        # get gene associations
        gpr = model.getGPRforReaction(r.id)
        if gpr:
            annotations[r.id]['GENE_ASSOCIATION'] = gpr.assoc
        else:
            annotations[r.id]['GENE_ASSOCIATION'] = None

    return annotations

def get_arc_paths(cx, cy, r, num):
    """ return circle segments (as pysvg path objects) """
    points = [(cx+r, cy)]
    for i in range(1, num):
        t = i*2*pi/num
        x = cx + r*cos(t)
        y = cy - r*sin(t)
        points.append((x,y))
    points.append((cx+r, cy))
    
    arcs = []
    for i in range(len(points)-1):
        p1 = points[i]
        p2 = points[i+1]
        d = 'M{},{} L{},{} A {},{} 0 0 0 {},{} Z'
        d = d.format(cx, cy, p1[0], p1[1], r, r, p2[0], p2[1])
        arcs.append(d)
    return arcs

def get_label_width(label, font, font_size):
    """" get width of text 'label' (in pixels) in svg file """
    font_size = float(font_size)
    width, _ = font.getsize(label)
    width = width*(font_size/font.size)
    return width

def get_svg_elements(svg_el, notextcontent=True):
    """" parse svg 'easy_edit' file. Return svg elements
    as a list of pysvg objects. """
    elements = []
    i=0
    exit_loop = False
    while not exit_loop: 
        try:
            elements.append(svg_el.getElementAt(i))
            i+=1
        except IndexError:
            exit_loop = True
            pass
        if i>1000000:
            break
    if notextcontent:
        # remove newline chars between elements
        elements = [e for e in elements if len(dir(e))>6] 
    return elements


def get_layout_from_easy_edit(svgdoc, model, r_suffix = 'r_', s_suffix = 's_'):
    # parse svg file
    svg_elements = get_svg_elements(svgdoc)
    
    # dictionaries holding layout information
    rxn_id = {}
    met_id = {}
    rxn_layout = {}
    labels = {}

    # regex for path ids
    pathpattern = re.compile('path_({}.+?)({}.+)'.format(r_suffix, s_suffix))
    # path_rxn_met = 'path_({}.+?)({}.+)'.format(r_suffix, s_suffix)
    # path_rxn_met_copy = 'path_({}.+?)({}.+?)(_.+$)'.format(r_suffix, s_suffix)
    # pathpattern = re.compile('({})|({})'.format(path_rxn_met_copy, path_rxn_met))

    for e in svg_elements:
        e_id = e.get_id()
        if not e_id:
            e_id = ''
        # path elements (a.k.a. reaction arrows)
        match_path = pathpattern.match(e_id)
        if match_path:
            rid = re.sub('_copy_[0-9]+', '', match_path.group(1))
            sid = re.sub('_copy_[0-9]+', '', match_path.group(2))
            e_id = match_path.group(1)
            # if match_path.group(1):
            #     # copy 
            #     rid = match_path.group(2)
            #     sid = match_path.group(3)
            #     e_id = rid + match_path.group(4)
            # else:
            #     rid = match_path.group(6)
            #     sid = match_path.group(7)
            #     e_id = rid

            rxn_id[e_id] = rid

            d = e.get_d() # get svg path ('d' element)
            # marker = re.search('marker-end:url\(#([a-z]*)\)', e.get_style()).group(1) # get marker
            marker = model.getReaction(rid).getReagentWithSpeciesRef(sid).role # get marker from model
            # put layout information in rxn_layout dictionary
            if e_id in rxn_layout:
                if 'paths' in rxn_layout[e_id]:
                    rxn_layout[e_id]['paths'].append({'d':d, 'marker':marker})
                else:
                    rxn_layout[e_id]['paths'] = [{'d':d, 'marker':marker}]
            else:
                rxn_layout[e_id] = {'paths': [{'d':d, 'marker':marker}]}

        # metabolite labels
        elif e_id.startswith(s_suffix):
            sid = re.sub('{}.+$'.format(r_suffix), '', re.sub('_copy_[0-9]+', '', e_id))
            met_id[e_id] = sid
            x = float(e.get_x())
            y = float(e.get_y())
            if 'font-size' in str(e.get_style()):
                font_size = re.search('font-size:([0-9]*)', e.get_style()).group(1)
            else:
                font_size = e.get_font_size()
            txt = e.getElementAt(0).content
            txt_hvr = model.getSpecies(sid).name
            labels[e_id] = {'x':x, 'y':y, 'label_text': txt, 'hover_text': txt_hvr, 'font_size': font_size}
        
        # reaction circles
        elif e_id.startswith(r_suffix):
            x = float(e.get_cx())
            y = float(e.get_cy())
            if e_id in rxn_layout:
                rxn_layout[e_id]['circle'] = {'x':x, 'y':y}
            else:
                rxn_layout[e_id]= {'circle': {'x':x, 'y':y}}
            if not e_id in rxn_id:
                rxn_id[e_id] = re.sub('_copy_[0-9]+', '', e_id)
        
        # flux value placeholders
        elif e_id.startswith('rval_'):
            x = float(e.get_x())
            y = float(e.get_y())
            e_id = e_id[5:]
            if e_id in rxn_layout:
                rxn_layout[e_id]['value'] = {'x':x, 'y':y}
            else:
                rxn_layout[e_id]= {'value': {'x':x, 'y':y}}

    # extract extra layout information from the model
    for e_id in rxn_layout:
        rid = rxn_id[e_id]

        # reaction labels (to display in tooltip)
        labels[e_id] = model.getReaction(rid).name 
        
        # gene associations
        GPR = model.getGPRforReaction(rid)
        if GPR:
            genes = GPR.generefs
            if len(genes)==1:
                rxn_layout[e_id]['genes'] = genes
            elif len(genes)>1:
                # if multiple genes, get circle segments
                arcs = get_arc_paths(rxn_layout[e_id]['circle']['x'], rxn_layout[e_id]['circle']['y'], 10, len(genes))
                rxn_layout[e_id]['genes'] = []
                for i in range(len(genes)):
                    rxn_layout[e_id]['genes'].append({'d':arcs[i], 'id':genes[i]})
            else:
                rxn_layout[e_id]['genes'] = []
        else:
            rxn_layout[e_id]['genes'] = []

    return [rxn_id, met_id, rxn_layout, labels]

### structure of input ###

# font : ImageFont.truetype object

# rxn_id     {node_id: cmod_id,
#             node_id: cmod_id,
#             etc.}

# met_id     {node_id: cmod_id
#             node_id: cmod_id,
#             etc.}

# rxn_layout {node_id : {
#                        'genes' : [{'d' : svg_path,
#                                    'id': cmod_id},
#                                   {'d' : svg_path,
#                                    'id': cmod_id}, etc.],
#                        'paths' : [{'d'     : svg_path,
#                                    'marker': 'substrate' or 'product'},
#                                   {'d'     : svg_path,
#                                    'marker': 'substrate' or 'product'}, etc.],
#                        'circle': {'x' : x-coordinate,
#                                   'y' : y-coordinate (svg!)},
#                        'value' : {'x' : x-coordinate,
#                                   'y' : y-coordinate (svg!)}
#                        }, 
#             etc.}

# labels {node_id: {'x'         : x-coordinate,
#                   'y'         : y-coordinate (svg!),
#                   'label_text': label,
#                   'hover_text': label on hover,
#                   'font_size' : font size},
#         etc.}

# annotations {cmod_id (reaction) : {'link': reference link,
#                                    'DBrefs': {'KEGG':[kegg_id1, kegg_id2], 'PUBMED':[pubmed_id1, pubmed_id2], etc.},
#                                    'GENE_ASSOSCIATION' : gene association},
#              cmod_id (species) : {'link': reference link,
#                                    'DBrefs': {'KEGG':[kegg_id1, kegg_id2], 'PUBMED':[pubmed_id1, pubmed_id2], etc.},
#                                    'formula' : chemical formula,
#                                    'stoichiometry': [(rid, coef), (rid, coef), etc.]},
#              etc. }

def wrap_svg_metabolic_map(svg, css, height, width, title):
    #<?xml-stylesheet type="text/css" href="svg.css" ?>
    return """<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg" height="{}" width="{}" version="1.1" xmlns:xlink="http://www.w3.org/1999/xlink">

  <title>{}</title>
  <text x="0" y="10" font-family = "Raleway" font-size="7" style="fill:black">2017, Amsterdam, Netherlands</text>
  <text x="0" y="18" font-family = "Raleway" font-size="7" style="fill:black">(c) Stijn Kok, stijn.kok@student.uva.nl</text>

  <style>
{}
  </style>

  <defs>
    <marker id="product" markerWidth="6" markerHeight="6" refX="0.75" refY="2.5" orient="auto" markerUnits="strokeWidth">
      <path
         d="M 1,1 Q 3,2.5 1,4 L 4,2.5 Z"
         style="fill:#000000;stroke:#000000;stroke-width:0.5;stroke-linecap:round;stroke-linejoin:round" />
    </marker>
    <marker id="substrate" markerWidth="6" markerHeight="6" refX="0.75" refY="2.5" orient="auto" markerUnits="strokeWidth">
      <circle
         cx="2.5999999"
         cy="2.5"
         r="1.1"
         style="fill:#000000" />
    </marker>
  </defs>

{}
</svg>""".format(height, width, title, css, svg)

  # <rect style="fill:none;stroke:#aaaaaa;stroke-width:50" rx="50" height="1015" width="2430" y="1600" x="645" id="rect5104"  />
  # <rect style="fill:#dddddd;fill-opacity:1" rx="10" height="37.187725" width="106.97604" y="2476.406" x="601.78857" id="rect5277-49"  />
  # <rect style="fill:#dddddd;fill-opacity:1" rx="10" height="70.745087" width="109.4573" y="2119.6267" x="601.78857" id="rect5277-4"  />
  # <rect style="fill:#dddddd;fill-opacity:1" rx="10" height="70.745087" width="109.4573" y="2009.6272" x="601.78857" id="rect5277-6"  />
  # <rect style="fill:#dddddd;fill-opacity:1" rx="10" height="70.745087" width="109.4573" y="2234.6272" x="601.78857" id="rect5277"  />
  # <text font-size="40" style="font-size:40px;font-family:Raleway;fill:#ffffff" id="text5106" y="2629.3992" x="1196.4995">Mitochondrion</text>

############# MAIN ###############
def main(args):

    svgdoc = parse(args.svg_easy_edit_file)
    model = cbm.CBRead.readSBML3FBC(args.SBML_file)
    font = ImageFont.truetype(args.font_file, 1000)
    # get layout infromation from svg file and model
    rxn_id, met_id, rxn_layout, labels = get_layout_from_easy_edit(svgdoc, model, args.r_suffix, args.s_suffix)

    # 'annotations' is a dictionary, 
    # with keys: the reaction/metabolite ids
    # with (in case of reaction) values: dictionary with keys 'link', 'DBrefs', 'GENE_ASSOCIATION'
    # with (in case of species) values: dictionary with keys 'link', 'DBrefs', 'formula', 'stoichiometry'
    if args.annotations:
        with open(args.annotations) as f:
            annotations = json.load(f)
    else:
        annotations = parse_annotations(model)

    
    css = """
    text {{font-family:"{}"}}
    .annotations       circle {{opacity: 0}}
    .annotations:hover circle {{opacity: 1}}
    .annotations a       {{opacity: 0}}
    .annotations:hover a {{opacity: 1}}
    .annotations       .tooltip {{visibility: hidden}}
    .annotations:hover .tooltip {{visibility: visible}}
""".format(font.getname()[0])
    
    if font.getname()[0] == 'Raleway':
        css += """
    @font-face {
      font-family: 'Raleway';
      src: url(data:application/font-woff2;charset=utf-8;base64,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_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) format('woff2'),
           url(data:application/font-woff;charset=utf-8;base64,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) format('woff');
      font-weight: normal;
      font-style: normal
      }
"""

    rxn_style = """
    #{r} {{stroke:#cccccc; stroke-width:1.0; stroke-dasharray:1.5}}
    #{r}.substrate {{marker-end:url(#substrate)}}
    #{r}.product   {{marker-end:url(#product)}}
    #{r}.fluxvalue_tooltip {{fill-opacity:0}}
    #{r}.fluxvalue {{fill-opacity:0}}
    #{r}.FVAspan {{fill-opacity:0}}
    #{r}.FVAmin {{fill-opacity:0}}
    #{r}.FVAmax {{fill-opacity:0}}
    #{r}.reversible {rev}
    #{r}.irreversible {irr}
    #{r}.inactive {ko}
"""

    met_style = """
    #{s} {{fill:black}}"""

    gen_style = """
    #{g} {{fill:none; stroke:#000000; stroke-width:1.0; stroke-opacity:0}}"""


    # start concatenating new svg file
    svg = ''

    # add reaction style (set reaction reversibility based on reaction bounds in model)
    for r in model.reactions:
        if r.getLowerBound() == 0 or r.getUpperBound() == 0:
            # irreversible
            rev = '{stroke-opacity:0}'
            irr = '{stroke-opacity:1}'
            ko = '{stroke-opacity:0}'
        elif r.getLowerBound() == 0 and r.getUpperBound() == 0:
            # knocked out/ inactive
            rev = '{stroke-opacity:0}'
            irr = '{stroke-opacity:0}'
            ko = '{stroke-opacity:1}'
        else:
            # reversible
            rev = '{stroke-opacity:1}'
            irr = '{stroke-opacity:0}'
            ko = '{stroke-opacity:0}'
        css += rxn_style.format(r=r.id, rev=rev, irr=irr, ko=ko)


    for r, rid in rxn_id.items():

        # get coordinates
        x = rxn_layout[r]['circle']['x']
        y = rxn_layout[r]['circle']['y']

        # add element for visualizing gene expression
        G_gn = g()
        G_gn.set_class('gene_expression')

        if len(rxn_layout[r]['genes']) == 1:
            # if single gene, make a circle
            C = circle(x,y, 10)
            C.set_id(rxn_layout[r]['genes'][0])
            G_gn.addElement(C)
            # add gene style
            css += gen_style.format(g=rxn_layout[r]['genes'][0]) 
        else:
            # if multiple genes, make path segments
            # (or do nothing if no genes)
            for layout in rxn_layout[r]['genes']:
                P = path(layout['d'])
                P.set_id(layout['id'])
                G_gn.addElement(P)
                # add gene style
                css += gen_style.format(g=layout['id'])
        svg += G_gn.getXML()

        # add reaction
        G_rn = g()
        G_rn.set_id(rid)
        G_rn.set_style('fill:none; stroke-linecap:round; stroke-linejoin:round')

        # add paths (reaction arrows) to reaction
        if 'paths' in rxn_layout[r]:
            for layout in rxn_layout[r]['paths']:
                P = path(layout['d'])
                P.set_id(rid)
                P.set_class(layout['marker'])
                G_rn.addElement(P)

        # add circle (displaying reversibility) to reaction
        C = circle(x, y, 5)
        C.set_id(rid) 
        C.set_style('fill:#ffffff; stroke-dasharray:none')
        G_rn.addElement(C)
        
        P_rev = path('m {},{} h 2 m -2,2 h 2'.format(x-1, y-1))
        P_rev.set_id(rid)
        P_rev.set_class('reversible')
        P_rev.set_style('stroke:#000000; stroke-width:1.0; stroke-dasharray:none')
        G_rn.addElement(P_rev)

        P_irr = path('m {},{} 2,1 -2,1'.format(x-1, y-1))
        P_irr.set_id(rid)
        P_irr.set_class('irreversible')
        P_irr.set_style('stroke:#000000; stroke-width:1.0; stroke-dasharray:none')
        G_rn.addElement(P_irr)

        P_ko = path('m {},{} 3,3 m -3,0 3,-3'.format(x-1.5, y-1.5))
        P_ko.set_id(rid)
        P_ko.set_class('inactive')
        P_ko.set_style('stroke:#ff0000; stroke-width:1.0; stroke-dasharray:none')
        G_rn.addElement(P_ko)

        # add flux value placeholder
        x = rxn_layout[r]['value']['x']
        y = rxn_layout[r]['value']['y']

        T = text('ReactionValue:abs2:'+rid, x, y) # ':abs2:' means display absolute value in 2 significant digits
        T.set_id(rid)
        T.set_class('fluxvalue')
        T.set_style('font-size:8px; fill:black; stroke-opacity:0')
        G_rn.addElement(T)

        # add span value placeholder
        T = text('ReactionSpan:2:'+rid, x, y) # :2: means display value in 2 significant digits
        T.set_id(rid)
        T.set_class('FVAspan')
        T.set_style('font-size:8px; fill:black; stroke-opacity:0')
        G_rn.addElement(T)

        svg += G_rn.getXML()
        
    css+='\n'

    for s, sid in met_id.items():
        # add metabolite labels
        x = labels[s]['x']
        y = labels[s]['y']
        T = text(labels[s]['label_text'], x, y)
        T.set_id(sid)
        T.set_style('font-size:{}px;'.format(labels[s]['font_size']))

        svg += T.getXML()
        css += met_style.format(s=sid) # add metabolite style

    # add reaction info on hover

    for r, rid in rxn_id.items():

        G_nfo = g()
        G_nfo.set_id(rid)
        G_nfo.set_class('annotations')

        x = rxn_layout[r]['circle']['x']
        y = rxn_layout[r]['circle']['y']

        C = circle(x, y, 5)
        C.set_style('fill:#000000')

        href = annotations[rid]['link']
        if href:
            A = a()
            A.set_target('_blank')
            A.set_xlink_href(href)
            A.addElement(C)
            G_nfo.addElement(A)
        else:
            G_nfo.addElement(C)

        DBrefs = [] # database references
        for DBname, reflist in annotations[rid]['DBrefs'].items():
            DBrefs += [DBname+':'+ref for ref in reflist]
        DBrefs = ', '.join(DBrefs)

        gene_assoc = annotations[rid]['GENE_ASSOCIATION']
        if not gene_assoc:
            gene_assoc = 'NA'

        G_ttp = g()
        G_ttp.set_class('tooltip')
        w1 = get_label_width(labels[r], font, 10)
        w2 = get_label_width('Genes: '+gene_assoc, font, 10)
        w3 = get_label_width(DBrefs, font, 10)
        w4 = get_label_width('ID: '+rid, font, 10)
        ttp_w = max([100, w1, w2, w3, w4]) +10

        shadow = rect(x+14, y-51, width= ttp_w, height= 110, rx=5)
        shadow.set_style('fill:#000000; stroke:#000000; stroke-width:1.0; opacity:0.1; stroke-dasharray:none')
        G_ttp.addElement(shadow)
        box = rect(x+10, y-55, width= ttp_w, height= 110, rx=5)
        box.set_id(rid)
        box.set_style('fill:#ffffff; stroke-width:1.0; stroke-dasharray:none')
        G_ttp.addElement(box)

        T = text(labels[r], x+15, y-42.5)
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text('ID: '+rid, x+15, y-27.5)
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text(DBrefs, x+15, y-12.5)
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text('Genes: '+gene_assoc, x+15, y+2.5)
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text('ReactionValue:6:'+rid, x+15, y+17.5)
        T.set_id(rid)
        T.set_class('fluxvalue_tooltip')
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text('min: ReactionMinValue:6:'+rid, x+15, y+32.5)
        T.set_id(rid)
        T.set_class('FVAmin')
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        T = text('max: ReactionMaxValue:6:'+rid, x+15, y+47.5)
        T.set_id(rid)
        T.set_class('FVAmax')
        T.set_style('fill:black; font-size:10px; stroke-opacity:0')
        G_ttp.addElement(T)

        G_nfo.addElement(G_ttp)

        svg += G_nfo.getXML()


    img_str = '<image class="tooltip" x="{}" y="{}" width="100" height="60" preserveAspectRatio="xMinYMid meet" xlink:href="http://www.genome.jp/Fig/compound/{}.gif"></image>'
    for s, sid in met_id.items():
        x = labels[s]['x']
        y = labels[s]['y']

        G_nfo = g()
        G_nfo.set_class('annotations')
        
        DBrefs = [] # database references
        for DBname, reflist in annotations[sid]['DBrefs'].items():
            DBrefs += [DBname+':'+ref for ref in reflist]
        DBrefs = ', '.join(DBrefs)

        chemform = annotations[sid]['formula']
        if not chemform:
            chemform='NA'
        
        w1 = get_label_width('Chemical formula: '+chemform, font, 10)
        w2 = get_label_width(DBrefs, font, 10)
        w3 = get_label_width('ID: '+sid, font, 10)
        w_rids=[]
        w_coef=[]
        stoichiometry = annotations[sid]['stoichiometry']
        for rid, coef in stoichiometry:
            w_rids.append(get_label_width(rid, font, 10))
            w_coef.append(get_label_width(str(coef), font, 10))
        w_col1 = max(w_rids)
        w_col2 = max(w_coef)

        w4 = w_col1 + w_col2 + 80
        
        ttp_w = max([110, w1, w2, w3, w4]) + 10
        ttp_h = 130+min(26, len(stoichiometry))*15
        
        # table with reactions
        G_ttp = g()
        G_ttp.set_class('tooltip')

        R = rect(x-2.5, y-12.5, width = get_label_width(labels[s]['hover_text'], font, 10)+5, height = 15)
        R.set_style('fill:#ffffff; stroke:none')
        G_ttp.addElement(R) 
        
        shadow = rect(x+4, y+9, width= ttp_w, height= ttp_h , rx=5)
        shadow.set_style('fill:#000000; stroke:#000000; opacity:0.1; stroke-dasharray:none')
        G_ttp.addElement(shadow)
        box = rect(x, y+5, width= ttp_w, height= ttp_h, rx=5)
        box.set_id(rid)
        box.set_style('fill:#ffffff; stroke:#000000; stroke-width:1.0; stroke-dasharray:none')
        G_ttp.addElement(box)

        if 'kegg.compound' in annotations[sid]['DBrefs']:
            kegg_ref = annotations[sid]['DBrefs']['kegg.compound']
            if kegg_ref:
                kegg_ref=kegg_ref[0]
                img = TextContent(img_str.format(x+5, y+15, kegg_ref))
            else:
                img=TextContent(' ')
        else:
            img=TextContent(' ')

        G_ttp.addElement(img)

        x = x+5
        
        T = text('Chemical formula: '+chemform, x, y + 97.5)
        T.set_style('fill:black; font-size:10px')
        G_ttp.addElement(T)

        T = text('ID: '+sid, x, y + 112.5)
        T.set_style('fill:black; font-size:10px')
        G_ttp.addElement(T)

        T = text(DBrefs, x, y+127.5)
        T.set_style('fill:black; font-size:10px')
        G_ttp.addElement(T)

        for i in range(len(stoichiometry)):
            rid = stoichiometry[i][0]
            coef= stoichiometry[i][1]
            x0 = x
            x1 = x + w_col1 + 5 + w_col2 - w_coef[i]
            x2 = x + w_col1 + w_col2 + 15
            
            T = text(rid, x0, y+142.5 + i*15)
            T.set_style('fill:black; font-size:10px')
            G_ttp.addElement(T)

            T = text(str(coef)+' x ', x1, y+142.5 + i*15)
            T.set_style('fill:black; font-size:10px')
            G_ttp.addElement(T)
            
            T = text('ReactionValue:6:'+rid, x2, y+142.5 + i*15)
            T.set_id(rid)
            T.set_class('fluxvalue_tooltip')
            T.set_style('fill:black; font-size:10px; stroke-opacity:0')
            G_ttp.addElement(T)

            if i == 24:
                T = text('...and {} more'.format(len(stoichiometry)-25), x, y+157.5 + i*15)
                T.set_style('fill:black; font-size:10px')
                G_ttp.addElement(T)
                break
        
        G_nfo.addElement(G_ttp)

        A = a()
        href = annotations[sid]['link']
        if href:
            A.set_target('_blank')
            A.set_xlink_href(href)
        else:
            A.set_xlink_href('#')
        
        T = text(labels[s]['hover_text'], x-5, y)
        T.set_style('fill:#999999; font-size:10px')
        A.addElement(T)
        
        G_nfo.addElement(A)

        svg += G_nfo.getXML()
    
    if args.height:
        height = args.height
    else:
        height = svgdoc.get_height()
    if args.width:
        width = args.width
    else:
        width = svgdoc.get_width()
    if args.title:
        title = args.title
    else:
        title = 'Graphical map of ' + args.SBML_file

    svg = wrap_svg_metabolic_map(svg, css, height, width, title)

    # save file
    if not os.path.exists(os.path.join(os.getcwd(), args.output_dir)):
        os.makedirs(args.output_dir)
    with open(os.path.join(os.getcwd(), args.output_dir, args.svg_name), 'wb') as f:
        f.write(svg)

    if args.open_browser:
        webbrowser.open_new_tab(os.path.join(os.getcwd(), args.output_dir, args.svg_name))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('svg_easy_edit_file', metavar= 'file_name.svg', help ="Editable svg-file")
    parser.add_argument('SBML_file', metavar= 'model.xml', help = "SBML model with annotations")
    parser.add_argument('r_suffix', help = "Reaction suffix")
    parser.add_argument('s_suffix', help = "Species suffix")
    parser.add_argument('--title', help = "Title of output svg")
    parser.add_argument('--height', help = "Height (pixels) of the output svg")
    parser.add_argument('--width', help = "Width (pixels) of the output svg")
    parser.add_argument('--font_file', default = 'C:\Users\User\Documents\Raleway\Raleway-Regular.ttf')
    parser.add_argument('--annotations', metavar = 'annotations.json', default= '')
    parser.add_argument('--svg_name', '-o', default = 'temp.svg')
    parser.add_argument('--output_dir', default = 'metabolic_maps')
    parser.add_argument('--open_browser', type= bool, default = True)
    parser.set_defaults(height = False)
    parser.set_defaults(width = False)
    args = parser.parse_args()
    main(args)
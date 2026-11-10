"""
*********
JSON data
*********
Generate and parse JSON serializable data for blosm_networkx graphs.

These formats are suitable for use with the d3.js examples https://d3js.org/

The three formats that you can generate with blosm_networkx are:

 - node-link like in the d3.js example https://bl.ocks.org/mbostock/4062045
 - tree like in the d3.js example https://bl.ocks.org/mbostock/4063550
 - adjacency like in the d3.js example https://bost.ocks.org/mike/miserables/
"""
from .node_link import *
from .adjacency import *
from .tree import *
from .cytoscape import *

#!/usr/bin/env python3
import graph_tool.all as gt
from collections import defaultdict
import sys

def check_arguments():
    if len(sys.argv) != 2:
        sys.exit("\nERROR: Must provide one argument (links file).\n")
    return sys.argv[1]

def read_graph(filename):
    fh        = open(filename, "r")
    graph     = gt.Graph(directed=False)
    vprop     = graph.new_vertex_property("string")
    scores    = defaultdict(lambda: defaultdict(float))
    node_dict = dict()
    node_idx  = 0
    for line in fh:
        line = line.strip()
        s, score, t = line.split(",")
        scores[s][t] = float(score)

        # Add vertices
        if s and s not in node_dict:
            source         = graph.add_vertex()
            node_dict[s]   = source
            vprop[source]  = s
        if t and t not in node_dict:
            target         = graph.add_vertex()
            node_dict[t]   = target
            vprop[target]  = t

        # Add edge
        if s and t:
            graph.add_edge(node_dict[s], node_dict[t])

    return graph, vprop, scores

filename = check_arguments()
graph, vprop, scores = read_graph(filename)
components, histogram = gt.label_components(graph)

for edge in graph.edges():
    v1, v2 = edge.source(), edge.target()
    print("{},{},{},{},{}".format(vprop[v1], scores[vprop[v1]][vprop[v2]],  vprop[v2], components[v1], components[v2]))
        
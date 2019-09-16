import sys
from collections import defaultdict

def check_arguments():
    if len(sys.argv) != 3:
        sys.exit("\nERROR: Must provide at least two arguments (Experiment:Filename Experiment:Filename).\n")
    return sys.argv[1:]

def read_graph(experiment, link_file, graph):
    with open(link_file, "r") as fh:
        next(fh) # skip header
        for line in fh:
            line = line.strip()
            d1, score, d2, h1, h2, g1, gn1, g2, gn2, reactome = line.split("\t")
            d1 = d1.replace("_0", "_0_1")
            d2 = d2.replace("_0", "_0_1")
            if d1 in graph and d2 in graph[d1]:
                graph[d1][d2]["BOTH"] = True
            else:
                graph[d1][d2]["BOTH"] = False
            
            graph[d1][d2][experiment] = [ float(score), reactome ]

def print_graphs(graph):
    for d1 in graph:
        for d2 in graph[d1]:
            for experiment in graph[d1][d2]:
                if experiment == "BOTH":
                    continue
                print("{}\t{}\t{}\t{}\t{}\t{}".format(experiment, d1, d2, round(graph[d1][d2][experiment][0], 3), graph[d1][d2][experiment][1], graph[d1][d2]["BOTH"]))


def main():
    sys.stderr.write("# Checking arguments.\n")
    link_files = check_arguments()
    

    sys.stderr.write("# Reading graphs.\n")
    graph = defaultdict(lambda: defaultdict(dict))
    for link in link_files:
        experiment, link_file = link.split(":")
        sys.stderr.write("## Reading graph for {}...\n".format(experiment))
        read_graph(experiment, link_file, graph)

    print_graphs(graph)
    sys.stderr.write("# All done.\n")


if __name__ == "__main__":
    main()

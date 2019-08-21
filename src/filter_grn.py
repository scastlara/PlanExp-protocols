#!/usr/bin/env python3

import sys
from enum import Enum
from collections import defaultdict

class Params(Enum):
    SCORE = 0.005
    MAX_TARGETS = 50
    MAX_REGULATORS = 50
        
def check_arguments():
    if len(sys.argv) != 2:
        sys.exit("\nERROR: Must provide one argument (links file).\n")
    return sys.argv[1]

def read_graph(links_file):
    graph = defaultdict(lambda: defaultdict(float))
    num_of_links = 0
    num_of_regulators = set()
    num_of_targets = set()
    with open(links_file, "r") as fh:
        for line in fh:
            line = line.strip()
            regulator, score, target = line.split(",")
            num_of_links += 1
            num_of_regulators.add(regulator)
            num_of_targets.add(target)         
            graph[regulator][target] = float(score)
    print_graph_stats(num_of_links, num_of_regulators, num_of_targets)
    return graph

def apply_filters(graph):
    best_targets = compute_best_targets(graph)
    best_regulators = compute_best_regulators(graph)
    filtered_number = 0
    for regulator, interactions in graph.items():
        for target, score in interactions.items():
            if not is_score_high_enough(score):
                continue
            if not is_target_among_best(regulator, target, best_targets):
                continue
            if not is_regulator_among_best(regulator, target, best_regulators):
                continue
            filtered_number += 1
            print("{},{},{}".format(regulator, score, target))
    sys.stderr.write("# Filtered interactions: {}\n".format(filtered_number))

def is_score_high_enough(score):
    return score >= Params.SCORE.value

def is_target_among_best(regulator, target, best_targets):
    return target in best_targets[regulator]

def is_regulator_among_best(regulator, target, best_regulators):
    return regulator in best_regulators[target]

def compute_best_targets(graph):
    best_targets = defaultdict(set)
    for regulator, interactions in graph.items():
        best_targets_for_regulator = sorted(interactions.keys(), key=lambda x: graph[regulator][x], reverse=True)
        best_targets[regulator] = set(best_targets_for_regulator[:Params.MAX_TARGETS.value])
    return best_targets

def compute_best_regulators(graph):
    best_regulators = defaultdict(set)
    regulators_for_target = defaultdict(set)
    for regulator, interactions in graph.items():
        for target in interactions.keys():
            regulators_for_target[target].add(regulator)
    
    for target, regulators in regulators_for_target.items():
        best_regulators[target] = set(sorted(regulators, key=lambda x: graph[x][target], reverse=True)[:Params.MAX_REGULATORS.value])
    return best_regulators

def print_graph_stats(num_of_links, num_of_regulators, num_of_targets):
    sys.stderr.write("\n# Graph stats:\n")
    sys.stderr.write("## Number of Links: {}\n".format(num_of_links))
    sys.stderr.write("## Number of Regulators: {}\n".format(len(num_of_regulators)))
    sys.stderr.write("## Number of Targets: {}\n\n".format(len(num_of_targets)))



def main():
    sys.stderr.write("# Checking arguments.\n")
    links_file = check_arguments()

    sys.stderr.write("# Reading graph.\n")
    graph = read_graph(links_file)

    sys.stderr.write("# Applying filters.\n")
    apply_filters(graph)
    sys.stderr.write("# All done.\n")


if __name__ == "__main__":
    main()
#!/usr/bin/env python3
import sys
from collections import defaultdict


def check_arguments():
    if len(sys.argv) != 3:
        sys.exit("\nERROR: Must provide two arguments (reactome file and links file).\n")
    return sys.argv[1], sys.argv[2]

def read_reactome(filename):
    interaction_to_reactome = defaultdict(lambda: defaultdict(list))
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            name, identifier, _, *members = line.split("\t")
            for member1 in members:
                if member1 == "NA":
                    continue
                for member2 in members:
                    if member2 == "NA":
                        continue
                    interaction_to_reactome[member1][member2].append((name, identifier))
    return interaction_to_reactome

def reactome_for_link(reactome_membership, a, b):
    if a in reactome_membership and b in reactome_membership[a]:
        return reactome_membership[a][b]
    else:
        return []


def format_reactome_pathways(reactome_pathways):
    reactome_pathways = set(reactome_pathways)
    if reactome_pathways:
        return [ "{}-{}".format(p[1], p[0]) for p in sorted(reactome_pathways) ]
    else:
        return ["NA"]


def print_links_with_reactome(links_filename, reactome_membership):

    with open(links_filename, "r") as fh:
         # skip header
        header = fh.readline()
        header = header.strip()
        print("{},ReactomePathways".format(header))
        for line in fh:
            reactome_pathways = []
            line = line.strip()
            columns = line.split(",")
            homolog1, homolog2 = columns[3], columns[4]
            name1, name2 = columns[6], columns[8]
            reactome_pathways.extend(reactome_for_link(reactome_membership, homolog1, homolog2))
            reactome_pathways.extend(reactome_for_link(reactome_membership, name1, name2))
            reactome_pathways = format_reactome_pathways(reactome_pathways)
            print("{},{}".format(",".join(columns), ";".join(reactome_pathways)))





def main():
    reactome_filename, links_filename = check_arguments()
    reactome_membership = read_reactome(reactome_filename)
    print_links_with_reactome(links_filename, reactome_membership)


if __name__ == "__main__":
    main()

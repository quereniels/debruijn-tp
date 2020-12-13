#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    i = -1
    for row in open(fastq_file, "r"):
        if i%4 == 0: 
            row = row.rstrip("\n")
            yield row
        i += 1        
                
def cut_kmer(read, kmer_size):
    for i in range(len(read)-kmer_size+1):
        yield read[i:kmer_size+i]


def build_kmer_dict(fastq_file, kmer_size):
    with open(fastq_file) as f:
        for r, l in enumerate(f):
            pass
    r += 1
    genS = read_fastq(fastq_file)
    dct = {}
    for i in range(r//4):
        S = next(genS)
        genK = cut_kmer(S , kmer_size)
        for j in range(len(S)-kmer_size+1):
            K = next(genK)
            dct[K] = dct.get(K,0) + 1
    
    return dct
            
            

def build_graph(kmer_dict):
    G = nx.Graph()
    H = nx.DiGraph(G)
    for k, v in kmer_dict.items():
        x = k[0:(len(k)-1)]
        y = k[1:len(k)]

        H.add_edge(x,y, weight = v)
    return H

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node is True:
            graph.remove_node(path[0])
        if delete_sink_node is True:
            graph.remove_node(path[-1])
        graph.remove_nodes_from(path[1:-1])
    return graph
    
def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):

    best = 0
    for i in range(1, len(path_list)):
        weight1 = weight_avg_list[best]
        weight2 = weight_avg_list[i]
        if weight1 < weight2:
            best = i
        elif weight1 == weight2:
            length1 = path_length[best]
            length2 = path_length[i]
            if length1 < length2:
                best = i
            elif length1 == length2:
                best = random.choice([best, i])


    path_list.remove(path_list[best])
    graph = remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph

    pass

def path_average_weight(graph, path):
    
    weight_list = []
    for i in path:
        if graph.get_edge_data(i,i+1) != None:
            weight_list.append(graph.get_edge_data(i,i+1).get("weight"))        
    return statistics.mean(weight_list)
    
def solve_bubble(graph, ancestor_node, descendant_node):
    paths = nx.all_simple_paths(graph,ancestor_node,descendant_node)
    path_list = []
    path_length = []
    avg_weight_list =[]
    for path in paths:
        path_list.append(path)
        path_length.append(len(path))
        avg_weight_list.append(path_average_weight(graph, path))
    solved = select_best_path(graph, path_list, path_length, avg_weight_list)
    return solved
    

def simplify_bubbles(graph):
    pass
            
        
    

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    return [n for n, d in graph.in_degree() if d == 0]

def get_sink_nodes(graph):
    return [n for n, d in graph.out_degree() if d == 0]

def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for nodeS in starting_nodes:
        for nodeE in ending_nodes:
            for paths in nx.all_simple_paths(graph, nodeS, nodeE):
                path = paths[0]    
                path2 = paths[1:] 
                for node in path2:
                    node2 = node[1:]
                    path = path + node2 
                contigs.append( (path, len(path)) )
    return contigs


def save_contigs(contigs_list, output_file):
    with open(output_file,'w') as f:
        c=0
        for seq,sz in contigs_list:
            print(seq,sz)
            c+=1
            f.write(f">contig_{c} len={sz}\n")
            f.write(f"{fill(seq)}\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    # Create graph
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)

    # Solve bubbles
    #graph = simplify_bubbles(graph)

    # Get and save contigs
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    save_contigs(contigs_list, args.output_file)

if __name__ == '__main__':
    main()

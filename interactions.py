

from objects import *
from util import *
import networkx as nx
from itertools import combinations
from collections import Counter
import math
import numpy as np


drugs,proteins = parse(
    '/home/stanislaw/Downloads/git_drugbank/download/full_database.xml')


def as_graph():
    G = nx.MultiDiGraph()
    for p in proteins:
        G.add_node(p,typ='target')
    for d in drugs.values():
        G.add_node(d.dbid,typ='drug')
        for e in d.edges:
            target = e[1] if isinstance(e[1],str) else e[1].upid
            G.add_edge(d.dbid,target,itype=e[0])
    return G




def node_degree(G):
    ddegrees = []
    pdegrees = {}
    oneedge = []
    fewedges = []
    isolated = 0
    

    print('Number of nodes = %d' % G.number_of_nodes())
    print('Number of edges = %d' % G.number_of_edges())
    degrees = [nx.degree(G,node) for node in G]
    print('Largest hub = '+ max([(node,nx.degree(G,node)) for node in G],key=lambda d: d[1])[0])
    drugand = sum(degrees) / len(degrees)
    print('Average node degree = %f' % drugand)
    degrees.sort(reverse=True)
    hubsvalues = degrees[:math.ceil(len(degrees) / 20)]
    print('Minimal degree of hub = %d' % hubsvalues[-1])
    print('Maximum degree of hub = %d' % hubsvalues[0])

    hubs = []
    for node in G:
        if nx.degree(G,node) >= hubsvalues[-1]:
            hubs.append(node)
    
    print('Number of hubs = %d' % len(hubs))

def stats(G)
    # Centrality metrics
    bc = nx.betweenness_centrality(G)
    print('Average betweeness centrality: %f' % np.mean(list(bc.values())))

    # Connected component
    maks = len(max(list(nx.connected_components(G)),key=lambda g: len(g)))
    print('Giant component size: %d' % maks)

    # Average path length
    for g in nx.connected_component_subgraphs(G):
        if nx.number_of_nodes(g) == maks:
            gc = g
            break
    spath = nx.average_shortest_path_length(gc)
    print('Average path length: %f' % spath)

    # Clustering coefficient
    ac = nx.average_clustering(G)
    print('Clustering coefficient: %f' % ac)

    # Network density
    nd = nx.density(G)
    print('Network density: %f' % nd)

    # Degree distribution, network entropy
    counter = Counter(list(dict(nx.degree(G)).values()))
    distribution = {}
    for k in counter.keys():
        distribution[k] = counter[k]/nx.number_of_nodes(G)
    entropy = 0
    for p in distribution.values():
        entropy += p*math.log(p)
    entropy *= -1
    print('Network entropy: %f' % entropy)


drugs,proteins = parse(
    '/home/stanislaw/Downloads/git_drugbank/download/full_database.xml')
G = as_graph()

H = nx.Graph()
G_ud = nx.Graph(G)
pairs = combinations(list(drugs.keys()),2)
for u,v in pairs:
    comm_neigh = list(nx.common_neighbors(G_ud,u,v))
    if len(comm_neigh)>0:
        H.add_weighted_edges_from([(u,v,len(comm_neigh))])


I =nx.Graph()
pairs = combinations(list(proteins.keys()),2)
for u,v in pairs:
    comm_neigh = list(nx.common_neighbors(G_ud,u,v))
    if len(comm_neigh)>0:
        I.add_weighted_edges_from([(u,v,len(comm_neigh))])

node_degree(H)
stats(H)


node_degree(I)
stats(I)


import networkx as nx
import util
import numpy as np
from collections import Counter
import math
import matplotlib.pylab as plt
import objects


def node_degree(drugs):
    ddegrees = []
    pdegrees = {}
    oneedge = []
    fewedges = []
    isolated = 0
    for drug in drugs.values():
        if len(drug.edges) == 0:
            isolated += 1
        for protein in [d[1] for d in drug.edges]:
            if not isinstance(protein, str):
                protein = protein.name
            if protein not in pdegrees:
                pdegrees[protein] = 1
            else:
                pdegrees[protein] += 1
        if len(drug.edges) == 1:
            oneedge.append(drug)
        if len(drug.edges) < 4:
            fewedges.append(drug)
        ddegrees.append(len(drug.edges))

    print('Number of nodes = %d + %d' % (len(drugs), len(pdegrees)))
    print('Number of isolated nodes = %d' % isolated)
    degrees = ddegrees + [v for v in pdegrees.values()]
    drugand = sum(ddegrees) / len(ddegrees)
    print('Average node degree of drugs = %f' % drugand)
    proteinand = sum(pdegrees.values()) / len(pdegrees)
    print('Average node degree of proteins = %f' % proteinand)
    alland = sum(degrees) / len(degrees)
    print('Average node degree of all network = %f' % alland)

    degrees.sort(reverse=True)
    hubsvalues = degrees[:math.ceil(len(degrees) / 20)]
    print('Minimal degree of hub = %d' % hubsvalues[-1])
    print('Maximum degree of hub = %d' % hubsvalues[0])

    drughubs = []
    for drug in drugs.values():
        if len(drug.edges) >= hubsvalues[-1]:
            drughubs.append(drug)
    prothubs = []
    for prot in pdegrees.keys():
        if pdegrees[prot] >= hubsvalues[-1]:
            prothubs.append(prot)
    print('Number of drug-hubs = %d' % len(drughubs))
    print('Number of target-hubs = %d' % len(prothubs))

    print('Number of drugs with only one connection = %d' % len(oneedge))
    print('Number of drugs with less than 4 connections = %d' % len(fewedges))

    return ddegrees, pdegrees, drughubs, prothubs, hubsvalues


def plot_degree_dist(counter):

    deg, cnt = zip(*counter.items())
    plt.bar(deg, cnt)
    plt.title("Degree distribution", fontsize=20)
    plt.ylabel('Number of nodes', fontsize=12)
    plt.xlabel('Degree', fontsize=12)
    plt.show()


def count_stats(drugs):
    # Average node degree, statistics of hubs
    node_degree(drugs)

    G = nx.Graph()
    for drug in drugs.values():
        for edge in drug.edges:
            target = edge[1] if isinstance(edge[1], str) else edge[1].upid
            G.add_edge(drug.name, target)

    # Centrality metrics
    bc = nx.betweenness_centrality(G)
    print('Average betweeness centrality: %f' % np.mean(list(bc.values())))

    # Connected component
    maks = 0
    for g in nx.connected_components(G):
        if len(g) > maks:
            maks = len(g)
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
    plot_degree_dist(counter)
    distribution = {}
    for k in counter.keys():
        distribution[k] = counter[k] / nx.number_of_nodes(G)
    entropy = 0
    for p in distribution.values():
        entropy += p * math.log(p)
    entropy *= -1
    print('Network entropy: %f' % entropy)

    # Modularity Index
    # modularity = nx.modularity_matrix(G)
    # mindex = np.trace(modularity) - sum((modularity ** 2).ravel())
    # print('Modularity index: %f' % mindex)


def target_subtype(drugs, type):

    subdrugs = {}
    for drug in drugs.values():
        for e in drug.edges:
            if e[0] == type:
                if drug.name not in subdrugs:
                    subdrugs[drug.name] = objects.Drug(drug.dbid, drug.name, drug.groups, drug.classification, drug.patent_date)
                subdrugs[drug.name].edges.append(e)

    return subdrugs


def drug_subtype(drugs, type):

    subdrugs = {}
    for drug in drugs.values():
        if type in drug.groups:
            subdrugs[drug.name] = drug
    return subdrugs


drugs, proteins = util.parse('/home/marni/Dokumenty/projektowanie_lekow/full_database.xml')

# statistics of all network
# count_stats(drugs)

# statistics of subtype of targets
'''
for type in ['carrier', 'transporter', 'enzyme']:
    print('\n%s-based network' % type)
    count_stats(target_subtype(drugs, type))
'''

# Statistics of subtype of drugs
for type in ['experimental', 'approved', 'investigational', 'illicit', 'withdrawn']:
    print('\nnetwork based on %s drugs' % type)
    count_stats(drug_subtype(drugs, type))




import igraph as ig
import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
import json


# function to load data SNAP file undirected file (assuming header has been removed)
# returns data in an iGraph graph
def loadData(filename):
    data = None
    data = ig.Graph.Read_Ncol(filename, directed=False)
    data = data.simplify()
    return data

# function to check if size of either of the partitions is less than 3
def checkZeroSize(group1, group2):
    num1 = len(group1)
    num2 = len(group2)
    if (num1 < 3) or (num2 < 3):
        return True

# function to check the ratio of sizes of the partition
# returns false if the condition is satisfied 
def ratioCondition(group1, group2):
    num1 = len(group1)
    num2 = len(group2)
    if num1 > num2:
        ratio = float(num2) / num1
    else:
        ratio = float(num1) / num2
    # if ratio is greater than 0.4 then split
    if ratio > 0.4:
        return False
    else:
        return True

# function to check if the partitioned clusters increase the modularity as compared to the unpartitioned graph
# returns false if the condition is satisfied
def modularityCondition(group1, group2, graph):
    oldMod = 0.0
    newMod = getMod(graph, group1)
    if (newMod > oldMod):
        oldMod = newMod
        return False
    else:
        return True

# function to check if group1 and group2 are good partitions for graph
# inputs: group1 are the id's of nodes in partition 1
#         group2 are the id's of nodes in partition 2
#         graph the unpartitioned graph
# returns false if the groups give a good partition
def recCondition(group1, group2, graph):
    # if size too small don't split
    if checkZeroSize(group1, group2):
        return True
    c1 = False
    # c1 = ratioCondition(group1, group2)
    c2 = modularityCondition(group1, group2, graph)
    if c1:
        return True
    else:
        if c2:
            return True
        return False


# calculates the modularity of a graph given the set of labels
# assumes an igraph graph input and an iterable storing graph label (membership) (integers)
def getMod(graph, labels):
    labeled = [1 if int(b) in labels else 0 for b in graph.vs['name']]
    return graph.modularity(labeled)


# writes graph data and labels into a .json file
# assumes an igraph graph input, an iterable storing graph labels (membership)
# and a set/list with the select cluster labels that you want to extract
def writeJSON(graph, labels, clusters, out='data.json'):
    nodes = [{'name': str(label[0]), 'group': label[1]} for label in enumerate(labels) if label[1] in clusters]
    cut_nl = {}
    for node_tup in enumerate(nodes):
        temp = int(node_tup[1]['name'])
        cut_nl[temp] = node_tup[0]
    el = graph.get_edgelist()
    keys = cut_nl.keys()
    links = [{"source": cut_nl[edge[0]], "target": cut_nl[edge[1]]} for edge in el if edge[0] in keys and edge[1] in keys]
    data = {}
    data['nodes'] = nodes
    data['links'] = links
    with open(out, 'w') as outfile:
        json.dump(data, outfile)


# function to process the graph
def processGraph(graph, l):
    # retrieve the laplacian and calculate the eigenvector
    laplacian = graph.laplacian(normalized=True)
    e_val, e_vec = np.linalg.eigh(laplacian)

    # retrieve the second and third eigenvectors
    idx = e_val.argsort()
    e_val = e_val[idx]
    if (len(e_val) < 2):
        l.append(graph.vs['name'])
        return
    e_vec = e_vec[:, idx]
    e_vec_2 = e_vec[:, 1]
    sorted_idx = e_vec_2.argsort()
    sorted_e_vec_2 = e_vec_2[sorted_idx]

    # separate graph vertices according the the second eigenvector
    x2_pos = np.where(sorted_e_vec_2 > 0)[0]
    x2_neg = np.where(sorted_e_vec_2 <= 0)[0]

    # determine whether the cut should be made
    tup = recCondition(x2_pos, x2_neg, graph)
    cond = tup

    if cond:
        l.append([graph.vs[i]['name'] for i in sorted_idx])
        return

    if (len(e_val) < 3):
        l.append([graph.vs[i]['name'] for i in x2_pos])
        l.append([graph.vs[i]['name'] for i in x2_neg])
        return

    # calculate cuts communities based on the thir eigenvector and determine
    # whether cut should be made
    sub_pos = graph.subgraph(x2_pos, implementation="create_from_scratch")
    sub_neg = graph.subgraph(x2_neg, implementation="create_from_scratch")

    sorted_e_vec_3 = e_vec[:, 2][sorted_idx]

    x3_pos = np.where(sorted_e_vec_3 > 0)[0]
    x3_neg = np.where(sorted_e_vec_3 <= 0)[0]
    pos_pos = np.intersect1d(x2_pos, x3_pos)
    pos_neg = np.intersect1d(x2_pos, x3_neg)
    neg_pos = np.intersect1d(x2_neg, x3_pos)
    neg_neg = np.intersect1d(x2_neg, x3_neg)

    tup1 = recCondition(pos_pos, pos_neg, sub_pos)
    tup2 = recCondition(neg_pos, neg_neg, sub_neg)
    c1 = tup1
    c2 = tup2

    # recursively continue making cuts until stopping conditions are met
    if c1 and c2:
        l.append([graph.vs[i]['name'] for i in x2_pos])
        l.append([graph.vs[i]['name'] for i in x2_neg])
    elif c1:
        # mod_neg = tup2[1]
        l.append([graph.vs[i]['name'] for i in x2_pos])
        sub_neg_pos = graph.subgraph(neg_pos, implementation="create_from_scratch")
        processGraph(sub_neg_pos, l)
        sub_neg_neg = graph.subgraph(neg_neg, implementation="create_from_scratch")
        processGraph(sub_neg_neg, l)
    elif c2:
        # mod_pos = tup1[1]
        l.append([graph.vs[i]['name'] for i in x2_neg])
        sub_pos_pos = graph.subgraph(pos_pos, implementation="create_from_scratch")
        processGraph(sub_pos_pos, l)
        sub_pos_neg = graph.subgraph(pos_neg, implementation="create_from_scratch")
        processGraph(sub_pos_neg, l)
    else:
        sub_pos_pos = graph.subgraph(pos_pos, implementation="create_from_scratch")
        sub_pos_neg = graph.subgraph(pos_neg, implementation="create_from_scratch")
        sub_neg_pos = graph.subgraph(neg_pos, implementation="create_from_scratch")
        sub_neg_neg = graph.subgraph(neg_neg, implementation="create_from_scratch")
        processGraph(sub_neg_pos, l)
        processGraph(sub_neg_neg, l)
        processGraph(sub_pos_pos, l)
        processGraph(sub_pos_neg, l)
    return


if __name__ == '__main__':
    # assumes that the header has been removed
    graph = loadData('fb_graph.txt')
    l = []
    processGraph(graph, l)

    # build membership
    nl = [0]*len(graph.vs)
    count = 0
    for label in l:
        for i in label:
            nl[int(i)] = count
        count += 1
    v = ig.VertexClustering(graph, membership=nl)

    # plot top 4 clustering coeffs in D3
    subs = v.subgraphs()
    clustering_coeffs = []
    for sub in enumerate(subs):
        clustering_coeffs += [(sub[0], sub[1].transitivity_undirected())]

    def getKey(item):
        return item[1]

    clustering_coeffs = sorted(clustering_coeffs, key=getKey, reverse=True)
    clusters = [item[1][0] for item in enumerate(clustering_coeffs) if item[0] < 4]
    writeJSON(graph, v.membership, clusters)

    # plot using igraph kk algorithm
    color_list = ig.known_colors.keys()
    shuffle(color_list)
    # color_list = [
    #     'red',
    #     'blue',
    #     'green',
    #     'cyan',
    #     'pink',
    #     'orange',
    #     'grey',
    #     'yellow',
    #     'white',
    #     'black',
    #     'purple',
    #     'pink',
    #     'maroon',
    #     'navy',
    #     'tan',
    #     'turquoise'
    # ]
    layout = graph.layout("kk")
    plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in v.membership], target="test2.png", vertex_size=5)

    # calculate louvain community and plot
    # dend = graph.community_multilevel()
    # plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in dend.membership], target="test3.png", vertex_size=5)

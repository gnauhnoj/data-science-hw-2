import igraph as ig
import numpy as np
import matplotlib.pyplot as plt
from random import shuffle
import json

def loadData(filename):
    data = None
    data = ig.Graph.Read_Ncol(filename, directed=False)
    data = data.simplify()
    return data

def checkZeroSize(group1,group2):
    num1 = len(group1)
    num2 = len(group2)
    if (num1<3) or (num2<3):
        return True

def ratioCondition(group1, group2):
    num1 = len(group1)
    num2 = len(group2)
    if num1 > num2:
        ratio = float(num2) / num1
    else:
        ratio = float(num1) / num2
    if ratio < 0.4:
        return True
    else:
        return False

def modularityCondition(group1, group2, graph):
    oldMod = 0.0
    newMod = getMod(graph,group1)
    if (newMod >= oldMod):
        print 'modularity splitting', oldMod, newMod
        oldMod = newMod
        return False
    else:
        print 'modularity not splitting', oldMod, newMod
        return True

def recCondition(group1, group2, graph):
    # if size too small don't split
    if checkZeroSize(group1,group2):
        return True
    c1 = ratioCondition(group1,group2)
    c2 = modularityCondition(group1,group2,graph)
    if c1:
        print 'not splitting'
        return True
    else:
        if c2:
            return True
        print 'splitting'
        return False

def getMod(graph, labels):
    # print len(labels)
    # g2 = graph.subgraph(labels, implementation="create_from_scratch")
    labeled = [1 if int(b) in labels else 0 for b in graph.vs['name']]
    # v = ig.VertexClustering(g2)
    # print v.membership
    return graph.modularity(labeled)


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
    print len(nodes)
    data['links'] = links
    with open(out, 'w') as outfile:
        json.dump(data, outfile)


def processGraph(graph, l):
    laplacian = graph.laplacian(normalized=True)
    e_val, e_vec = np.linalg.eigh(laplacian)
    idx = e_val.argsort()
    e_val = e_val[idx]

    if (len(e_val) < 2):
        print 'Done - not enough eval 2', len(e_val)
        l.append(graph.vs['name'])
        return

    e_vec = e_vec[:, idx]
    e_vec_2 = e_vec[:, 1]
    sorted_idx = e_vec_2.argsort()
    sorted_e_vec_2 = e_vec_2[sorted_idx]

    x2_pos = np.where(sorted_e_vec_2 > 0)[0]
    x2_neg = np.where(sorted_e_vec_2 <= 0)[0]

    tup = recCondition(x2_pos, x2_neg, graph)
    cond = tup

    if cond:
        print 'Done - x2 splitt', len(x2_pos), len(x2_neg)
        l.append([graph.vs[i]['name'] for i in sorted_idx])
        return

    if (len(e_val) < 3):
        print 'Done - not enough eval 3', len(e_val)
        l.append([graph.vs[i]['name'] for i in x2_pos])
        l.append([graph.vs[i]['name'] for i in x2_neg])
        return

    sub_pos = graph.subgraph(x2_pos, implementation="create_from_scratch")
    sub_neg = graph.subgraph(x2_neg, implementation="create_from_scratch")

    sorted_e_vec_3 = e_vec[:, 2][sorted_idx]

    x3_pos = np.where(sorted_e_vec_3 > 0)[0]
    x3_neg = np.where(sorted_e_vec_3 <= 0)[0]
    pos_pos = np.intersect1d(x2_pos, x3_pos)
    pos_neg = np.intersect1d(x2_pos, x3_neg)
    neg_pos = np.intersect1d(x2_neg, x3_pos)
    neg_neg = np.intersect1d(x2_neg, x3_neg)

    print 'lengths l3 cuts', len(pos_pos), len(pos_neg), len(neg_pos), len(neg_neg)

    tup1 = recCondition(pos_pos, pos_neg, sub_pos)
    tup2 = recCondition(neg_pos, neg_neg, sub_neg)
    c1 = tup1
    c2 = tup2

    if c1 and c2:
        print 'Done - x3 split', len(pos_pos), len(pos_neg), len(neg_pos), len(neg_neg)
        l.append([graph.vs[i]['name'] for i in x2_pos])
        l.append([graph.vs[i]['name'] for i in x2_neg])
    elif c1:
        # mod_neg = tup2[1]
        print 'Done - x3 neg split', len(pos_pos), len(pos_neg)
        l.append([graph.vs[i]['name'] for i in x2_pos])
        sub_neg_pos = graph.subgraph(neg_pos, implementation="create_from_scratch")
        processGraph(sub_neg_pos, l)
        sub_neg_neg = graph.subgraph(neg_neg, implementation="create_from_scratch")
        processGraph(sub_neg_neg, l)
    elif c2:
        # mod_pos = tup1[1]
        print 'Done - x3 pos split', len(neg_pos), len(neg_neg)
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
    nl = [0]*len(graph.vs)
    count = 0
    for label in l:
        # if (len(label)<3):
        #     community = 0
        # elif (len(label)<15):
        #     community = 1
        # else:
        #     community = count
        #     count += 1
        for i in label:
            nl[int(i)] = count
        count += 1

    v = ig.VertexClustering(graph, membership=nl)
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
    # print 'starting layout'
    layout = graph.layout("kk")
    # print 'finishing layout'
    plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in v.membership], target="test2.png", vertex_size=5)

    # plot any communities of size < 1000
    subs = v.subgraphs()
    clustering_coeffs = []
    for sub in subs:
        clustering_coeffs += [sub.transitivity_undirected()]
    clusters = [item[0] for item in enumerate(l) if len(item[1]) < 1000]
    writeJSON(graph, v.membership, clusters)

    # louvain
    # dend = graph.community_multilevel()
    # plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in dend.membership], target="test3.png", vertex_size=5)

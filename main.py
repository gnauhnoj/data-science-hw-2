import igraph as ig
import numpy as np
import matplotlib.pyplot as plt


def loadData(filename):
    data = None
    data = ig.Graph.Read_Ncol(filename, directed=False)
    data = data.simplify()
    return data


def recCondition(group1, group2, graph, oldMod):
    # num1 = len(group1)
    # num2 = len(group2)
    newMod = getMod(graph, group1)
    if (newMod > oldMod):
        print 'splitting', oldMod, newMod
        oldMod = newMod
        return (False, oldMod)
    else:
        print 'not splitting', oldMod, newMod
        return (True, oldMod)
    # if num1 > num2:
    #     ratio = float(num2) / num1
    # else:
    #     ratio = float(num1) / num2
    # if ratio < 0.1:
    #     return True
    # else:
    #     return False


def getMod(graph, labels):
    # print len(labels)
    # g2 = graph.subgraph(labels, implementation="create_from_scratch")
    labeled = [1 if int(b) in labels else 0 for b in graph.vs['name']]
    # v = ig.VertexClustering(g2)
    # print v.membership
    return graph.modularity(labeled)


def processGraph(graph, l, mod=0):
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
    cond = tup[0]

    if cond:
        print 'Done - x2 splitt', len(x2_pos), len(x2_neg)
        l.append([graph.vs[i]['name'] for i in sorted_idx])
        return
    else:
        mod = tup[1]
        print 'mod change', mod
    # else:
    #     mod = tup[1]
    #     print mod
    #     sub_x2_pos = graph.subgraph(x2_pos, implementation="create_from_scratch")
    #     processGraph(sub_x2_pos, l, mod)
    #     sub_x2_neg = graph.subgraph(x2_neg, implementation="create_from_scratch")
    #     processGraph(sub_x2_neg, l, mod)
    # return

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

    mod_pos = mod
    mod_neg = mod
    tup1 = recCondition(pos_pos, pos_neg, sub_pos, mod_pos)
    tup2 = recCondition(neg_pos, neg_neg, sub_neg, mod_neg)
    c1 = tup1[0]
    c2 = tup2[0]

    if c1 and c2:
        print 'Done - x3 split', len(pos_pos), len(pos_neg), len(neg_pos), len(neg_neg)
        l.append([graph.vs[i]['name'] for i in x2_pos])
        l.append([graph.vs[i]['name'] for i in x2_neg])
    elif c1:
        mod_neg = tup2[1]
        print 'Done - x3 neg split', len(pos_pos), len(pos_neg)
        l.append([graph.vs[i]['name'] for i in x2_pos])
        sub_neg_pos = graph.subgraph(neg_pos, implementation="create_from_scratch")
        processGraph(sub_neg_pos, l)
        sub_neg_neg = graph.subgraph(neg_neg, implementation="create_from_scratch")
        processGraph(sub_neg_neg, l)
    elif c2:
        mod_pos = tup1[1]
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
        for i in label:
            nl[int(i)] = count
        count += 1
    v = ig.VertexClustering(graph, membership=nl)
    color_list = [
        'red',
        'blue',
        'green',
        'cyan',
        'pink',
        'orange',
        'grey',
        'yellow',
        'white',
        'black',
        'purple',
        'pink',
        'maroon',
        'navy',
        'tan',
        'turquoise'
    ]
    # print 'starting layout'
    layout = graph.layout("kk")
    # print 'finishing layout'
    plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in v.membership], target="test.png", vertex_size=5)

    # dend = graph.community_multilevel()
    # plot = ig.plot(graph, layout=layout, vertex_color=[color_list[x] for x in dend.membership], target="test3.png", vertex_size=5)

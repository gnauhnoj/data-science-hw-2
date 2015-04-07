import igraph as ig
import numpy as np
import matplotlib.pyplot as plt


def loadData(filename):
    data = None
    data = ig.Graph.Read_Ncol(filename, directed=False)
    data = data.simplify()
    return data


if __name__ == '__main__':
    # assumes that the header has been removed
    # graph = loadData('yt_graph2.txt')
    # graph = loadData('fb_graph.txt')
    graph = loadData('enron_graph.txt')
    # laplacian = graph.laplacian(normalized=True)
    # e_val, e_vec = np.linalg.eigh(laplacian)
    # idx = e_val.argsort()
    # e_val = e_val[idx]
    # e_vec = e_vec[:, idx]
    # e_vec_2 = e_vec[:, 1]

    # sorted_idx = e_vec_2.argsort()
    # sorted_e_vec_2 = e_vec_2[sorted_idx]
    # group1 = []
    # group2 = []
    # for a in sorted_idx:
    #     if e_vec_2[a] > 0:
    #         group1.append(a)
    #     else:
    #         group2.append(a)
    # sorted_e_vec_3 = e_vec[:, 2][sorted_idx]
    # sorted_e_vec_4 = e_vec[:, 3][sorted_idx]
    # sorted_e_vec_5 = e_vec[:, 4][sorted_idx]

    # plt.plot(sorted_e_vec_2)
    # plt.plot(sorted_e_vec_3)
    # plt.plot(sorted_e_vec_4)
    # plt.plot(sorted_e_vec_5)

    dend = graph.community_multilevel()
    ig.plot(graph, "test.png")

    # x = np.array([[3, -1, -1, 0, -1, 0], [-1, 2, -1, 0, 0, 0], [-1, -1, 3, -1, 0, 0], [0, 0, -1, 3, -1, -1], [-1, 0, 0, -1, 3, -1], [0, 0, 0, -1, -1, 2]])
    # e_val, e_vec = np.linalg.eig(x)

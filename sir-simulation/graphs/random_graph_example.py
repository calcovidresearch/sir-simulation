import networkx as nx
import matplotlib.pyplot as plt

def random_graph(N, E):
    '''
    Creates adjacency list representation of a random graph with N nodes and E edges
    '''
    G = nx.gnm_random_graph(N, E)
    draw_graph(G)
    return nx.convert.to_dict_of_dicts(G)
def draw_graph(G):
    '''
    Draws networkx graph G
    '''
    nx.draw_networkx(G)
    plt.show()
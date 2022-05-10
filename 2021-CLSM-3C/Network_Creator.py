#creates the networks used in this project

import networkx as nx
import numpy as np

##############################################################################
#NETWORK CREATION FUNCTIONS
##############################################################################

def random(n,p):
    #random network
    #n = number of nodes
    #p = probability of connection between two nodes
    return nx.gnp_random_graph(n,p)
    
def ba(n,m):
    #Barabasi-Albert network
    #n = number of nodes
    #m = number of nodes each new node connects to
    return nx.barabasi_albert_graph(n,m)

def caveman(l,k,p):
    #relaxed caveman network
    #l = number of cliques
    #k = number of members per clique
    #p = probability of a node redrawing a connection to connect with another clique
    return nx.relaxed_caveman_graph(l,k,p)

##############################################################################
#NETWORK GENERATOR
##############################################################################

def gen_network(network_func, save=False, filename="network", draw=False, **kwargs):
    """
    Arguments:
        network_func (functions): network generating function
        filename (string): file directory/name for resulting network
        save (bool, optional): If true, saves as .npy file
        draw (bool, optional): If true, plots the network
        **kwargs: parameters for network generating function

    Outputs:
        network_array (numpy array): Array representing network nodes/edges
        edge_density (number): Edge density of network = number of edges/number of nodes
    """

    #get the network
    network = network_func(**kwargs)

    #find the edge density
    edge_density = nx.density(network)

    #convert to numpy array
    network_array = nx.to_numpy_array(network)

    #plot the network
    if (draw == True):
        nx.draw(network, node_size=10)

    #save network as numpy array
    if (save == True):
        np.save(filename + ".npy",network_array)

    return network_array, edge_density

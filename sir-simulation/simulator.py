import networkx as nx
import matplotlib.pyplot as plt
import random
import argparse

def draw_graph(G):
    '''
    Draws networkx graph G
    '''
    nx.draw_networkx(G)
    plt.show()
def getNeighbors(G, n):
    '''
    Returns set/dict of neighbors for node n
    '''
    return G[n]
def random_graph(N, E):
    '''
    Creates adjacency list representation of a random graph with N nodes and E edges
    '''
    G = nx.gnm_random_graph(N, E)
    draw_graph(G)
    return nx.convert.to_dict_of_dicts(G)

def simulate(G, model_params, policy):
    '''
    Randomly simulate spread of a disease using SIR model.
    model_params contains the following:
        S: initial number of susceptible people
        I: initial number of infected
        R: initial removed/dead
        t_i: time spent infected before removal
        T: total timesteps for the simulation
        W: edge weights dictionary
    '''

    V = list(G.keys())
    states = {}
    infected = []
    s_states = [model_params['S']]
    i_states = [model_params['I']]
    r_states = [model_params['R']]
    i = 0
    while i < model_params['S']:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = -1
            i+=1
    i = 0
    while i < model_params['I']:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = model_params['t_i']
            infected.append(random_person)
            i+=1
    i = 0
    while i < model_params["R"]:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = 0
            i+=1
    for t in range(model_params['T']): # TODO: for large N and T this would be slow, percolation?
        for infected_person in infected:
            states[infected_person] -= 1
            if states[infected_person] == 0:
                model_params['R'] += 1
                infected.remove(infected_person) # TODO: this removal is potentially slow
                model_params['I'] -= 1
            for neighbor in getNeighbors(G, infected_person):
                if states[neighbor] == -1 and random.random() < model_params['W'][(infected_person, neighbor)]:
                    model_params['S'] -= 1
                    model_params['I'] += 1
                    states[neighbor] = model_params['t_i']
                    infected.append(neighbor)
        s_states.append(model_params['S'])
        i_states.append(model_params['I'])
        r_states.append(model_params['R'])
        policy(model_params, G, t)
    return (s_states, i_states, r_states)



if __name__ == "__main__":


    G = random_graph(100, 90)
    weights = {}
    for v in G:
        for n in getNeighbors(G,v):
            weights[(v, n)] = 0.5

    model_params = {
        "N" : 100,
        "S" : 95,
        "I" : 5,
        "R" : 0,
        "t_i" : 22,
        "T" : 30,
        "W" : weights
    }
    def no_policy(model_params, G, t):
        return
    states = simulate(G, model_params, no_policy)




    time_states = [i for i in range(model_params["T"] + 1)]

    plt.plot(time_states, states[0], label='susceptible')
    plt.plot(time_states, states[1], label='infected')
    plt.plot(time_states, states[2], label='recovered')
    plt.legend()
    plt.show()
    model_params = {
        "N" : 100,
        "S" : 95,
        "I" : 5,
        "R" : 0,
        "t_i" : 22,
        "T" : 30,
        "W" : weights
    }

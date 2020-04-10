import networkx as nx
import matplotlib.pyplot as plt
import random

def random_graph(N, E):
    '''
    Creates adjacency list representation of a random graph with N nodes and E edges
    '''
    G = nx.gnm_random_graph(N, E)
    return nx.convert.to_dict_of_dicts(G)

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


def simulate(N, E, S, I, R, t_i, T, p):
    '''
    Randomly simulate spread of a disease using SIR model.
    N: number of people
    E: number of contacts
    S: initial number of susceptible people
    I: initial number of infected
    R: initial removed/dead
    t_i: time spent infected before removal
    T: total timesteps for the simulation
    p: probability of transmission across a given contact/edge
    '''
    G = (random_graph(N, E))
    V = list(G.keys())
    print (V)
    states = {}
    infected = []
    s_states = [S]
    i_states = [I]
    r_states = [R]
    i = 0
    while i < S:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = -1
            i += 1
    i = 0
    while i < I:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = t_i
            infected.append(random_person)
            i += 1
    i = 0
    while i < R:
        random_person = random.choice(V)
        if random_person not in states:
            states[random_person] = 0
            i += 1
    print (states)
    for t in range(T): # TODO: for large N and T this would be slow, percolation?
        for infected_person in infected:
            states[infected_person] -= 1
            if states[infected_person] == 0:
                R += 1
                infected.remove(infected_person) # TODO: this removal is potentially realy slow
                I -= 1
            for neighbor in getNeighbors(G, infected_person):
                if states[neighbor] == -1 and random.random() < p:
                    S -= 1
                    I += 1
                    states[neighbor] = t_i
                    infected.append(neighbor)
        s_states.append(S)
        i_states.append(I)
        r_states.append(R)
    return (s_states, i_states, r_states)



if __name__ == "__main__":
    N = 1000  # number of people/nodes

    E = 900  # number of contacts (arbitrary)

    S = 995  # number of people initially in the S state

    I = 5  # number of people initially in the I state

    R = 0  # number of people initially in the R state

    p = 0.6

    t_i = 22

    T = 80
    time_states = [i for i in range(T)]
    time_states = [0] + time_states
    states = simulate(N, E, S, I, R, t_i, T, p)

    plt.plot(time_states, states[0], label='susceptible')
    plt.plot(time_states, states[1], label='infected')
    plt.plot(time_states, states[2], label='recovered')
    plt.legend()
    plt.show()
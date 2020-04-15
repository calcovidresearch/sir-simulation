import networkx as nx
import matplotlib.pyplot as plt
import random
import matplotlib.lines as mlines

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
def show_graph(G, pos):
    mapping = {'S': 'blue', 'I': 'green', 'A': 'orange', 'R': 'red', 'D': 'purple'}
    nodes = G.nodes()
    colors = [mapping[G.nodes[n]['state']] for n in nodes]



    ec = nx.draw_networkx_edges(G, pos, alpha=0.2)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=nodes, node_color=colors,
                                with_labels=False, node_size=100)
    S_legend = mlines.Line2D([], [], color='blue', marker='o', markersize=10,linewidth=0,label='susceptible')
    A_legend = mlines.Line2D([], [], color='orange', marker='o',markerSize=10, linewidth=0, label='asymptomatic')
    I_legend = mlines.Line2D([], [], color='green',marker='o', markerSize=10, linewidth=0, label='infected')
    R_legend = mlines.Line2D([], [], color='red',marker='o', markerSize=10, linewidth=0, label='recovered')
    D_legend = mlines.Line2D([], [], color='purple',marker='o', markerSize=10, linewidth=0, label='dead')

    plt.legend(handles=[S_legend, A_legend, I_legend, R_legend, D_legend], loc=0)
    plt.axis('off')
    plt.show(block=False)
    plt.pause(.5)
    plt.close()

def simulate_corona(G, S, A, I, R, D, T, a_i, a_r, die, policy):
    '''
        Randomly simulate spread of the Coronavirus.
            G
            S
            A
            I
            R
            D
            T: timesteps
            a_i(people_obj): returns a probability that an a node transitions to the I state
            a_r(people_obj): returns a probability that a node transitions to the R state
            die(people_obj): returns a probability that this person will die
            policy : policy function to be performed every epoch/timestep

    '''
    D_counts = [D]
    S_counts = [S]
    A_counts = [A]
    I_counts = [I]
    R_counts = [R]
    pos = nx.spring_layout(G, dim=2, k=None, pos=None, fixed=None, iterations=50, weight='weight', scale=1.0)

    for t in range(1, T+1):
        seen = set()
        show_graph(G, pos)
        for node in G.nodes():

            person = G.nodes[node]
            for neighbor in G[node]:
                if node not in seen:

                    if G.nodes[neighbor]['state'] =='S' and (person["state"] == "I" or person["state"] == "A") and random.random() < G[node][neighbor]['weight']:
                        G.nodes[neighbor]['state'] = 'A'
                        A+=1
                        S-=1
                        seen.add(neighbor)

            death_flip = random.random()
            removal_flip = random.random()
            infection_flip = random.random()
            if person['state'] == "I":
                if death_flip < die(person):
                    person['state'] = "D"
                    D += 1
                    I -= 1
                elif removal_flip < a_r(person):
                    R += 1
                    I -= 1
                    person['state'] = "R"

            elif person['state'] == 'A':
                if removal_flip < a_r(person):
                    R += 1
                    A -= 1
                    person['state'] = "R"
                elif infection_flip < a_i(person):
                    I += 1
                    A -= 1
                    person['state'] = "I"

        D_counts.append(D)
        S_counts.append(S)
        A_counts.append(A)
        I_counts.append(I)
        R_counts.append(R)





        policy(G)
    return (S_counts, A_counts, I_counts, R_counts, D_counts)




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

def initialize_state():
    '''
            G
            S
            A
            I
            R
            D
            T: timesteps
            a_i(people_obj): returns a probability that an a node transitions to the I state
            a_r(people_obj): returns a probability that a node transitions to the R state
            die(people_obj): returns a probability that this person will die
            policy : policy function to be performed every epoch/timestep
    '''
    G = nx.relaxed_caveman_graph(10, 5, .2)

    attributes = {}
    for node in G.nodes():
        attributes[node] = {
            'age' : int(random.normalvariate(38, 13)),
            'state' : "S",
            'underlying_condition' : False,
            'time_infection' : -1,
            'time_removal' : -1,
            'contacts' : [],
            'dead' : False
        }
    nx.set_node_attributes(G, attributes)



    for node in random.sample(G.nodes(), 3):
        G.nodes[node]["state"] = "I"

    for node in G.nodes():
        for neighbor in G[node]:
            G.add_edge(node, neighbor, weight=random.random())



    T = 25
    def a_i(person):
        return 0.2
    def a_r(person):
        return 1.0/14.0
    def die(person):
        return 0.03
    def policy(G):
        return

    lines = simulate_corona(G, 47, 0, 3, 0, 0, T, a_i, a_r, die, policy)
    time_steps = [i for i in range(T+1)]

    plt.plot(time_steps, lines[0], label="susceptible")
    plt.plot(time_steps, lines[1], label="asymptomatic")
    plt.plot(time_steps, lines[2], label="infected")
    plt.plot(time_steps, lines[3], label="removed")
    plt.plot(time_steps, lines[4], label="dead")
    plt.legend()
    plt.show()



if __name__ == "__main__":


    initialize_state()
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
class Person:
    def __init__(self, ageFunction=random.random, underlying_condition=False, time_of_infection = -1, time_of_removal = -1, contacts={}, state="S", dead=False):
        self.age = ageFunction()
        self.state = state
        self.underlying_condition = underlying_condition
        self.time_of_infection = time_of_infection
        self.time_of_removal = time_of_removal
        self.contacts = contacts
        self.dead = dead
    def age(self):
        return self.age
    def has_underlying_condition(self):
        return self.underlying_condition
    def time_of_infection(self):
        return self.time_of_infection()
    def getContacts(self):
        return self.contacts
    def time_of_removal(self):
        return self.time_of_removal
    def getState(self):
        return self.state
    def setState(self, new):
        self.state = new
    def isDead(self):
        return self.dead
    def setDead(self):
        self.dead = True


def simulate_corona(model_state, T, a_i, a_r, die, policy):
    '''
        Randomly simulate spread of the Coronavirus.
            model_state: dictionary of all the state variables
                G: contact graph
                R_state_list: list of people obj's in the R state
                S_state_list: list of people obj's in the S state
                I_state_list: list of people obj's in the I state
                A_state_list: list of people obj's in the A state
                D_state_list: list of people obj's in the that are dead
                people: dictionary of people objects {vector : people obj}
                deaths : list of dead people obj's
                s_a: dictionary of probabilities an s node will become an a node given a contact
            T: timesteps
            a_i(people_obj): returns a probability that an a node transitions to the I state
            a_r(people_obj): returns a probability that a node transitions to the R state
            die(people_obj): returns a probability that this person will die
            policy : policy function to be performed every epoch/timestep

    '''
    deathCounts = [len(model_state["D_state_list"])]
    S_counts = [len(model_state["S_state_list"])]
    A_counts = [len(model_state["A_state_list"])]
    I_counts = [len(model_state["I_state_list"])]
    R_counts = [len(model_state["R_state_list"])]
    D_counts = [len(model_state["D_state_list"])]
    for t in range(1, T+1):
        seen = set()
        for vertex in model_state["people"].keys():
            person = model_state["people"][vertex]
            contacts = []
            for neighbor in getNeighbors(G, vertex):

                if person.getState() == "A" or person.getState() == "I" and person not in seen:
                    if neighbor.getState() == "S" and random.random() < model_state["s_a"][neighbor][vertex]:
                        neighbor.setState("A")
                        seen.add(neighbor)
                        model_state["S_state_list"].remove(neighbor)
                        model_state["A_state_list"].append(neighbor)

                contacts.append(neighbor)
            probability_of_death = die(person)
            if person.getState() == "A" and random.random() < probability_of_death:
                model_state['A_state_list'].remove(vertex)
                person.setDead()
                model_state['D_state_list'].append(vertex)
            if person.getState() == "I" and random.random() < probability_of_death:
                model_state['I_state_list'].remove(vertex)
                model_state['D_state_list'].append(vertex)

                person.setDead()
            if person.getState() == "A" or person.getState() == 'I' and random.random() < a_r(person):
                model_state["A_state_list"].remove(vertex)
                model_state["R_state_list"].append(vertex)
            if person.getState() == "A" and random.random() < a_i(person):
                model_state["A_state_list"].remove(vertex)
                model_state["I_state_list"].append(vertex)
            person.contacts[t] = contacts
        deathCounts = deathCounts + [len(model_state["D_state_list"])]
        S_counts = S_counts + [len(model_state["S_state_list"])]
        A_counts = A_counts + [len(model_state["A_state_list"])]
        I_counts = I_counts + [len(model_state["I_state_list"])]
        R_counts = R_counts + [len(model_state["R_state_list"])]
        D_counts = D_counts + [len(model_state["D_state_list"])]
        policy(model_state)




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

def initialize_state(graph_making_function):
    '''

    :param graph_making_function:
    :return:
    '''

if __name__ == "__main__":


    G = random_graph(100, 90)

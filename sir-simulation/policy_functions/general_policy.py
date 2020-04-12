def policy_function(model_params, G, t):
    '''
    Takes in a dictionary model_params:
        S: Number of susceptible
        I: Number of infected
        R: Number of recovered
        W: Edge weights for probability of transmission along a edge
        t_i: Duration of infection
        T: Number of timesteps of simulation
    Mutates the dictionary and graph G
    '''
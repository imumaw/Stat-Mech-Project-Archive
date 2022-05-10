#Simulates pandemic on network

import numpy as np

################################################################################
#Death and recovery probability curves
################################################################################

def recovery(inf_length, k_r, x_r, **kwargs):
    """Logistic curve model of recovery
    
    Parameters:
        inf_length (int): length of infection
        k_r (number): sharpness of recovery curve
        x_r (number): average time at which recovery occurs
        **kwargs: additional arguments from death curve
    """
    return 1/(1+(np.exp(-1*k_r*(inf_length-x_r))))

def death(inf_length, k_d, **kwargs):
    """Constant function model of death
    
    Parameters:
        inf_length (int): length of infection
        k_d (number): probability of death
        **kwargs: additional arguments from recovery curve
    """
    return k_d

################################################################################
#Pandemic simulation
################################################################################

def import_network(npy_file):
    """Converts a .npy file to an array. 
    npy_file must be a string containing the directory of the file"""
    return np.load(npy_file)


def pandemic(network, cycles, N_inf, p_infect, p_death=None, p_recovery=None, N_immune=0, **kwargs):
    """Simulates pandemic spreading over a network
    
    Parameters:
        network (array): network in array form. If working directly from a file, first convert using import_network
        cycles (positive int): number of cycles through the network
        N_inf (positive int): number of initially infected nodes
        p_infect (number 0 to 1): probability of infection spreading to uninfected node
        p_death (func, optional): probability of infected node dying. If None, skips this step
        p_recovery (func, optional): probability of infected node recovering. If None, skips this step
        N_immune (positive int, optional): number of initially immune nodes. Defaults 0.
        **kwargs: any additional arguments required for probability curves
        
    Returns:
        active (array)
    """

    N = len(network[0])

    #build active array
    active = np.zeros((cycles+1,N),dtype=int)

    #assign initial infections
    rng = np.random.default_rng()
    infected_indices = rng.choice(N,N_inf,replace=False)
    for i in infected_indices:
        active[0,i] = 1
       
    #assign random immunity
    uninfected = np.where(active[0]==0)[0]
    immune_indices = rng.choice(uninfected,N_immune,replace=False)
    for i in immune_indices:
        active[0,i] = -1
    
    for cycle in range(cycles):
        
        for i in range(N):
            
            column = network[:,i]       #get the person's connections
            activated = active[cycle,i] #check the person's status
            
            #infection process
            if activated == 0:
                for j,connect in enumerate(column):
                    if connect==1 and active[cycle,j]==1 and np.random.rand()<p_infect:
                        active[cycle+1,i] = 1
                        continue
            
            #death/recovery process
            elif activated == 1:
                inf_length = np.sum(column)
                if p_recovery != None:
                    pr = p_recovery(inf_length,**kwargs)
                    if np.random.rand() <= pr:
                        active[cycle+1,i] = -1
                        continue
                if p_death != None:
                    pd = p_death(inf_length,**kwargs)
                    if np.random.rand() <= pd:
                        active[cycle+1,i] = -2
                        continue
                active[cycle+1,i] = 1
            
            #if the node is already marked as recovered or dead
            else:
                active[cycle+1,i] = activated
                
    return network, active

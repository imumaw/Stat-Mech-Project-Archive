#contains all visualizations for pandemic simulation

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
import Network_Creator as create
import network_sim

##############################################################################

def plot_data(pandemic_array, title="", xlabel="Days", ylabel="Population"):
    """Plots results of pandemic spread over time
    
    Parameters:
        pandemic_array (array): results from pandemic()
        title, xlabel, ylabel (strings, optional): plotting data
    """
    
    #get information about network
    population = len(pandemic_array[0])
    days = len(pandemic_array[:,0])
    days_array = np.arange(0, days)
    
    #each status gets a row with index for each day
    data_array = np.zeros((4,days))
    
    #get counts for each category for each day
    for i,row in enumerate(pandemic_array):
        data_array[0][i] = np.count_nonzero(row == 0)    #healthy
        data_array[1][i] = np.count_nonzero(row == 1)    #infected
        data_array[2][i] = np.count_nonzero(row == -1)   #recovered/immune
        data_array[3][i] = np.count_nonzero(row == -2)   #dead
    
    #plot results
    plt.fill_between(days_array,population-data_array[0],population,color="tab:blue",label="Healthy")
    plt.fill_between(days_array,data_array[2]+data_array[1]+data_array[3],data_array[1]+data_array[3],color="tab:cyan",label="Recovered/Immune")
    plt.fill_between(days_array,data_array[1]+data_array[3],data_array[3],color="tab:red",label="Infected")
    plt.fill_between(days_array,data_array[3],color="dimgrey",label="Dead")
    
    #make it look nice
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid()
   
    return data_array

def plot_curves(p_death, p_recovery, length=20, **kwargs):
    """Plots death and recovery curves vs days
    
    Parameters:
        p_death (function): chance of death
        p_recovery (function): chance of recovery
        length (tuple, optional): maximum number of days in plot
        **kwargs (optional): any additional arguments for death/recovery
    """
    
    days = np.arange(0,length,1)
    
    #plot data
    plt.plot(days, p_death(days,**kwargs), label="Chance of death")
    plt.plot(days, p_recovery(days,**kwargs), label="Chance of recovery")
    
    #convert to integer axes
    fig = plt.gcf()
    ax = fig.gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    
    #make it pretty
    plt.legend()
    plt.grid()
    
    return

def edges_vs_infections(p=np.linspace(0,.25,10), cycles=100, N_inf=5, p_infect=.01, n_nodes=500, title="Infection curves for varying edge densities"):
    """Plots the infection curve over time while varying the number of edges (used for testing)
    
    Parameters:
        p (array, opt): array of connection probabilities to cycle through
        cycles (int, opt): number of cycles throught the network
        N_inf (int, opt): number of initially infected nodes
        p_inf (float, 0 to 1, opt): probabability of infection of healthy node
        n_nodes (int, opt): number of nodes in network
        title (string, opt): plot title
        
    Returns: 
        None
    """
    
    #cycle through all probabilities
    for p_i in p:
        
        #creates network and runs simulation
        network, density = create.gen_network(create.random, n=n_nodes, p=p_i)
        network, active = network_sim.pandemic(network, cycles, N_inf, p_infect)
        
        #gets infected nodes at each day
        infected = np.zeros(cycles+1)
        for i,row in enumerate(active):
            infected[i] = np.count_nonzero(row == 1)   #infected
        
        #plot it
        plt.plot(infected,label="{:.3f}".format(density))
    
    #make it pretty
    plt.legend(title="Edge density")
    plt.xlabel("Days")
    plt.ylabel("Number of infected nodes")
    plt.title(title)
    
    return


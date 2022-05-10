import network_sim
import visualize as viz
import Network_Creator as create

import matplotlib.pyplot as plt
import numpy as np

##############################################################################
##############################################################################

#Test plots

"""
plt.figure()
viz.edges_vs_infections(p=np.linspace(0,.25,10), cycles=100, N_inf=5, p_infect=.01, n_nodes=500, title="Infection curves for varying edge densities")

plt.figure()
viz.edges_vs_infections(p=np.linspace(0,.25,10), cycles=100, N_inf=5, p_infect=.005, n_nodes=500, title="Infection curves for varying edge densities")

plt.figure()
viz.edges_vs_infections(p=np.linspace(0,.25,10), cycles=100, N_inf=5, p_infect=.02, n_nodes=500, title="Infection curves for varying edge densities")
"""

##############################################################################
##############################################################################

#Generate random network based on known data

#random
"""
plt.figure()
n, e = create.gen_network(create.random, save=True, filename="/Users/isaiahmumaw/Desktop/random_network", draw=True, n=2496, p=.00641026)
print(e)
"""

#caveman
"""
plt.figure()
n, e = create.gen_network(create.caveman, save=True, filename="/Users/isaiahmumaw/Desktop/caveman_network", draw=True, l=156, k=16, p=0.02)
print(e)
"""

##############################################################################
##############################################################################

#COVID-19 plots

#death/recovery curves
"""
plt.figure()
viz.plot_curves(network_sim.death, network_sim.recovery, length=20, k_r=0.4, x_r=8, k_d=.0917, x_d=6)
"""

#get networks
#random = network_sim.import_network("/Users/isaiahmumaw/Documents/GitHub/2021-CLSM-3C/random_network.npy")
#caveman = network_sim.import_network("/Users/isaiahmumaw/Documents/GitHub/2021-CLSM-3C/caveman_network.npy")

#Sim 1 (normal)
"""
plt.figure()
network, active = network_sim.pandemic(random, 100, N_inf=10, p_infect=0.05, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=0.4, x_r=8, k_d=.000917)
viz.plot_data(active, title="Random Network", xlabel="Days", ylabel="Population")

plt.figure()
network, active = network_sim.pandemic(caveman, 100, N_inf=10, p_infect=0.05, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=0.4, x_r=8, k_d=.000917)
viz.plot_data(active, title="Caveman Network", xlabel="Days", ylabel="Population")
"""

#Sim 2 (higher infectivity)
"""
plt.figure()
network, active = network_sim.pandemic(random, 100, N_inf=10, p_infect=0.25, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=1, x_r=8, k_d=.000917)
viz.plot_data(active, title="Random Network", xlabel="Days", ylabel="Population")

plt.figure()
network, active = network_sim.pandemic(caveman, 100, N_inf=10, p_infect=0.25, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=1, x_r=8, k_d=.000917)
viz.plot_data(active, title="Caveman Network", xlabel="Days", ylabel="Population")
"""

#Sim 3 (higher mortality, lower recovery)
"""
plt.figure()
network, active = network_sim.pandemic(random, 100, N_inf=10, p_infect=0.25, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=1, x_r=10, k_d=.03)
viz.plot_data(active, title="Random Network", xlabel="Days", ylabel="Population")

plt.figure()
network, active = network_sim.pandemic(caveman, 100, N_inf=10, p_infect=0.25, p_death=network_sim.death, p_recovery=network_sim.recovery, N_immune=0, k_r=1, x_r=10, k_d=.03)
viz.plot_data(active, title="Caveman Network", xlabel="Days", ylabel="Population")
"""

plt.show()
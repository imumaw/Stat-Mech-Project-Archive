#Isaiah Mumaw
#test code for ising_model. 
#Does not include all plots, some plots were made by modifying arguments

#import modules
import ising_model as ising
import matplotlib.pyplot as plt
import numpy as np
import math

#supress divide by 0 error output. Sometimes fails but this is the easiest workaround
np.seterr(divide='ignore')

#set up main variables
kb=1
mu=1
T = 0
J= [1, 1, 1]
B = 0
p_init = [0.5, 0.5]

#shapes
shape1D = (1,1,2500)
shape2D = (1,50,50)
shape3D = (13,14,14)    #as close as possible to 2500 particles

#############################################################################
#############################################################################
#PART 1: Testing for equilbrium and plotting magnetization
#############################################################################
#############################################################################

#known probability of up spin for test cases of (T,J,B):
    #(0,0,1): 0.00
    #(0,1,0): 0.00 or 1.00 - takes a long time, will form regions of similar spin first
    #(0,1,-1): 1.00
    #(0,1,1): 0.00
    #(>0,0,1): mixed
    #(inf,0,1): 0.50

#total number of tests (higher number increases accuracy, massively increases runtime)
numberOfTests = 1

#cycles per test
cycles = 20

#magnetizations
avgMag1D = np.empty((numberOfTests,cycles*shape1D[0]*shape1D[1]*shape1D[2]))
avgMag2D = np.empty((numberOfTests,cycles*shape2D[0]*shape2D[1]*shape2D[2]))
avgMag3D = np.empty((numberOfTests,cycles*shape3D[0]*shape3D[1]*shape3D[2]))

#running tests for each dimension
for i in range(numberOfTests):
    
    ogState, endState, magArray, eArray = ising.flipStates(shape1D, cycles=cycles, J=J, B=B, kb=kb, T=T, mu=mu, p=p_init) 
    avgMag1D[i,:] = magArray
    
    ogState, endState, magArray, eArray = ising.flipStates(shape2D, cycles=cycles, J=J, B=B, kb=kb, T=T, mu=mu, p=p_init) 
    avgMag2D[i,:] = magArray
    
    ogState, endState, magArray, eArray = ising.flipStates(shape3D, cycles=cycles, J=J, B=B, kb=kb, T=T, mu=mu, p=p_init) 
    avgMag3D[i,:] = magArray

if numberOfTests > 0:

    print(40*"█"+" 100.0%")

    plt.plot(np.transpose(avgMag1D),color="tab:blue",linewidth=.25)
    plt.plot(np.transpose(avgMag2D),color="tab:orange",linewidth=.25)
    plt.plot(np.transpose(avgMag3D),color="tab:green",linewidth=.25)
    plt.title("Magnetizations for each trial")
    plt.ylabel("Magnetization")
    plt.xlabel("Number of steps")
    
    avgMag1D = np.mean(avgMag1D,axis=0)
    avgMag2D = np.mean(avgMag2D,axis=0)
    avgMag3D = np.mean(avgMag3D,axis=0)
    
    plt.figure()
    ising.plotMagVsTime(avgMag1D)
    ising.plotMagVsTime(avgMag2D)
    ising.plotMagVsTime(avgMag3D)
    
    plt.legend(["1D System","2D System","3D System"])
    plt.title("Magnetization vs. Number of Steps")
   

#############################################################################
#############################################################################
#PART 2: Visualizer
#############################################################################
#############################################################################

cycles = 20

ogState, endState2, magArray, eArray = ising.flipStates(shape2D, cycles=cycles, J=J, B=B, kb=kb, T=T, mu=mu, p=p_init) 
ogState, endState3, magArray, eArray = ising.flipStates(shape3D, cycles=cycles, J=J, B=B, kb=kb, T=T, mu=mu, p=p_init) 

plt.figure()
ising.plotter2D(endState2, layer=0, cmap="hot", title="2D Ising Model Domains")

plt.figure()
ising.plotter3D(endState3, alpha=0.6, upcolor="white", downcolor="black", title="3D Ising Model Domains")
"""

#############################################################################
#############################################################################
#PART 3: Other tests
#############################################################################
#############################################################################

#not used in report due to the way I arranged things, but there is also a magnetization vs B field test

equil = 4
cycles = 10
numTests = 25

#plot magnetization vs Temperature

plt.figure()
ising.plotMagVsTemp(shape1D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, B=B, kb=kb, mu=mu, p=p_init)
ising.plotMagVsTemp(shape2D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, B=B, kb=kb, mu=mu, p=p_init)
ising.plotMagVsTemp(shape3D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, B=B, kb=kb, mu=mu, p=p_init)

#plot with theoretical values
#T = np.linspace(0, 15, 100)
#plt.plot(T,abs(np.tanh(mu*B/kb/T)))
#plt.legend(["1D System","2D System","3D System","Theoretical Magnetization"])

#regular plot
plt.legend(["1D System","2D System","3D System"])


#plot magnetic susceptibility vs Temperature

plt.figure()
ising.plotSusVsTemp(shape1D, 0, 15, B-0.5, B, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, mu=mu, p=p_init)
ising.plotSusVsTemp(shape2D, 0, 15, B-0.5, B, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, mu=mu, p=p_init)
ising.plotSusVsTemp(shape3D, 0, 15, B-0.5, B, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, mu=mu, p=p_init)

plt.legend(["1D System","2D System","3D System"])


#plot energy vs temperature

plt.figure()
e1, s1, t1 = ising.plotEVsTemp(shape1D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, B=B, mu=mu, p=p_init)
e2, s2, t2 = ising.plotEVsTemp(shape2D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, B=B, mu=mu, p=p_init)
e3, s3, t3 = ising.plotEVsTemp(shape3D, 0, 15, equil=equil, numTests=numTests, cycles=cycles, J=J, kb=kb, B=B, mu=mu, p=p_init)
plt.legend(["1D System","2D System","3D System"])

#plot specific heat vs temperature
plt.figure()
plt.plot(t1, s1)
plt.plot(t2, s2)
plt.plot(t3, s3)

plt.xlabel("Temperature")
plt.ylabel("Specific Heat")
plt.title("Specific Heat vs. Temperature")
plt.legend(["1D System","2D System","3D System"])
"""
plt.show()

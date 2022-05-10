import numpy as np
import matplotlib.pyplot as plt

#############################################################################

def isingEnergy(array, B, J, mu):
    """Gets energy of ising model
    
    Parameters:
        array (array): current system of spins
        B, J, mu (numbers): constants involved in calculation (see flipStates for full description)
   
    Returns:
        energy (number): energy of system
    """

    energy = mu * B * np.sum(array)
    
    #accounts for interaction strength J along each axis
    for i in range(len(J)):
        energy += -J[i] * np.sum(np.multiply(array, np.roll(array, 1, axis=i)))
        
    return energy

#############################################################################

def boltzProb(E,kb,T):
    """Gets Boltzmann probability given energy, kb, and T"""
    return np.exp(-1*E/(kb*T))

#############################################################################

def flipStates(shape, cycles=100, energyFunc=isingEnergy, probFunc=boltzProb, J=[1, 1, 1], B=1, kb=1, T=1, mu=1, p=[0.5,0.5]):
    """Probabilistically flips the spin of a two-state system to reach equilibrium.
    Begins using a randomly oriented system, however the probability can be adjusted in the code as needed
    
    Parameters:
        shape (3-tuple of ints): shape of system (for 1 and 2 dimensions, set unwanted dimensions to 1)
        cycles (int, optional): number of times to cycle through system. Defaults to 100.
        energyFunc (function, optional): function for energy given conditions
        probFunc (function, optional): function for probability of each proposed energy transition
        J (list, optional): interaction strength along each axis. Defaults to [1,1,1].
        B (number, optional): magnetic field. Defaults to 1
        kb (number, optional): Boltzmann constant. Defaults to 1
        T (number, optional): temperature. Defaults to 1
        mu (numbers, optional): amount of effect due to magnetic field. Defaults to 1
        p (list of two numbers, optional): probabilities of down/up spin. Should total 1.
    
    Returns:
        ogState (array): initial system
        endState (array): equilibrium system
        magArray (array): magnetization at each step
        eArray (array): energy at each step
    """

    #set up starting and ending arrays. Adjust p values for alternate starting points
    ogState = np.random.choice([-1,1], size=shape, p=p)
    endState = np.copy(ogState)
    
    #set up array of random values
    randArray = np.random.rand(cycles, shape[0], shape[1], shape[2])

    #set up arrays for magnetization and energy at each step
    magArray = np.empty(cycles*shape[0]*shape[1]*shape[2])
    eArray = np.empty(cycles*shape[0]*shape[1]*shape[2])
    i_me = 0

    #cycle through all states 100 times
    for i in range(cycles):
        for j in range(shape[0]):
            for k in range(shape[1]):
                for l in range(shape[2]):

                    #find overall energy of current state
                    og_E = energyFunc(endState, B, J, mu)

                    #find overall energy of flipped state
                    endState[j][k][l] = -endState[j][k][l]
                    flipped_E = energyFunc(endState, B, J, mu)

                    #get difference in energies, flip back to original
                    dE = flipped_E - og_E
                    endState[j][k][l] = -endState[j][k][l]

                    #get probability of transition
                    probability = probFunc(dE,kb,T)
                    val = randArray[i][j][k][l]

                    #check if we can make transition
                    if val < probability:
                        endState[j][k][l] = -endState[j][k][l]

                    #update magnetization, energy arrays
                    magArray[i_me] = np.mean(endState)
                    eArray[i_me] = energyFunc(endState, B, J, mu)
                    i_me += 1
    
    return ogState, endState, magArray, eArray

#############################################################################
#PLOTTING FUNCTIONS
#############################################################################

def plotMagVsTime(array,title="Magnetization vs. Time",xlabel="Number of Steps",ylabel="Magnetization"):
    """Plots the magnetization array from flipStates
    
    Parameters:
        array (array): magnetization array
        title, xlabel, ylabel (strings, optional): plotting information
    
    Returns:
        None
    """
    
    plt.plot(array)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return

#############################################################################

def plotMagVsB(shape, Bmin, Bmax, equil=50, numTests=10, cycles=100, energyFunc=isingEnergy, probFunc=boltzProb, J=[1, 1, 1], kb=1, T=1, mu=1, p=[0.5,0.5], title="Magnetization vs. Magnetic Field",xlabel="B-field",ylabel="Magnitude of Magnetization"):
    """Plots magnetization while varying magnetic field
    
    Parameters:
        shape (tuple): shape of system
        Bmin, Bmax (numbers): maximum and minimum magnetic fields over which to calculate
        equil (number): cycle after which equilbrium is reached (must be less than cycles)
        numTests (int, optional): number of data points to produce
        cycles (int, optional): number of cycles through system (ideally will be enough to reach equilibrium)
        energyFunc (function, optional): energy function
        probFunc (function, optional): probability function
        J (number or list, optional): value of J for each dimension
        kb, T, mu (numbers, optional): constants used in calculation
        title, xlabel, ylabel (strings, optional): plotting information
    
    Returns:
        magArray (array): final magnetization for each test (y axis)
        BArray (array): magnetic field for each test (x axis)
    """
    
    #generate arrays
    BArray = np.linspace(Bmin, Bmax, num=numTests)
    magArray = np.empty(numTests)
    
    equil *= shape[0]*shape[1]*shape[2]
        
    #run tests for all B values
    for i,B in enumerate(BArray):
        ogState, endState, m, e = flipStates(shape, cycles, energyFunc, probFunc, J, B, kb, T, mu, p)
        magArray[i] = np.mean(m[equil:]) #get magnetization for final state
    
    #plot results
    plt.plot(BArray,abs(magArray))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return magArray, BArray

#############################################################################

def plotMagVsTemp(shape, Tmin, Tmax, equil=50, numTests=10, cycles=100, energyFunc=isingEnergy, probFunc=boltzProb, J=[1, 1, 1], B=1, kb=1, mu=1, p=[0.5,0.5], title="Magnetization vs. Temperature",xlabel="Temperature",ylabel="Magnitude of Magnetization"):
    """Plots magnetization while varying temperature
    
    Parameters:
        shape (tuple): shape of system
        Tmin, Tmax (numbers): maximum and minimum temperatures over which to calculate. Must be at least 0.
        equil (number): cycle after which equilbrium is reached (must be less than cycles)
        numTests (int, optional): number of data points to produce
        cycles (int, optional): number of cycles through system (ideally will be enough to reach equilibrium)
        energyFunc (function, optional): energy function
        probFunc (function, optional): probability function
        J (number or list, optional): value of J for each dimension
        B, kb, mu (numbers, optional): constants used in calculation
        title, xlabel, ylabel (strings, optional): plotting information
    
    Returns:
        magArray (array): final magnetization for each test (y axis)
        TArray (array): temperature for each test (x axis)
    """
    
    #generate arrays
    TArray = np.linspace(Tmin, Tmax, num=numTests)
    magArray = np.empty(numTests)
    
    equil = equil*shape[0]*shape[1]*shape[2]
        
    #run tests for all T values
    for i,T in enumerate(TArray):
        ogState, endState, m, e = flipStates(shape, cycles, energyFunc, probFunc, J, B, kb, T, mu, p)
        magArray[i] = np.mean(m[equil:]) #get magnetization for final state
    
    #plot results
    plt.plot(TArray,abs(magArray))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return magArray, TArray

#############################################################################

def plotSusVsTemp(shape, Tmin, Tmax, B0, B1, equil=50, numTests=10, cycles=100, energyFunc=isingEnergy, probFunc=boltzProb, J=[1, 1, 1], kb=1, mu=1, p=[0.5,0.5], title="Magnetic Susceptibility vs. Temperature",xlabel="Temperature",ylabel="Magnetic Susceptibility"):
    """Plots magnetic susceptibility while varying temperature
    
    Parameters:
        shape (tuple): shape of system
        Tmin, Tmax (numbers): maximum and minimum temperatures over which to calculate. Must be at least 0.
        B0, B1 (numbers): initial and final magnetic field
        equil (number): cycle after which equilbrium is reached (must be less than cycles)
        numTests (int, optional): number of data points to produce
        cycles (int, optional): number of cycles through system (ideally will be enough to reach equilibrium)
        energyFunc (function, optional): energy function
        probFunc (function, optional): probability function
        J (number or list, optional): value of J for each dimension
        kb, mu (numbers, optional): constants used in calculation
        title, xlabel, ylabel (strings, optional): plotting information
    
    Returns:
        susArray (array): magnetic susceptibility values (y axis)
        TArray (array): temperature for each test (x axis)
        magArray (array): final magnetization for each test
    """
    
    #generate arrays
    TArray = np.linspace(Tmin, Tmax, num=numTests)
    magArray = np.empty((2,numTests))
    
    equil *= shape[0]*shape[1]*shape[2]
    
    #get magnetization for initial field
    for i,T in enumerate(TArray):
        ogState, endState, m, e = flipStates(shape, cycles, energyFunc, probFunc, J, B0, kb, T, mu, p)
        magArray[0,i] = np.mean(m[equil:]) #get magnetization for current endState
      
    #get magnetization for final field
    for i,T in enumerate(TArray):
        ogState, endState, m, e = flipStates(shape, cycles, energyFunc, probFunc, J, B1, kb, T, mu, p)
        magArray[1,i] = np.mean(m[equil:]) #get magnetization for current endState
    
    #get susceptibility
    susArray = (magArray[1,:]-magArray[0,:])/(B1-B0)
    
    #plot results
    plt.plot(TArray,susArray)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return susArray, TArray, magArray

#############################################################################

def plotEVsTemp(shape, Tmin, Tmax, equil=50, numTests=10, cycles=100, energyFunc=isingEnergy, probFunc=boltzProb, J=[1, 1, 1], B=1, kb=1, mu=1, p=[0.5,0.5], title="Energy vs. Temperature",xlabel="Temperature",ylabel="Energy"):
    """Plots energy while varying temperature. Returns specific heat for future plotting
    
    Parameters:
        shape (tuple): shape of system
        Tmin, Tmax (numbers): maximum and minimum temperatures over which to calculate. Must be at least 0.
        equil (number): cycle after which equilbrium is reached (must be less than cycles)
        numTests (int, optional): number of data points to produce
        cycles (int, optional): number of cycles through system (ideally will be enough to reach equilibrium)
        energyFunc (function, optional): energy function
        probFunc (function, optional): probability function
        J (number or list, optional): value of J for each dimension
        B, kb, mu (numbers, optional): constants used in calculation
        title, xlabel, ylabel (strings, optional): plotting information
    
    Returns:
        energyArray (array): total energy for each test (y axis)
        specHeatArray (array): specific heat for each test (y axis)
        TArray (array): magnetic field for each test (x axis)
    """
    
    #generate arrays
    TArray = np.linspace(Tmin, Tmax, num=numTests)
    energyArray = np.empty(numTests)
    
    #calculate energy for each test
    for i,T in enumerate(TArray):
        ogState, endState, m, e = flipStates(shape, cycles, energyFunc, probFunc, J, B, kb, T, mu, p)
        energyArray[i] = np.mean(e[equil:])
        
    #get specific heat
    specHeatArray = np.gradient(energyArray)
    
    #plot results
    plt.plot(TArray, energyArray)

    #add details
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    
    return energyArray, specHeatArray, TArray

#############################################################################

def plotter2D(matrix, layer=0, cmap="hot", title="2D Ising Model"):

    """Generates heatmap of 2D ising model
    
    Parameters:
        matrix (array): resulting data from flipStates, trimmed to a 2D array if needed
        layer (int, optional): specific layer to plot. For 2D, assuming the shape is in 
                     the form (1,x,y), this can just be 0.
        cmap (string, optional): colormap of plot. Defaults hot (black and white)
        title (string, optional): title of plot
    """
    
    plt.imshow(matrix[layer], cmap=cmap)
    
    plt.title(title)
    
    return

#############################################################################

def plotter3D(matrix, alpha=0.5, upcolor="white", downcolor="black", title="3D Ising Model"):

    """Generates voxel plot of 3D ising model
    
    Parameters:
        matrix (array): resulting data from flipStates
        alpha (number between 0 and 1): opacity of voxels. Higher number means more opaque voxels
        upcolor, downcolor (strings or RGB tuple): color of up and down spins
        title (string, optional): title of plot
    """

    #generate 3D plot
    fig = plt.gcf()
    ax = fig.gca(projection='3d')
    
    #plot up spins (set down spins to 0)
    matrix = matrix*0.5+0.5
    ax.voxels(matrix, facecolor=upcolor, edgecolor=None, alpha=alpha)
    
    #plot down spins (set up spins to 0)
    matrix = matrix-1
    ax.voxels(matrix, facecolor=downcolor, edgecolor=None, alpha=alpha)
    
    plt.title(title)
    
    return

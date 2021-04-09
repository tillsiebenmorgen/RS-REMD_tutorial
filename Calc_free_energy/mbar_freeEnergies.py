from __future__ import division
import os
import numpy as np
import cPickle as pickle
import pymbar




def get_data(data):
    with open(data, "rb") as pfile:
        energies = pickle.load(pfile)
    return energies


def get_potential(energies, data_set1, Ecomponent, framePercentage):
    """
    We calculate the bias that we introduced in each trajectory and save it in the appropiate shape for mbar. 
    Bias = Energy_ReplicaBias - Energy_ReplicaRef

    """
    N_k = []
    TrajBiases = []
    for replikaRef in range(1,replicaNum+1):
        trajNum = replikaRef -1
        biases = []
        for replikaBias in range(1,replicaNum+1):
            # Bias = Energy_ReplicaBias - Energy_ReplicaRef
            bias = energies[trajNum][data_set1][str(replikaBias)][Ecomponent]-energies[trajNum][data_set1][str(replikaRef)][Ecomponent]
            bias = bias[int(framePercentage[0]*len(bias)):int(framePercentage[1]*len(bias))]
            if trajNum == 0:
                N_k.append(len(bias))
            for i in range(len(bias)):
                biases.append(bias[i])
        TrajBiases.append(biases)
    u_kn = np.asarray(TrajBiases)
    return u_kn, N_k

def get_percentages():
    """
    We calculate the biases/free energy values according to the splitting defined in this function in order to calculate uncertainties later 
    """
    Percentages = []
    for start in range(10):
        Percentages.append((start/10., (start+1)/10.))
    return sorted(Percentages, key = lambda x: x[1])

def calculate_free_energies(path, mbarFile, framePercentage):
    """
    We calculate the biases that we introduced by the parameter scaling: U_kn
    And plug these values into mbar
    """
    energies = get_data(path+mbarFile)
    U_kn, N_k = get_potential(energies, data_set1 = "bothUnbiased", Ecomponent = "vdw", framePercentage = framePercentage)
    u_kn = U_kn/kt
    mbar = pymbar.mbar.MBAR(u_kn, N_k)
    result = mbar.getFreeEnergyDifferences()
    free_energy = (kt*result[0],kt*result[1])
    return free_energy

def write_individual_values(framePercentage, path, free_energy):
    """
    Write the individual values to files and pickle them
    """
    calcValues = {}
    scoref = open("free_energies/mbar_"+str(framePercentage[0])+"-"+str(framePercentage[1])+".score", "w")
    scoref.write("{}      {}      {}\n".format(path, str(free_energy[0][0][-1])[0:7],str(free_energy[1][0][-1])[0:7]))
    pfile = path+"free_energies_mbar_"+str(framePercentage[0])+"-"+str(framePercentage[1])+".pickle"
    with open(pfile, "wb") as f:
        pickle.dump(free_energy, f)
    calcValues[path] = (free_energy[0][0],free_energy[1][0])
    pfile = "free_energies/mbar_"+str(framePercentage[0])+"-"+str(framePercentage[1])+".pickle"
    with open(pfile, "wb") as f:
        pickle.dump(calcValues, f)
    return calcValues


def write_all_values(allPercentageCalculated):
    pfileAll = "free_energies/mbar_allPercentages.pickle"
    with open(pfileAll, "wb") as f:
        pickle.dump(allPercentageCalculated, f)

########
# MAIN
########



kt_j =300*1.38064852*10**(-23)
kt=kt_j*0.000239006*6.022140 * 10**(23)


if not os.path.isdir("free_energies"):
    os.makedirs("free_energies")


percentages = get_percentages()
replicaNum = 16
mbarFile = "energies.pickle"
path = "trajectories/"
allPercentageCalculated = {}
Affinities = []
for framePercentage in percentages:
    print 'Fraction of traj that is being evaluated: ', framePercentage
    free_energy = calculate_free_energies(path, mbarFile, framePercentage)
    calcValues = write_individual_values(framePercentage, path, free_energy)
    allPercentageCalculated[framePercentage] = calcValues
    Affinities.append(calcValues[path][0][-1])

write_all_values(allPercentageCalculated)

print 'Binding free energy', np.mean(Affinities[5:]), '+/-', np.std(Affinities[5:])


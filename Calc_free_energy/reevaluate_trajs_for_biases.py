#!/usr/bin/env python

from __future__ import division
import pytraj as pt
import cPickle as pickle
import sander
import parmed as pmd
import parmed.tools as t


def load_trajs(files):
    trajs = []
    for crd in files:
        num = crd.split(".")[1]
        top = "strippedTopologies/system_{}.top".format(int(num))
        traj = pt.iterload(crd, top=str(top), stride = stride)
        trajs.append(traj)
    print "loaded {} trajectories".format(len(trajs))
    return trajs

def get_mask():
	with open("./mask.info", "r") as maskf:
		masks = []
		for l in maskf:
			masks.append(l[:-1])
	return masks[0], masks[1]

def pickle_open(pickle_file):
	picklef = open(pickle_file, "wb")
	return picklef

def calc_energy(traj, mask=None):
	options = sander.pme_input()
	options.cut = 9.0
	files_top_new = ["strippedTopologies/system_"+str(i)+".top" for i in range(1,1+replicas)]
	unbias_energies = {}
	for file_top_new in files_top_new:
		ene = pt.energy_decomposition(traj, top=file_top_new, mm_options=options)
		unbias_energies[file_top_new.split("_")[1].split(".")[0]] = ene
	return unbias_energies

def calc_energies(traj, receptor_mask, ligand_mask):
    filename = traj.filename
    print "Energy calculation of ", filename
    bothUnbiased = calc_energy(traj, "{}|{}".format(receptor_mask, ligand_mask))
    return {"name": traj.filename,
            "bothUnbiased": bothUnbiased}



#########
# MAIN
#########

ligand_mask, receptor_mask = get_mask()
path = "trajectories"
replicas = 16
NumberOfEvaluatedFrames = 50
traj_stride = pt.iterload(path+"/stripped.1.nc", top = "strippedTopologies/system_1.top")
stride = (len(traj_stride) + NumberOfEvaluatedFrames-1) // NumberOfEvaluatedFrames
energies = [] 
pickleFile = pickle_open(path+"/energies.pickle")
trajectoryNames = [path+"/stripped."+str(i)+".nc" for i in range(1,1+replicas)]
trajs = load_trajs(trajectoryNames)
def fn(traj):
	return calc_energies(traj, receptor_mask, ligand_mask)	
energies = map(fn, trajs)
pickle.dump(energies, pickleFile)











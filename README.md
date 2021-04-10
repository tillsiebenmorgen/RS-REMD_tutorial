# Tutorial for RS-REMD

## Tutorial to perform RS-REMD simulations as introduced in the following articles:
---
**Prediction of protein–protein complexes using replica exchange with repulsive scaling**

**Till Siebenmorgen** till.siebenmorgen@tum.de, **Michael Engelhard**, **Martin Zacharias**

https://doi.org/10.1002/jcc.26187

---

**Efficient Refinement and Free Energy Scoring of Predicted Protein–Protein Complexes Using Replica Exchange with Repulsive Scaling**

**Till Siebenmorgen**, **Martin Zacharias**

https://doi.org/10.1021/acs.jcim.0c00853

---

Scripts are written in python 2.7 (for Amber users)

### Example_preparation

Following these steps you can set up a RS-REMD simulation

Perform a setup of the groupfile, modify the Lennard-Jones parameters and start the simulation:  


```
python prepare_files_rs.py
python lj.py parm.top 4 -d 0.00 0.1 0.2 0.5 -e 1.00 0.95 0.90 0.80 -c heated.rst -r :1-66 -l :67-123 -hMassRepartitioning -o system.top | tee lj.log
source rs.run
```





**lj.py (python lj.py -help):**

*parm.top* = input topology file that you have created in advance of your system

*4* = Number of replicas that you want to use, in our case 4; for proper systems between 8 and 16

*-d* = distance parameter for LJ scaling of the individual replica; for my systems with 16 replicas:  `-d 0.00 0.01 0.02 0.04 0.08 0.12 0.16 0.20 0.24 0.28 0.32 0.38 0.44 0.50 0.58 0.68`

*-e* = epsilon parameter for LJ scaling of the individual replica;  for my systems with 16 replicas: `-e 1.00 0.99 0.98 0.97 0.96 0.94 0.92 0.90 0.88 0.86 0.84 0.82 0.80 0.78 0.76 0.74`

*-c* = heated restart file

*-r* = receptor mask, specifies the receptor residues whose interaction should be scaled

*-l* = ligand mask, specifies the receptor residues whose interaction should be scaled

*-o* = outputfile name of the generated replicas


### Calc_free_energy

From generated trajectories we can now calculate free energy differences. The simulations were performed in explicit solvent and the water molecues were stripped.   


```
python reevaluate_trajs_for_biases.py 

python mbar_freeEnergies.py
```




### Installation
You might be able to use your python installation from ambertools. You can also install the necessary packages using anaconda
```
conda install -c conda-forge mdtraj
conda install -c ambermd ambertools
```

For the `pymbar` installation (in order to calculate free energies) you can use
`pip install pymbar`


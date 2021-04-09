#!/usr/bin/env python2

from __future__ import division
import sys
import os
from math import sqrt
import subprocess
import argparse
import numpy as np

import parmed as pmd
import parmed.tools as t
from parmed.amber import AmberMask



def mask_from_atoms(atoms):
    return "@" + ",".join(map(lambda a: str(a.idx + 1), atoms))


def group_by(atoms, keyfunc):
    """
    Return a dictionary with the following format: {keyfunc(atom): [atoms]}, where the
    atoms are grouped by the return value of keyfunc(atom).
    """
    groups = {}
    for atom in atoms:
        key = keyfunc(atom)
        group = groups.setdefault(key, [])
        group.append(atom)
    return groups


def group_by_ljtype(parm, mask):
    """Group atoms in mask by LJ type index."""
    atoms = parm.view[mask]
    groups = group_by(atoms, lambda atom: atom.nb_idx)
    return groups


def add_lj_types(parm, not_ligand_groups, ligand_groups):
    """
    Add new LJ type for every atom type in the ligand, so that the ligand has no LJ type
    in common with the rest.
    """
    assert len(ligand_groups) > 0
    
    for nb_idx, atoms in ligand_groups.iteritems():
        if not nb_idx in not_ligand_groups:
            # don't create a new type when we don't need to i.e. when this type only 
            # appears in the ligand; T: wenn Bedingung erfuellt macht keinen neuen lj type
            continue
        t.addLJType(parm, mask_from_atoms(atoms)).execute()


def change_lj_pairs(parm, receptor_groups, ligand_groups, interaction_func):
    """
    Change LJ pair interaction parameters between receptor and ligand.
    interaction_func(atom1, atom2):
        a function that returns the new interaction parameters (rmin, epsilon)
    """
    assert len(ligand_groups) > 0 and len(receptor_groups) > 0



    for group_r in receptor_groups.itervalues():
        for group_l in ligand_groups.itervalues():           
            # enough to look at first atom (only nb_idx is relevant)
            rmin, epsilon = interaction_func(group_r[0], group_l[0])

            mask_r = "@" + str(group_r[0].idx + 1)
            mask_l = "@" + str(group_l[0].idx + 1)

            t.changeLJPair(parm, mask_r, mask_l, rmin, epsilon).execute()


class InterpolateBase(object):
    def __init__(self, strength, **kwargs):
        self.strength = strength

    method = "no interpolation"
    def __call__(self, atom1, atom2):
        rmin = atom1.rmin + atom2.rmin
        epsilon = sqrt(atom1.epsilon * atom2.epsilon)
        return rmin, epsilon



class InterpolateVolumeScaling(InterpolateBase):
    method = "volume scaling of e"
    def __init__(self, strength, interpolator_prev, d_max, e_min, d_fixed, e_fixed, **kwargs):
        if d_fixed:
            self.d = d_max
        else:
            self.d = strength * d_max
        if e_fixed:
            self.e = e_min
        else:
            self.e = interpolate_lin(1, e_min, strength)

    def __call__(self, atom1, atom2):
        # sometimes happens between two H atoms for example
        if atom1.rmin + atom2.rmin == 0.0:
            return 0.0, 0.0

        new_rmin = (atom1.rmin + atom2.rmin) + self.d
        
        e_scaling = (atom1.rmin + atom2.rmin)**3 / new_rmin**3
        epsilon = sqrt(atom1.epsilon * atom2.epsilon) * e_scaling * self.e

        return new_rmin, epsilon


def interpolate_lin(a, b, fraction):
    return a + (b-a) * fraction

def HMassRepartition(parm):
    t.HMassRepartition(parm).execute()
    
    return parm

def write_topology(parm, args, n=None):
    filename = args.output_file
    if n != None:
        filename = filename.replace(".top", "_{}.top".format(n+1))
    print filename
    parm.write_parm(filename)
    with open(filename + ".txt", "w") as f:
        f.write(str(t.printLJMatrix(parm, "*")))




if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", help="amber topology file")
    parser.add_argument("num_replicas", type=int)
    parser.add_argument("-d", "--distances", nargs="*", type=float, help="max distance that gets added to rmin or list of floats", default=[0])
    parser.add_argument("-e", "--epsilons", nargs="*", type=float, help="min fraction of epsilon or list of floats", default=[1])
    parser.add_argument("--d_fixed", action="store_true")
    parser.add_argument("--e_fixed", action="store_true")
    parser.add_argument("-o", "--output_file", required=False, help="amber topology file", default="fk_replica.top")
    parser.add_argument("-r", "--receptor_mask", required=False, help="receptor mask", default="@1-1984")
    parser.add_argument("-l", "--ligand_mask", required=False, help="ligand mask", default=":LIG")
    parser.add_argument("-c", "--coordinate_file", required=False, help="input .inpcrd or .rst file for translating or solvating", default=None)
    parser.add_argument("-hMassRepartitioning", "--hMassRepartitioning", action='store_true', required=False, help="Performs hMassRepartitioning on the topopoly files", default=None)
    args = parser.parse_args()


    args.not_ligand_mask = "!" + args.ligand_mask
    print "receptor mask", args.receptor_mask
    print "ligand mask", args.ligand_mask
    print "not-ligand mask", args.not_ligand_mask

    interpolation_methods = [InterpolateVolumeScaling]
    Method = interpolation_methods[0]

    ds = args.distances
    es = args.epsilons
    assert(len(ds) == 1 or len(ds) == args.num_replicas)
    assert(len(es) == 1 or len(es) == args.num_replicas)
    if len(ds) > 1:
        args.d_fixed = True
    if len(es) > 1:
        args.e_fixed = True

    print "interpolation method:", Method.method
    print "max distance that gets added:", args.distances
    print "min fraction of epsilon:", args.epsilons
    print "d fixed / not interpolated:", args.d_fixed
    print "e fixed / not interpolated:", args.e_fixed
    if args.hMassRepartitioning is not None:
        print 'Perform HMassRepartitioning'
    # load parameter and possibly coordinate file
    print args.input_file, args.coordinate_file
    parm = pmd.amber.AmberParm(args.input_file, xyz=args.coordinate_file)


    # change the lj types
    receptor_groups = group_by_ljtype(parm, args.receptor_mask)
    ligand_groups = group_by_ljtype(parm, args.ligand_mask)
    not_ligand_groups = group_by_ljtype(parm, args.not_ligand_mask)
    add_lj_types(parm, not_ligand_groups, ligand_groups)

    # go through each replica and change lj parameters
    interpolator_prev = None
    parms = []
    orig_charge = parm.parm_data['CHARGE']
    for n, fraction in enumerate(np.linspace(0, 1, args.num_replicas)):
        d = ds[n] if len(ds) > 1 else ds[0]
        e = es[n] if len(es) > 1 else es[0]

        interpolator = Method(fraction, interpolator_prev=interpolator_prev, 
                              d_max=d, e_min=e,
                              d_fixed=args.d_fixed, e_fixed=args.e_fixed)
        interpolator_prev = interpolator
        print "{}: d={}, e={}".format(n+1, interpolator.d, interpolator.e)
        # change lj parameters
        change_lj_pairs(parm, receptor_groups, ligand_groups, interpolator)
        if args.hMassRepartitioning is not None:
            HMassRepartition(parm)
        write_topology(parm, args, n)


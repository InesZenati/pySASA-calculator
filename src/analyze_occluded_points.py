"""Script to analyze solvent occlusion of a parsed protein structure."""
import math
import os
import sys
import time
from datetime import UTC, datetime

import click
from loguru import logger
from tqdm import tqdm

from parse_pdb import load_radius_json, parse_pdb

def compute_atom_surface(atom_object):
    """
    Compute point counts and surface areas (total, occluded, accessible) of one atom.
    
    """
    atom_infos  = {}
    number_of_points = len(atom_object.sphere.pointlist)
    number_of_occluded_atoms = sum(atom_object.sphere.occluded_points) 
    number_of_accessible = number_of_points - number_of_occluded_atoms

    total_surface = 4 * math.pi * atom_object.sphere.radius ** 2
    occluded_surface = (total_surface / number_of_points) * number_of_occluded_atoms
    accessible_surface = total_surface - occluded_surface

    atom_infos = {
        "number_of_points": number_of_points,
        "number_of_occluded_atoms": number_of_occluded_atoms,
        "number_of_accessible": number_of_accessible,
        "total_surface": total_surface,
        "occluded_surface": occluded_surface,
        "accessible_surface": accessible_surface,
    }
    
    return atom_infos

def compute_residu_surface(residue):
    result = {}
    atom_list = residue.atomlist
    for atom in atom_list:
        atom_result = compute_atom_surface(atom)
        for key, value in atom_result.items():
            result[key] = result.get(key, 0) + value
    return result


def compute_chain_surface(residuelist, chain_id):

    chain_atoms = []
    result = {}
    for residue in residuelist:
        if residue.chain_id == chain_id:
            chain_atoms.extend(residue.atomlist)
    for atom in chain_atoms:
        atom_result = compute_atom_surface(atom)
        for key, value in atom_result.items():
            result[key] = result.get(key, 0) + value
    return result
    
def get_chain_ids(residuelist):
    chain_ids = set()
    for residue in residuelist:
        chain_ids.add(residue.chain_id)
    return sorted(chain_ids)


if __name__ == "__main__":
    radius_table = load_radius_json("data/radius.json")
    protein = parse_pdb("data/1CRN.pdb", "1CRN", radius_table)
    chain_ids = get_chain_ids(protein.residuelist)
    print(f"chain id : {chain_ids}")

    for residue in protein.residuelist:
        res_surface = compute_residu_surface(residue)

    for chain_id in chain_ids:
        chain_surface = compute_chain_surface(protein.residuelist, chain_id)
        print(
            f"Chain {chain_id} total surface : {chain_surface.get('total_surface', 0)} A"
        )

    if protein.residuelist and protein.residuelist[0].atomlist:
        first_atom = protein.residuelist[0].atomlist[0]
        atom_surface = compute_atom_surface(first_atom)
        print(
            f"First atom total surface : {atom_surface['total_surface']} A"
        )
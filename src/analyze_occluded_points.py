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

def compute_atom_surface_infos(atom):
    """
    Compute point counts and surface areas (total, occluded, accessible) of one atom.
    
    """
    atom_infos  = {}
    number_of_points = len(atom.sphere.pointlist)
    number_of_occluded_atoms = sum(atom.sphere.occluded_points) 
    number_of_accessible = number_of_points - number_of_occluded_atoms

    total_surface = 4 * math.pi * atom.sphere.radius ** 2
    occluded_surface = (total_surface / number_of_points) * number_of_occluded_atoms
    accessible_surface = total_surface - occluded_surface

    atom_infos = {
        "number_of_points": number_of_points,
        "number_of_occluded_points": number_of_occluded_atoms,
        "number_of_accessible": number_of_accessible,
        "total_surface": total_surface,
        "occluded_surface": occluded_surface,
        "accessible_surface": accessible_surface,
    }
    
    return atom_infos

def compute_all_atoms_surface_infos(atomlist):
    """Sum point counts and surfaces over a list of atoms."""
    for atom in atomlist:
        atom_result = compute_atom_surface_infos(atom)
    return atom_result


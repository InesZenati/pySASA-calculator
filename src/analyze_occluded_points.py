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


def compute_atom_surface(atom):
    """
    Compute point counts and surface areas (total, occluded, accessible) of one atom.
    
    
    """
    number_points = len(atom.sphere.pointlist)
    number_of_occluded_atoms = sum(atom.sphere.occluded_points) 
    n_accessible = number_points - number_of_occluded_atoms

    total_surface = 4 * math.pi * atom.sphere.radius ** 2
    occluded_surface = (total_surface / number_points) * number_of_occluded_atoms
    accessible_surface = total_surface - occluded_surface

    return {
        "n_points": number_points,
        "n_occluded": number_of_occluded_atoms,
        "n_accessible": n_accessible,
        "total_surface": total_surface,
        "occluded_surface": occluded_surface,
        "accessible_surface": accessible_surface,
    }


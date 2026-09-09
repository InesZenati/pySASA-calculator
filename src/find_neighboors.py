"""Script to find neighboring atoms with a manual double loop."""
import math


def compute_atom_atom_distance(atom_one, atom_two):
    return math.sqrt((atom_one.x - atom_two.x) ** 2
                      + (atom_one.y - atom_two.y) ** 2
                      + (atom_one.z - atom_two.z) ** 2)


def get_max_cutoff(all_atoms):
    max_radius = max(atom.sphere.radius for atom in all_atoms)
    return 2 * max_radius

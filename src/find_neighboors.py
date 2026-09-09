"""Script to find neighboring atoms with a manual double loop."""
import math


def compute_atom_atom_distance(atom_one, atom_two):
    """
    Compute the Euclidean distance between the centers of two atoms.

    Parameters
    ----------
    atom_one : Atom
        The first atom.
    atom_two : Atom
        The second atom.

    Return
    ------
    float :
        The distance between the two atom centers, in Angstroms.
    """
    return math.sqrt((atom_one.x - atom_two.x) ** 2
                      + (atom_one.y - atom_two.y) ** 2
                      + (atom_one.z - atom_two.z) ** 2)


def get_max_cutoff(all_atoms):
    """
    Compute the maximum distance at which two atoms can still occlude
    each other.

    Two atoms can only occlude each other if the sum of their sphere radius
    is greater than the distance between them, so the cutoff is twice the
    largest sphere radius found in the protein.

    Parameters
    ----------
    all_atoms : list
        A list of all Atom objects in the protein.

    Return
    ------
    float :
        The cutoff distance, in Angstroms.
    """
    max_radius = max(atom.sphere.radius for atom in all_atoms)
    return 2 * max_radius
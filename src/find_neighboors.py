"""Script to find neighboring atoms using Biopython's NeighborSearch."""
from Bio.PDB import NeighborSearch


def build_neighbor_search(all_atoms):
    """
    Build a spatial search structure from a list of atoms.

    Parameters
    ----------
    all_atoms : list
        A list of all Atom objects in the protein.

    Return
    ------
    NeighborSearch :
        A NeighborSearch instance built from all_atoms.
    """
    return NeighborSearch(all_atoms)

def get_max_radius(all_atoms):
    """
    Get the largest sphere radius found among all atoms of the protein.

    Parameters
    ----------
    all_atoms : list
        A list of all Atom objects in the protein.

    Return
    ------
    float :
        The largest sphere radius, in Angstroms.
    """
    return max(atom.sphere.radius for atom in all_atoms)

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


def get_atom_cutoff(atom, max_radius):
    """
    Compute the cutoff distance for one atom's neighbor search.

    Two atoms i and j can only occlude each other if the distance between
    them is smaller than the sum of their sphere radius (Ri + Rj). Since we
    do not know Rj in advance, we use the largest possible radius found in
    the protein as an upper bound for Rj.

    Parameters
    ----------
    atom : Atom
        The atom object for which the cutoff is computed.
    max_radius : float
        The largest sphere radius found among all atoms of the protein.

    Return
    ------
    float :
        The cutoff distance for this atom's neighbor search, in Angstroms.
    """
    return atom.sphere.radius + max_radius


def get_atom_neighbors(atom, neighbor_search, max_radius):
    """
    Get the atoms within occlusion distance from an atom.

    Parameters
    ----------
    atom : Atom
        The atom object around which to search neighbors.
    neighbor_search : NeighborSearch
        A NeighborSearch instance built from all the atoms of the protein.
    max_radius : float
        The largest sphere radius found among all atoms of the protein.

    Return
    ------
    list :
        A list of Atom objects that could occlude atom, excluding the atom
        itself.
    """
    atom_cutoff = get_atom_cutoff(atom, max_radius)
    neighbor_atoms = neighbor_search.search(atom.get_coord(), atom_cutoff)
    neighbor_atoms_list = []
    for neighbor_atom in neighbor_atoms:
        if neighbor_atom is not atom:
            neighbor_atoms_list.append(neighbor_atom)
    return neighbor_atoms_list
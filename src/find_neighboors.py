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


def get_atom_neighbors(atom, neighbor_search, cutoff):
    """
    Get the atoms within cutoff distance from an atom.

    Parameters
    ----------
    atom : Atom
        The atom object around which to search neighbors.
    neighbor_search : NeighborSearch
        A NeighborSearch instance built from all the atoms of the protein.
    cutoff : float
        The maximum distance (in Angstroms) to consider an atom a neighbor.

    Return
    ------
    list :
        A list of Atom objects within cutoff distance from atom, excluding
        the atom itself.
    """
    neighbor_atoms_list = []
    neighbor_atoms = neighbor_search.search(atom.get_coord(), cutoff)
    for neighbor_atom in neighbor_atoms:
        if neighbor_atom is not atom:
            neighbor_atoms_list.append(neighbor_atom)
    return neighbor_atoms_list
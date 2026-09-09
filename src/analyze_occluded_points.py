"""Script to analyze solvent occlusion of a parsed protein structure."""
import math

def compute_atom_surface(atom_object):
    """
    Compute point counts and surface areas (total, occluded, accessible) of one atom.
    
    Parameters
    ----------
    atom_object : Atom
        An instance of the Atom class containing information about the atom and its 
        sphere.
        
    Returns
    -------
    dict : 
        A dictionary containing the number of points, number of occluded atoms, number 
        of accessible atoms,
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
        "occluded_points": number_of_occluded_atoms,
        "accessible_point": number_of_accessible,
        "total_surface (A)": round(total_surface,2),
        "occluded_surface (A)": round(occluded_surface,2),
        "accessible_surface (A)": round(accessible_surface,2)
    }
    
    return atom_infos

def compute_residu_surface(residue_object):
    """ 
    Compute the surface area of a residue by summing the surface areas of its atoms.
    
    Parameters
    ----------
    residue_object : Residues
        An instance of the Residues class containing information about the residue and 
        its atoms.
    
    Return
    -------
    dict : 
        A dictionary containing the total number of points, number of occluded atoms, 
        number of accessible atoms,
    """
    result = {}
    atom_list = residue_object.atomlist
    for atom in atom_list:
        atom_result = compute_atom_surface(atom)
        for key, value in atom_result.items():
            result[key] = result.get(key, 0) + value
    return result


def compute_chain_surface(protein, chain_id):
    """
    Compute the surface area of a chain by summing the surface areas of its residues.
    
    Parameters
    ----------
    protein : Protein
        An instance of the Protein class containing information about the protein and 
        its residues.
    
    chain_id : str
        The identifier of the chain for which to compute the surface area.
    
    Return
    ------
    dict :  
        A dictionary containing the total number of points, number of occluded atoms, 
        number of accessible atoms,
    """

    chain_atoms = []
    result = {}
    for residue in protein.residuelist:
        if residue.chain_id == chain_id:
            chain_atoms.extend(residue.atomlist)
    for atom in chain_atoms:
        atom_result = compute_atom_surface(atom)
        for key, value in atom_result.items():
            result[key] = result.get(key, 0) + value
    return result
    
def get_chain_ids(residuelist):
    """
    Extract unique chain identifiers from a list of residues.
    
    Parameters
    ----------
    residuelist : list
        A list of Residue objects from which to extract chain identifiers.
    
    Return
    ------
    set :
        A sorted list of unique chain identifiers present in the provided list of 
        residues.
    """
    chain_ids = set()
    for residue in residuelist:
        chain_ids.add(residue.chain_id)
    return sorted(chain_ids)


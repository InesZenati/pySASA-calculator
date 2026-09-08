"Script to parse a PDB file and create the corresponding objects."
import json


from Bio.PDB import PDBParser
from loguru import logger
from object import Protein, Residues, Atom, Sphere


def load_radius_json(radius_json_file):
    """
    Load the van der Waals radius tables from a JSON file.
    
    Parameters:
    -----------
    radius_json_file : str
        Path to the JSON file containing the radius tables.
        
    Returns:
    --------
    element_radius : dict
        Dictionary mapping Van Der Waals radius for element like N, O, S, Zn.
    carbon_radius : dict
        Dictionary mapping backbone carbon atom names to their van der Waals radius.
    carbon_radius_by_residue : dict
        Dictionary mapping residue names to dictionaries of sidechain carbon atom names 
        and their van der Waals radius.
    """
    with open(radius_json_file) as file:
        radius_data = json.load(file)
    element_radius = radius_data["element_radius"]
    carbon_radius = radius_data["carbon_radius"]
    carbon_radius_by_residue = radius_data["carbon_radius_by_residue"]
    return element_radius, carbon_radius, carbon_radius_by_residue
    
def get_atom_radius_based_on_residue(element_radius, carbon_radius,
                                     carbon_radius_by_residue, residue, atomname):
    """Get the vdW radius of an atom.
    
    Parameters:
    -----------
    element_radius: dict
        Dictionary mapping Van Der Waals radius for element like N, O, S, Zn.
    carbon_radius: dict
        Dictionary mapping backbone carbon atom names to their van der Waals radius.
    carbon_radius_by_residue: dict
        Dictionary mapping residue names to dictionaries of sidechain carbon atom names
        and their van der Waals radius.    
    residue: str
        The name of the residue (e.g., "ALA", "GLY").
    atomname: str
        The name of the atom found in the PDB (e.g., "CA", "CB").
    
    Returns:
    --------
    int | None :
        The van der Waals radius of the atom if found, otherwise None.  
    """
    if atomname in element_radius :
            return element_radius[atomname]
    elif atomname[0] in element_radius:
            return element_radius[atomname[0]]
    elif atomname in carbon_radius:
            return carbon_radius[atomname]
    elif residue in carbon_radius_by_residue:
            if atomname in carbon_radius_by_residue[residue]:
                return carbon_radius_by_residue[residue][atomname]
    else : 
        logger. warning(f"Atom {atomname} in residue {residue}" 
                        "not found in radius tables.")
    return None


def parse_pdb(filename, pdbname, radius_json_file):
    """
    Parse the PDB file and create the corresponding objects
    
    Parameters:
    -----------
    element_radius: dict
        Dictionary mapping Van Der Waals radius for element like N, O, S, Zn.
    backbone_carbon_radius: dict
        Dictionary mapping backbone carbon atom names to their van der Waals radius.
    sidechain_carbon_radius_by_residue: dict
        Dictionary mapping residue names to dictionaries of sidechain carbon atom names
        and their van der Waals radius.
    pdbname: str
        The name of the PDB file (without path).
    filename: str
        The path to the PDB file.
        
    Returns:
    --------
    Protein : 
        An instance of the Protein class containing the parsed structure.
    """
    
    element_radius = radius_json_file[0]
    backbone_carbon_radius = radius_json_file[1]
    sidechain_carbon_radius_by_residue = radius_json_file[2]

    parser = PDBParser()
    structure = parser.get_structure(pdbname, filename)
    
    protein = Protein(pdbname=filename)
    # We work on the first structure that has been parsed
    model = structure[0] 
    
    for chain in model:
        logger.debug(f"Working on chain {chain.id}")
        for residue in chain:
            logger.debug(f"Working on residue {residue.get_resname()}")
            logger.debug(f" Chain {chain.id} | Residue : {residue.get_resname()}")
            my_residue =  Residues(name=residue.get_resname(),number=residue.id[1],
                                   chain_id=chain.id)
            logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()}"
                            f"created")
            for atom in residue:
                x, y, z = atom.get_coord()
                radius = get_atom_radius_based_on_residue(element_radius,
                                                          backbone_carbon_radius,
                                                          sidechain_carbon_radius_by_residue,
                                         residue.get_resname(), atom.get_name())
                my_sphere = Sphere(radius=radius, nbpoints=92)
                logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()} |" 
                            f"Atom : {atom.get_name()} | Sphere created")
                my_atom = Atom(type=atom.get_name(), atomres=my_residue, 
                               x=float(x), y=float(y), z=float(z), sphere=my_sphere)
                logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()} | "
                            f"Atom : {atom.get_name()} | Radius : {radius}") 
                
                my_residue.atomlist.append(my_atom)
            protein.residuelist.append(my_residue)

    return protein

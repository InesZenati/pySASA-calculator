import json
from Bio.PDB import PDBParser
from loguru import logger
from object import Protein, Residues, Atom, Sphere


def load_radius_json(radius_json_file):
    with open(radius_json_file) as file:
        atom_radius_value_by_residue = json.loads(file)
    return atom_radius_value_by_residue["aminoacids_atom_radius"]
    
def get_atom_radius_based_on_residue(atom_radius_value_by_residue, residue, atomname):
    residue_atom =  atom_radius_value_by_residue.get(residue, "Residue not found")
    if atomname in residue_atom:
        logger.info(f"Atom {atomname} found in residue {residue}"
                    f"with radius {residue_atom[atomname]}")
        return residue_atom[atomname]
    else:
        logger.warning(f"Atom {atomname} not found in residue {residue}")


def parse_pdb(atom_radius_value_by_residue, pdbname, filename):
    """Parse the PDB file and create the corresponding objects"""

    parser = PDBParser()
    structure = parser.get_structure(pdbname, filename)
    
    protein = Protein(pdbname=filename)
    # We work on the first structure that has been parsed
    model = structure[0] 
    
    for chain in model:
        logger.info(f"Working on chain {chain.id}")
        for residue in chain:
            logger.info(f"Working on residue {residue.get_resname()}")
            # logger.info(f" Chain {chain.id} | Residue : {residue.get_resname()}")
            my_residue =  Residues(name=residue.get_resname(),number=residue.id[1])
            logger.info(f"Chain {chain.id} | Residue : {residue.get_resname()}"
                            f"created")
            for atom in residue:
                x, y, z = atom.get_coord()
                radius = get_atom_radius_based_on_residue(atom_radius_value_by_residue,
                                                          residue.get_resname(), 
                                                          atom.get_name())           
                my_sphere = Sphere(radius=radius, nbpoints=400)
                logger.info(f"Chain {chain.id} | Residue : {residue.get_resname()} |" 
                            f"Atom : {atom.get_name()} | Sphere created")
                my_atom = Atom(type=atom.get_name(), atomres=my_residue, 
                               x=float(x), y=float(y), z=float(z), sphere=my_sphere)
                my_atom.translate_points_on_atom()
                logger.info(f"Chain {chain.id} | Residue : {residue.get_resname()} | "
                            f"Atom : {atom.get_name()} | Radius : {radius}") 
                
                my_residue.atomlist.append(my_atom)
            protein.residuelist.append(my_residue)

    return protein

if __name__ == "__main__":
    structure = parse_pdb("1CRN","data/1CRN.pdb")
    print(structure)
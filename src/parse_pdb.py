from Bio.PDB import PDBParser
from loguru import logger
from object import Protein, Residues, Atom


def parse_pdb(pdbname, filename):
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
            for atom in residue:
                logger.info(f" Chain {chain.id} | Residue : {residue.get_resname()} | Atom : {atom.get_name()}")
                x, y, z = atom.get_coord()
                my_atom = Atom(type=atom.get_name(), atomres=my_residue, x=float(x), y=float(y), z=float(z), sphere=None)
                my_residue.atomlist.append(my_atom)
            protein.residuelist.append(my_residue)

    return protein

if __name__ == "__main__":
    structure = parse_pdb("1CRN","data/1CRN.pdb")
    print(structure)
"Scrippt to parse a PDB file and create the corresponding objects."
import json
import click
import math
import time

from datetime import UTC, datetime
from tqdm import tqdm
import os
import sys

from Bio.PDB import PDBParser
from loguru import logger
from object import Protein, Residues, Atom, Sphere


def load_radius_json(radius_json_file):
    """Load the vdW radius tables from a JSON file."""
    with open(radius_json_file) as file:
        radius_data = json.load(file)
    element_radius = radius_data["element_radius"]
    backbone_carbon_radius = radius_data["backbone_carbon_radius"]
    sidechain_carbon_radius_by_residue = radius_data["sidechain_carbon_radius_by_residue"]
    return element_radius, backbone_carbon_radius, sidechain_carbon_radius_by_residue
    
def get_atom_radius_based_on_residue(element_radius, backbone_carbon_radius, sidechain_carbon_radius_by_residue,
                     residue, atomname):
    """Get the vdW radius of an atom."""
    if atomname == "OXT":
        return element_radius["O"]

    if atomname in backbone_carbon_radius:
        return backbone_carbon_radius[atomname]

    is_carbon = atomname[0] == "C"
    if is_carbon:
        residue_carbons = sidechain_carbon_radius_by_residue.get(residue, {})
        if atomname in residue_carbons:
            return residue_carbons[atomname]
        logger.warning(f"Carbon {atomname} of residue {residue} not in table, using default 2.0")
        return 2.0

    element = atomname[0]
    if element in element_radius:
        return element_radius[element]

    logger.error(f"Atom {atomname} of residue {residue} unknown, no matching element found")
    return None


def parse_pdb(element_radius, backbone_carbon_radius, sidechain_carbon_radius_by_residue,
              pdbname, filename):
    """Parse the PDB file and create the corresponding objects"""

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
            my_residue =  Residues(name=residue.get_resname(),number=residue.id[1])
            logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()}"
                            f"created")
            for atom in residue:
                x, y, z = atom.get_coord()
                radius = get_atom_radius_based_on_residue(element_radius, backbone_carbon_radius,
                                         sidechain_carbon_radius_by_residue,
                                         residue.get_resname(), atom.get_name())           
                my_sphere = Sphere(radius=radius, nbpoints=92)
                logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()} |" 
                            f"Atom : {atom.get_name()} | Sphere created")
                my_atom = Atom(type=atom.get_name(), atomres=my_residue, 
                               x=float(x), y=float(y), z=float(z), sphere=my_sphere)
                my_atom.translate_points_on_atom()
                logger.debug(f"Chain {chain.id} | Residue : {residue.get_resname()} | "
                            f"Atom : {atom.get_name()} | Radius : {radius}") 
                
                my_residue.atomlist.append(my_atom)
            protein.residuelist.append(my_residue)

    return protein

if __name__ == "__main__":
    timestamp = datetime.now(UTC).strftime("%Y%m%d_%H%M%S")
    os.makedirs("logs", exist_ok=True)
    logger_format = (
        "{time:YYYY-MM-DD HH:mm:ss} "
        "| <level>{level:<8}</level> "
        "| <level>{message}</level>"
    )
    logger.remove()
    logger.add(sys.stdout, format=logger_format, level="INFO")
    logger.add(
        f"logs/normalize_simulation_time{timestamp}.log",
        level="INFO",
        format=logger_format,
    )
    element_radius, backbone_carbon_radius, sidechain_carbon_radius_by_residue = load_radius_json("data/radius.json")
    structure = parse_pdb(element_radius, backbone_carbon_radius, sidechain_carbon_radius_by_residue,
                          "1CRN", "data/1CRN.pdb")

    all_atoms = structure.get_all_atom()
    logger.info(f"Total number of atoms: {len(all_atoms)}")

    start_time = time.time()
    for atom in tqdm(all_atoms, desc="Detecting occluded points"):
        atom.detect_occluded_point(all_atoms)
    elapsed_time = time.time() - start_time
    
total_surface = 0
total_occluded_surface = 0

for atom in all_atoms:
    npoints = len(atom.sphere.pointlist)
    n_occluded = sum(atom.sphere.occluded_points)

    surface = 4 * math.pi * atom.sphere.radius**2

    total_surface += surface
    total_occluded_surface += (n_occluded / npoints) * surface

accessible_surface = total_surface - total_occluded_surface

logger.info(f"Total surface: {total_surface:.2f} Å²")
logger.info(f"Occluded surface: {total_occluded_surface:.2f} Å²")
logger.info(f"Accessible surface: {accessible_surface:.2f} Å²")
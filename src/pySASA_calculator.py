"""Script to analyze solvent occlusion of a parsed protein structure."""
import os
import sys
import time
from pathlib import Path
from datetime import UTC, datetime

import click
from loguru import logger
from tqdm import tqdm

from parse_pdb import load_radius_json, parse_pdb
from write_surface_results import write_surface_results

@click.command()
@click.option("--pdb-name", 
              type=str,
              help="Name of the protein (PDB id).")                      
@click.option("--pdb-file", 
              type=click.Path(file_okay=True, path_type=Path), 
              help="Path to the PDB file.")
@click.option("--radius-json",
              type=click.Path(exists=True),
              help="Path to the JSON file with van der Waals radius tables.")
def analyze_protein_sasa(pdb_name, pdb_file, radius_json):
    radius_table = load_radius_json(radius_json)
    protein = parse_pdb(pdb_file, pdb_name, radius_table)

    all_atoms = protein.get_all_atom()
    logger.info(f"Total number of atoms: {len(all_atoms)}")

    start_time = time.time()
    for atom in tqdm(all_atoms, desc="Detecting occluded points"):
        atom.detect_occluded_point(all_atoms)
    elapsed_time = time.time() - start_time
    logger.success(f"Occlusion detection took {elapsed_time:.2f} seconds")

    write_surface_results(protein, "results")
    logger.success("Surface results written to results/")
    
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
    logger.add(f"logs/analyze_occlusion_{timestamp}.log", format=logger_format, 
               level="INFO")
    
    analyze_protein_sasa()
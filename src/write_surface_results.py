import csv
import os
from analyze_occluded_points import compute_atom_surface, compute_residu_surface\
                                    ,get_chain_ids, compute_chain_surface


def write_tsv(rows, column_names, output_file_path):
    """
    Write a list of rows to a TSV file.

    Parameters
    ----------
    rows : list
        A list of dictionaries, one per row, with keys matching column_names.
    column_names : list
        The column names, in the order they should appear in the file.
    output_file_path : str
        Path to the output TSV file.
    """
    with open(output_file_path, "w", newline="") as tsv_file:
        tsv_writer = csv.DictWriter(tsv_file, column_names, delimiter="\t")
        tsv_writer.writeheader()
        for row in rows:
            tsv_writer.writerow(row)


def get_atom_surface_rows(protein):
    """
    Compute the surface area of every atom of the protein.

    Parameters
    ----------
    protein : Protein
        An instance of the Protein class containing the parsed structure.

    Return
    ------
    list :
        A list of dictionaries, one per atom.
    """
    atom_rows = []
    for atom in protein.get_all_atom():
        atom_row = compute_atom_surface(atom)
        atom_row["chain_id"] = atom.atomres.chain_id
        atom_row["residue_name"] = atom.atomres.name
        atom_row["residue_number"] = atom.atomres.number
        atom_row["atom_type"] = atom.type
        atom_rows.append(atom_row)
    return atom_rows


def get_residue_surface_rows(protein):
    """
    Compute the surface area of every residue of the protein.

    Parameters
    ----------
    protein : Protein
        An instance of the Protein class containing the parsed structure.

    Return
    ------
    list :
        A list of dictionaries, one per residue.
    """
    residue_rows = []
    for residue in protein.residuelist:
        residue_row = compute_residu_surface(residue)
        residue_row["chain_id"] = residue.chain_id
        residue_row["residue_name"] = residue.name
        residue_row["residue_number"] = residue.number
        residue_rows.append(residue_row)
    return residue_rows


def get_chain_surface_rows(protein):
    """
    Compute the surface area of every chain of the protein.

    Parameters
    ----------
    protein : Protein
        An instance of the Protein class containing the parsed structure.

    Return
    ------
    list :
        A list of dictionaries, one per chain.
    """
    chain_rows = []
    for chain_id in get_chain_ids(protein.residuelist):
        chain_row = compute_chain_surface(protein, chain_id)
        chain_row["chain_id"] = chain_id
        chain_rows.append(chain_row)
    return chain_rows


def get_protein_surface_rows(protein):
    """
    Compute the total surface area of the protein.

    Parameters
    ----------
    protein : Protein
        An instance of the Protein class containing the parsed structure.

    Return
    ------
    list :
        A single-item list with the total surface area of the protein.
    """
    protein_row = {}
    for atom in protein.get_all_atom():
        atom_row = compute_atom_surface(atom)
        for key, value in atom_row.items():
            protein_row[key] = protein_row.get(key, 0) + value
    protein_row["pdb_name"] = protein.pdbname
    return [protein_row]


def write_all_surface_results_in_tsv(protein, output_dir):

    surface_field_names = [
        "number_of_points", "number_of_occluded_atoms", "number_of_accessible",
        "total_surface", "occluded_surface", "accessible_surface",
    ]
    os.makedirs(output_dir, exist_ok=True)

    atom_field_names = ["chain_id", "residue_name", "residue_number", "atom_type"] + \
        surface_field_names
    write_tsv(get_atom_surface_rows(protein), atom_field_names,
              os.path.join(output_dir, "atom_surface.tsv"))
    residue_field_names = ["chain_id", "residue_name", "residue_number"] + \
        surface_field_names
    write_tsv(get_residue_surface_rows(protein), residue_field_names,
              os.path.join(output_dir, "residue_surface.tsv"))

    chain_field_names = ["chain_id"] + surface_field_names
    write_tsv(get_chain_surface_rows(protein), chain_field_names,
              os.path.join(output_dir, "chain_surface.tsv"))

    protein_field_names = ["pdb_name"] + surface_field_names
    write_tsv(get_protein_surface_rows(protein), protein_field_names,
              os.path.join(output_dir, "protein_surface.tsv"))





# pySASA-calculator

## Introduction

**pySASA-calculator** computes the Solvent Accessible Surface Area (SASA) of a protein from a PDB structure. For each non-hydrogen atom, a sphere of points is placed at its van der Waals radius; points occluded by neighboring atoms are detected to derive the accessible surface at the atom, residue, chain, and whole-protein level.

## Setup environment

We use [uv](https://docs.astral.sh/uv/getting-started/installation/) to manage dependencies and the project environment.

Clone the GitHub repository:

```sh
git clone https://github.com/InesZenati/pySASA-calculator.git
cd pySASA-calculator
```

Sync dependencies:

```sh
uv sync
```

## Van der Waals radius file

To build the sphere of each atom, the program needs a JSON file describing the standard van der Waals radius of each atom type. This file is not included in the repository and must be provided by the user (a `data/radius.json` example can be created following the structure below).

| Key | Description |
|---|---|
| `element_radius` | Radius per element (N, O, S, Zn, etc.) |
| `carbon_radius` | Radius of backbone carbon atoms |
| `carbon_radius_by_residue` | Radius of side-chain carbon atoms, per residue |

## Usage

Run the SASA calculation on a PDB file using:

```sh
uv run src/pySASA_calculator.py \
    --pdb-name 168L \
    --pdb-file data/168L.pdb \
    --radius-json data/radius.json \
    --point_number 100
```

| Option | Description |
|---|---|
| `--pdb-name` | Name (PDB code) of the protein |
| `--pdb-file` | Path to the PDB file to analyze |
| `--radius-json` | Path to the van der Waals radius JSON file |
| `--point_number` | Number of points used to sample each atom's sphere |

Two example structures are provided in `data/`: `1CRN.pdb` and `1MH1.pdb`.

## Project structure

```
pySASA-calculator/
├── data/                          # Example PDB structures (1CRN.pdb, 1MH1.pdb)
├── src/
│   ├── pySASA_calculator.py       # Entry point: CLI, orchestrates the full pipeline, sets up logging
│   ├── parse_pdb.py               # Loads the van der Waals radius JSON and parses the PDB file into Protein/Residues/Atom/Sphere objects
│   ├── protein_structure.py       # Data model: Protein, Residues, Atom and Sphere classes (sphere points generated with a Fibonacci sphere algorithm)
│   ├── find_neighboors.py         # Builds a Biopython NeighborSearch index and returns, for each atom, the neighbors that could occlude it
│   ├── analyze_occluded_points.py # Computes surface areas (total/occluded/accessible) at the atom, residue and chain level
│   ├── write_surface_results.py   # Assembles the per-atom/residue/chain/protein rows and writes them to the TSV files in results/
│   └── visualize_sphere.py        # Standalone script to plot a single atom's sphere of points in 3D (matplotlib), for debugging/illustration
├── pyproject.toml                 # Project metadata and dependencies (managed with uv)
└── uv.lock                        # Locked dependency versions
```

## Output

Results are written to `results/<pdb-name>/` as four TSV files, one per level of granularity:

| File | Content |
|---|---|
| `atom_surface.tsv` | Surface area per atom |
| `residue_surface.tsv` | Surface area per residue |
| `chain_surface.tsv` | Surface area per chain |
| `protein_surface.tsv` | Total surface area of the protein |

Each file reports, per entry: `number_of_points`, `occluded_points`, `accessible_point`, `total_surface (A)`, `occluded_surface (A)`, `accessible_surface (A)`, and `relative_surface` [WIP].

## Notes

A log file is created after every execution in `logs/analyze_occlusion_<timestamp>.log`.

By default, to gain running time, the logger is set to the `INFO` level, so only progress messages are recorded. To get detailed step-by-step information in the log (e.g. per-point sphere coordinates, per-atom radius lookups), change the level passed to `logger.add(...)` in `src/pySASA_calculator.py` from `"INFO"` to `"DEBUG"`.
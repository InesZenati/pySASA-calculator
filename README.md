# pySASA-calculator
## Introduction
This project compute the Solvent Accessible Surface Area (SASA) of a protein. 

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

To create the sphere of each atom we use a JSON file that contains the standard Van der Waals radius of each atoms

| Clé | Description |
|---|---|
| `element_radius` | Rayon par élément (N, O, S, Zn, etc.) |
| `carbon_radius` | Rayon des carbones du squelette |
| `carbon_radius_by_residue` | Rayon des carbones de chaîne latérale, par résidu |

## Usage
 
Run the SASA calculation on a PDB file using :
 
```sh
uv run analyze_occlusion.py \
    --pdb-name 1ABC \
    --pdb-file data/1abc.pdb \
    --radius-json data/vdw_radius.json
```

| Option | Description |
|---|---|
| `--pdb-name` | Nom (code PDB) de la protéine |
| `--pdb-file` | Chemin vers le fichier PDB à analyser |
| `--radius-json` | Chemin vers le fichier JSON des rayons de van der Waals |

## Notes
A log file will be created after every execution in logs/analyze_occlusion_<timestamp>.log`
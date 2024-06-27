# ProFlex
Explore ProFlex pipeline for interface diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE. Optimize protein interactions dynamically.

For now the project includes scripts for preparation of inputs and outputs into the RFdiffusion-FastRelax-PELE pipeline. Soon it will be automatized.

## Installation (WIP)
### Environment manual setup
```bash
$ git clone git@github.com:ibalafkir/ProFlex.git
$ conda create --name proflex python=3.11
$ conda activate proflex
$ if [ -f requirements.txt ]; then pip install -r ./ProFlex/requirements.txt; fi
```

### Automated setup
For developers:
```bash
$ python setup.py develop
```
For users:
```bash
$ python setup.py install
```

## Usage
### Manual system preparation
The system should be prepared before running the protocol. In this context, this preparation includes: removing of water/metal/ligand atoms of no relevance for antibody-antigen interactions; correction of missing atoms, side chains and loops; calculation of protonation states at pH 7.4; and optimization of the hydrogen-bonding network. The protein preparation wizard tool implemented in the Schrödinger suite is generally used to do so. 

Assuming ```$SCHRODINGER``` contains the path to the protein preparation wizard.

```bash
$ $SCHRODINGER input.pdb output.pdb -rehtreat -disulfides -fillloops -fillsidechains -propka_pH 7.4 -minimize_adj_h -f OPLS_2005
```

Warning: each system has its own peculiarities and further adjustments might be needed

### Code organization

#### Main Scripts
- **proflex/contigs.py**
  - Gets the contigs for RFdiffusion

- **proflex/fix_diffusion_models.py**
  - Fixes diffusion modals for correct PDB format (chains ID, END, TER...)

- **main/interface_analyzer.py**
  - Gets interface metrics (Ab-Ag interface, Ab-Ab interface) for PELE control files

- **proflex/merge_antibody_chains.py**
  - Merges antibody chains under a single chain name for DockQ calculations

- **proflex/pdb_preprocessor.py**
  - Preproccesses a PDB for entering the pipeline. Generally for fixing and renumbering antibody insertion codes

- **proflex/pele_postprocess.py**
  - Changes certain residue names of PELE outputs

- **proflex/remove_sidechains.py**
  - Removes side chains of a given PDB (for repacking analysis)

- **proflex/repack_assessment_preprocess.py**
  - Merges all residues under the same chain name for repacking analysis using AttnPacker scripts

- **proflex/resnum.py**
  - Calculates the number of residue names

#### Modules
- **proflex/interface**
  - Classes related to interface analysis (used in contigs calculations, PELE metrics...)

- **proflex/pdb**
  - Classes related to pdb processing (as dataframes or pdb files)

- **proflex/rfdiffusion**
  - Classes related to RFdiffusion contigs generation and outputs postprocessing

#### Future improvements (WIP)

The code will be improved with the next features:
- Contigs
- Interface and neighbourhood amplification (use of biopython and increase the selected regions)
- More preprocessing/postprocessing adjustments if needed
- Tests
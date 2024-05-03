# ProFlex
Explore ProFlex protocol for interface diffusion and Monte Carlo rotations/translations using RFdiffusion and PELE. Optimize protein interactions dynamically.

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
In general, Schrodinger with Prime license through its protein preparation wizard program is needed for the system preparation steps. We'll assume ```$SCHRODINGER``` contains the path to the protein preparation wizard.

Generally, for filling missing loops and side chains, creating disulfide bonds between proximal sulfurs adding and optimizing hydrogens and minimizing in the OPLS force field run:
```bash
$ $SCHRODINGER input.pdb output.pdb -rehtreat -disulfides -fillloops -fillsidechains -propka_pH 7.0 -minimize_adj_h -f OPLS_2005
```
Although bear in mind that each system has its own peculiarities and more preparation actions might need to be run.

To keep the desired chains of a system and remove insertion codes, hetatms and unnecessary PDB lines run ./ProFlex/proflex/preprocessor.py.
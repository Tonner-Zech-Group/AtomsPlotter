# Atoms_plotter
Package for plotting ASE atoms-objects using matplotlib. Can be used to create vector graphics, or used interactively in a jupyter-notebook.

## Installation
clone repository into Atoms_plotter
install with 
```
cd Atoms_plotter
python3 -m pip install -e . 
```
Ideally with a development install (-e) so changes are easily updated.

## Usage
* use installed entry point `plot_atoms`. See `plot_atoms --help` for options
### In a python script
Initialize the class; set options either when initialising or modify class attributes.
```python
from atoms_plotter.atoms_plotter import atoms_plotter
from ase.io import read
from sys import argv
atoms=read(argv[1])
ap=atoms_plotter(atoms,show=True)
ap.plot()
```
### In a jupyter-notebook

```python
import ase
import ase.io
import matplotlib.pyplot
atoms=ase.io.read('CONTCAR')
%matplotlib widget #for interactive use
from atoms_plotter.atoms_plotter import atoms_plotter

plotter=atoms_plotter(atoms=atoms, constraints=True)
plotter.dimension='3D'
plotter.atoms=atoms
plotter.plot()
```
#### some options
```python
plotter.scale=70
plotter.bondlinewidth=2
plotter.colorbonds=True
plotter.constraints=False #draw constraints as x
plotter.draw_outline=True
plotter.repeat=[2,2]
plotter.use_bondorders=True #double/triple bonds
plotter.dimension='2D' #2D projection
plotter.plot()

#plot lewis-like structures. bonds as simple black lines, no atoms.

plotter.dimension='2D'
plotter.lewis=True
plotter.scale=0
plotter.atoms.rotate(90, 'z',rotate_cell=True)
plotter.view=2
plotter.plot()
```

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
from AtomsPlotter.atoms_plotter import atoms_plotter
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
from AtomsPlotter.atoms_plotter import atoms_plotter

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

## Changelog

### 0.3.3

* Fix `plot_atoms_2D` / `plot_atoms_3D` colour assignment: elements
  missing from `color_dict` (e.g. `Au`, `Pt`) used to inherit the
  last-iterated atom's jmol colour because `self.COLORS` was rebuilt
  inside the per-atom loop with the outer-loop atomic number as the
  fallback. Now built once with a per-atom jmol fallback.
* Add `Au` (gold) to the default `color_dict`.
* Promote `plot_atom_colorscaling` to a real constructor keyword
  (`plot_atom_colorscaling=False` by default). Previously it was
  referenced from both render paths but the `__init__` assignment was
  commented out, so any code that hit the `color_dict=None` branch
  raised `AttributeError`.
* Guard the z-fraction colour scaling against `z_norm == 0` so flat
  (2-D) structures no longer divide by zero.
* Tidy four `ruff E702` (semicolon-chained statements) in `check_pbc`.

*(0.3.2 was prepared but never published — its scope is included
in 0.3.3 above.)*

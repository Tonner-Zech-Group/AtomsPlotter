from AtomsPlotter.atoms_plotter import atoms_plotter
from ase.io import read
from sys import argv
atoms=read(argv[1])
ap=atoms_plotter(atoms,show=True)
ap.plot()

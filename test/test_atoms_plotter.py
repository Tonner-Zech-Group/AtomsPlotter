import sys

import matplotlib

matplotlib.use('Agg')

import numpy as np
import pytest
from ase.build import bulk, molecule

from AtomsPlotter import atoms_plotter
from AtomsPlotter.main import main


@pytest.fixture
def in_tmp(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_molecule_2d_saves(in_tmp):
    ap = atoms_plotter(molecule('C6H6'), show=False, name='benzene')
    ap.plot()
    assert (in_tmp / 'benzene.svg').is_file()


def test_molecule_arbitrary_placement_is_framed(in_tmp):
    mol = molecule('C6H6')
    mol.translate([100.0, -50.0, 7.0])
    before = mol.positions.copy()
    ap = atoms_plotter(mol, show=False, name='displaced')
    ap.plot()
    assert (in_tmp / 'displaced.svg').is_file()
    # plotting must not mutate the user's atoms
    assert np.allclose(mol.positions, before)
    assert mol.cell.rank == 0


def test_periodic_2d_with_repeat(in_tmp):
    si = bulk('Si', 'diamond', a=5.43, cubic=True)
    ap = atoms_plotter(si, show=False, name='si_2x2', repeat=[2, 2],
                       show_unit_cell=True)
    ap.plot()
    assert (in_tmp / 'si_2x2.svg').is_file()
    # repeat works on a copy
    assert len(ap.atoms) == len(si)


def test_figure_scales_with_cell(in_tmp):
    si = bulk('Si', 'diamond', a=5.43, cubic=True)
    # both sizes chosen large enough to stay clear of the figsize clamps
    ap_small = atoms_plotter(si, show=False, name='small', repeat=[2, 2])
    ap_small.plot()
    ap_large = atoms_plotter(si, show=False, name='large', repeat=[4, 4])
    ap_large.plot()
    w_small = ap_small.fig.get_size_inches()[0]
    w_large = ap_large.fig.get_size_inches()[0]
    assert w_large == pytest.approx(2 * w_small, rel=0.05)


def test_nonorthogonal_cell_frame_covers_all_corners(in_tmp):
    import ase
    a = 2.46
    graphene = ase.Atoms('C2',
                         positions=[[0, 0, 5], [a / np.sqrt(3), 0, 5]],
                         cell=[[a, 0, 0],
                               [a / 2, a * np.sqrt(3) / 2, 0],
                               [0, 0, 10]],
                         pbc=True)
    ap = atoms_plotter(graphene, show=False, name='graphene', repeat=[3, 3])
    ap.plot()
    cell = graphene.cell.array * 3
    # x extent must include the sheared corner a1+a2, not just a1
    assert ap.frame[0][1] >= cell[0][0] + cell[1][0] - 1e-8


def test_color_fallback_per_atom(in_tmp):
    atoms = molecule('H2O')
    atoms.symbols[2] = 'Pt'  # not in color_dict -> jmol fallback
    ap = atoms_plotter(atoms, show=False, name='fallback')
    ap.plot()
    from ase.data import colors
    assert tuple(ap.COLORS[0]) == ap.color_dict['O']
    pt_index = list(ap.atoms.symbols).index('Pt')
    assert np.allclose(ap.COLORS[pt_index], colors.jmol_colors[78])


def test_legacy_scale_path(in_tmp):
    ap = atoms_plotter(molecule('C6H6'), show=False, name='legacy', scale=100)
    ap.plot()
    assert ap.sizes[0] == pytest.approx(100 * 0.76)  # C covalent radius


def test_lewis_mode(in_tmp):
    ap = atoms_plotter(molecule('C6H6'), show=False, name='lewis', lewis=True)
    ap.plot()
    assert (in_tmp / 'lewis.svg').is_file()


def test_3d(in_tmp):
    ap = atoms_plotter(molecule('CH3CH2OH'), show=False, name='ethanol3d',
                       dimension='3D', elev=10, azim=30)
    ap.plot()
    assert (in_tmp / 'ethanol3d.svg').is_file()


def test_linear_molecule_auto_orients():
    ap = atoms_plotter(molecule('C2H2'))
    ap._orient_principal_axes()
    spans = ap.atoms.positions.max(axis=0) - ap.atoms.positions.min(axis=0)
    assert spans[0] == pytest.approx(spans.max())
    assert spans[2] == pytest.approx(0, abs=1e-8)


def test_cli(in_tmp, monkeypatch):
    from ase.io import write
    # '@' in the name checks that the CLI does not treat it as file@index
    write('mol@test.xyz', molecule('NH3'))
    monkeypatch.setattr(sys, 'argv', ['plot_atoms', 'mol@test.xyz'])
    main()
    assert (in_tmp / 'mol@test.svg').is_file()


def test_cli_repeat_and_options(in_tmp, monkeypatch):
    from ase.io import write
    write('POSCAR', bulk('Si', 'diamond', a=5.43, cubic=True), format='vasp')
    monkeypatch.setattr(sys, 'argv', [
        'plot_atoms', 'POSCAR', '-r', '2', '2', '1', '-u', '-N', 'si_cli',
        '-D', '3D', '-f', 'png'])
    main()
    assert (in_tmp / 'si_cli.png').is_file()

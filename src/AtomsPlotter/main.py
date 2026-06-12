import os

from ase.io import read

from .atoms_plotter import atoms_plotter


def main():
    import argparse

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n'):
            return False
        else:
            return v
    parser = argparse.ArgumentParser(
        description='Plotter using utility of ASE, in matplotlib')
    parser.add_argument(
        'atoms', type=str, help='atoms object, readable by ASE .poscar, .vasp, .xyz, ...')
    parser.add_argument('-r', '--repeat', default=[
                        1, 1, 1], help='repetition of the unit cell in [x,y,z] direction', type=int, nargs=3)
    parser.add_argument('-b', '--bondlinewidth',
                        help='line width of the bonds (multiplier of the default stick width)', default=1, type=float)
    parser.add_argument(
        '-B', '--bonds', help='plot the atoms but not the bonds', default=True, type=str2bool)
    parser.add_argument(
        '-u', '--unitcell', help='plot unit cell. Buggy if atoms are rotated', action='store_true')
    parser.add_argument(
        '-S', '--show', help='decide if you only want to save or show the figure', action='store_true')
    parser.add_argument('-N', '--name', default=None,
                        help='name for the saved plot if show == False', type=str)
    parser.add_argument('-f', '--format', default='svg',
                        help='output format (svg, png, pdf, ...)', type=str)
    parser.add_argument('-s', '--scale', default=None, type=float,
                        help='legacy fixed atom size (points^2 per Angstrom covalent radius); '
                             'by default sizes are derived from the structure')
    parser.add_argument('--atom_scale', default=0.5, type=float,
                        help='fraction of the covalent radius the atoms are drawn with')
    parser.add_argument('--inches_per_angstrom', default=0.4, type=float,
                        help='figure size per Angstrom of cell/structure extent')
    parser.add_argument(
        '-c', '--colorbonds', help='use the colors of the atoms for the bonds; or give value between 1 and 0; or False', type=str2bool, default=True)
    parser.add_argument('-O', '--outline', help='draw outline',
                        type=str2bool, default=True)
    parser.add_argument(
        '-j', '--jmol', help='use the jmol colors instead of default-dict', action='store_true')
    parser.add_argument('-D', '--dimension', help='plotting in 2 or 3D',
                        choices=['2D', '3D'], default='2D')
    parser.add_argument(
        '--elev', help='elevation angle of the 3D view', default=90, type=float)
    parser.add_argument(
        '--azim', help='azimuth angle of the 3D view', default=90, type=float)
    parser.add_argument('--unit_cell_linestyle', default='dashed', type=str)
    parser.add_argument('-C', '--cutoff_radius', default=1,
                        help='cutoff radius for drawing the bonds', type=float)
    parser.add_argument(
        '--no_bondorders', help='call if you don\'t want double or triple bonds', action='store_false')
    parser.add_argument('--triple_bonds', default=1.34,
                        help='distance below which bonds are read as triple bonds', type=float)
    parser.add_argument('--double_bonds', default=1.4,
                        help='distance below which bonds are read as double bonds', type=float)
    parser.add_argument('--double_bond_offset', default=0.1,
                        help='offset b will be used for plotting double bonds [b,b,0]', type=float)
    parser.add_argument('--triple_bond_offset', default=0.07,
                        help='offset t will be used for plotting triple bonds [t,t,0]', type=float)
    parser.add_argument(
        '--bond_gradient', help='use gradient for coloring the bonds', type=str2bool, default=True)
    parser.add_argument(
        '--constraints', help='draw constraints', type=str2bool, default=True)
    parser.add_argument(
        '--auto_orient', type=str2bool, default=None,
        help='align principal axes with x,y,z before plotting; '
             'default: on for molecules, off for periodic structures')
    parser.add_argument(
        '--rot', help='rotation of the atoms in [x,y,z] direction', type=float, nargs=3)
    parser.add_argument(
        '-v', '--view', help="view direction for z-order, use 0,1,2 for x,y,z", default=2, type=int)
    parser.add_argument(
        '-l', '--plot_lewis', help='plot lewis structure instead of colored', action='store_true')
    a = parser.parse_args()
    # an existing path containing '@' is a filename, not ASE's file@index
    if os.path.exists(a.atoms) and '@' in a.atoms:
        atoms = read(a.atoms, do_not_split_by_at_sign=True)
    else:
        atoms = read(a.atoms)
    plotter = atoms_plotter(atoms=atoms,
                            show=a.show,
                            lewis=a.plot_lewis,
                            show_unit_cell=a.unitcell,
                            unit_cell_linestyle=a.unit_cell_linestyle,
                            repeat=a.repeat,
                            bondlinewidth=a.bondlinewidth,
                            bond_cutoff=a.cutoff_radius,
                            show_bonds=a.bonds,
                            colorbonds=a.colorbonds,
                            draw_outline=a.outline,
                            constraints=a.constraints,
                            scale=a.scale,
                            atom_scale=a.atom_scale,
                            inches_per_angstrom=a.inches_per_angstrom,
                            use_bondorders=a.no_bondorders,
                            triple_bond=a.triple_bonds,
                            double_bond=a.double_bonds,
                            double_bond_offset=[
                                a.double_bond_offset, a.double_bond_offset, 0],
                            triple_bond_offset=[
                                a.triple_bond_offset, a.triple_bond_offset, 0],
                            bond_gradient=a.bond_gradient,
                            dimension=a.dimension,
                            view=a.view,
                            auto_orient=a.auto_orient,
                            azim=a.azim,
                            elev=a.elev,
                            format=a.format,
                            )
    if a.rot:
        plotter.atoms.wrap(pretty_translation=True)
        plotter.atoms.rotate(a.rot[0], 'x', center='COM', rotate_cell=True)
        plotter.atoms.rotate(a.rot[1], 'y', center='COM', rotate_cell=True)
        plotter.atoms.rotate(a.rot[2], 'z', center='COM', rotate_cell=True)
    if a.name is None:
        plotter.name = os.path.splitext(a.atoms)[0]
    else:
        plotter.name = a.name
    if a.jmol is True:
        plotter.color_dict = None
    plotter.plot()


if __name__ == "__main__":
    main()

""" Convert VASP input files to QE """

from math import ceil

INCAR = 'INCAR'
# TAGS = ['EDIFF', 'EDIFFG', 'ENCUT', 'IBRION', 'ISIF', 'ISMEAR', 'ISPIN',
#         'ISYM', 'LDAU', 'LDAUTYPE', 'LDAUU', 'LWAVE', 'MAGMOM', 'METAGGA',
#         'NELM', 'NSW', 'SIGMA']
# TYPES = [float, float, int, int, int, int, int,
#          int, bool, list, list, bool, list, str,
#          int, int, float]
# UNITS = ['eV', 'eV/A', 'eV', None, None, None, None,
#          None, None, None, None, None, None, None,
#          None, None, 'eV']

TAGS = ['EDIFF', 'EDIFFG', 'ENCUT', 'IBRION', 'ISIF', 'ISMEAR', 'ISPIN',
        'ISYM',
        'NELM', 'NSW', 'SIGMA']
TYPES = [float, float, int, int, int, int, int,
         int,
         int, int, float]
UNITS = ['eV', 'eV/A', 'eV', None, None, None, None,
         None,
         None, None, 'eV']
POSCAR = 'CONTCAR'
KPOINTS = 'KPOINTS'
PAW_PATH = '/u/rwexler/PAWs'


class Incar:
    """
Class for representing VASP INCAR files
    """

    def __init__(self, incar, tags, types, units):
        self.incar = incar
        self.tags = tags
        self.types = types
        self.units = units
        self.values = {}

    def read_incar(self):
        """
Read INCAR tags
        """
        for tag, typ, unit in zip(self.tags, self.types, self.units):
            with open(self.incar) as file:
                for line in file:
                    if tag in line:
                        break
                value = line.split('=')[-1].replace('\n', '').strip()
                if typ == float:
                    if unit == 'eV':
                        value = float(value) / 13.605693122994
                    elif unit == 'eV/A':
                        value = abs(
                            float(value)) / 13.605693122994 * 0.529177210903
                elif typ == int:
                    if unit == 'eV':
                        value = ceil(float(value) / 13.605693122994)
                    else:
                        value = int(value)
                elif typ == bool:
                    value = bool(value)
                self.values[tag] = value

    def parse_magmom(self):
        """
Parse MAGMOM
        """
        self.values['MAGMOM'] = [x.split('*') for x in
                                 self.values['MAGMOM'].split()]

    def clean_ldauu_values(self):
        """
Clean LDAUU values
        """
        self.values['LDAUU'] = [float(x) for x in self.values['LDAUU'].split()]


class Poscar:
    """
Class for representing VASP POSCAR and CONTCAR files
    """

    def __init__(self, poscar):
        self.poscar = poscar
        self.lattice_vectors = []
        self.symbols = []
        self.numbers = []
        self.number_of_atoms = None
        self.coordinates = []
        self.values = {}

    def read_poscar(self):
        """
Read POSCAR/CONTCAR
        """
        with open(self.poscar) as file:
            next(file)
            next(file)

            # get lattice vectors
            for i in range(3):
                for line in file:
                    self.lattice_vectors.append(line.split())
                    break

            # get element symbols and number of each element
            for i in range(2):
                for line in file:
                    if i == 0:
                        self.symbols = line.split()
                        break
                    else:
                        self.numbers = [int(x) for x in line.split()]
                        break
            self.number_of_atoms = sum(self.numbers)

            # get coordinates
            next(file)
            for i in range(self.number_of_atoms):
                for line in file:
                    self.coordinates.append([float(x) for x in line.split()])
                    break
        self.values['lattice_vectors'] = self.lattice_vectors
        self.values['symbols'] = self.symbols
        self.values['number'] = self.numbers
        self.values['number_of_atoms'] = self.number_of_atoms
        self.values['coordinates'] = self.coordinates


class Kpoints:
    """
Class for representing VASP KPOINTS files
    """

    def __init__(self, kpoints):
        self.kpoints = kpoints
        self.center = None
        self.grid = []
        self.values = {}

    def read_kpoints(self):
        """
Read KPOINTS
        """
        with open(self.kpoints) as file:
            next(file)
            next(file)
            for line in file:
                self.center = line.replace('\n', '')
                break
            for line in file:
                self.grid = [int(x) for x in line.split()]
        self.values['center'] = self.center
        self.values['grid'] = self.grid


class Qe:
    """
Class for representing QE input files
    """

    def __init__(self, values, paw_path):
        self.values = values
        self.paw_path = paw_path
        self.control = None
        self.qe = None

    def write_qe(self):
        """
Write QE input file
        """
        control_tags = ['calculation', 'nstep', 'forc_conv_thr', 'disk_io',
                        'pseudo_dir']
        system_tags = ['ibrav', 'nat', 'ntyp', 'starting_magnetization',
                       'ecutwfc', 'ecutrho', 'nosym', 'occupations',
                       'degauss', 'smearing', 'nspin', 'input_dft',
                       'lda_plus_u', 'lda_plus_u_kind', 'Hubbard_U']
        electrons_tags = ['electron_maxstep', 'conv_thr']
        # for index, tag in enumerate(control_tags):
        #     if index == 0:
        #         self.qe = '&CONTROL\n'
        #     if tag == 'calculation':
        #         if self.values['ISIF'] == 3:
        #             self.qe += '  {} = {},\n'.format(tag, '"vc-relax"')
        #         elif self.values['ISIF'] < 3:
        #             self.qe += '  {} = {},\n'.format(tag, '"relax"')
        #         else:
        #             raise ValueError('ISIF is greater than 3')
        #     if tag == 'nstep':
        #         self.qe += '  {} = {},\n'.format(tag, self.values['NSW'])
        #     if tag == 'forc_conv_thr':
        #         self.qe += '  {} = {},\n'.format(tag, self.values['EDIFFG'])
        #     if tag == 'disk_io':
        #         if self.values['LWAVE']:
        #             pass
        #         else:
        #             self.qe += '  {} = {},\n'.format(tag, '"none"')
        #     if tag == 'pseudo_dir':
        #         self.qe += '  {} = {},\n/\n\n'.format(
        #             tag, '"' + self.paw_path + '"')
        # for index, tag in enumerate(system_tags):
        #     if index == 0:
        #         self.qe += '&SYSTEM\n'
        #     if tag == 'ibrav':
        #         self.qe += '  {} = {},\n'.format(tag, '0')
        #     if tag == 'nat':
        #         self.qe += '  {} = {},\n'.format(tag, self.values[
        #             'number_of_atoms'])
        #     if tag == 'ntyp':
        #         self.qe += '  {} = {},\n'.format(tag,
        #                                          len(self.values['number']))
        #     if tag == 'starting_magnetization':
        #         magmom = [float(x[-1]) for x in self.values['MAGMOM']]
        #         maximum_magmom = max(magmom)
        #         for i, x in enumerate(magmom):
        #             self.qe += '  {}({}) = {},\n'.format(tag, i + 1, float(
        #                 x) / maximum_magmom)
        #     if tag == 'ecutwfc':
        #         self.qe += '  {} = {},\n'.format(tag, self.values['ENCUT'])
        #     if tag == 'ecutrho':
        #         self.qe += '  {} = {},\n'.format(tag,
        #                                          self.values['ENCUT'] * 12)
        #     if tag == 'nosym':
        #         if self.values['ISYM'] == 0:
        #             self.qe += '  {} = {},\n'.format(tag, '.TRUE.')
        #     if tag == 'occupations':
        #         if self.values['ISMEAR'] == 0:
        #             self.qe += '  {} = {},\n'.format(tag, '"smearing"')
        #         else:
        #             raise ValueError('ISMEAR is not 0')
        #     if tag == 'degauss':
        #         self.qe += '  {} = {},\n'.format(tag, self.values['SIGMA'])
        #     if tag == 'smearing':
        #         if self.values['ISMEAR'] == 0:
        #             self.qe += '  {} = {},\n'.format(tag, '"gaussian"')
        #         else:
        #             raise ValueError('ISMEAR is not 0')
        #     if tag == 'nspin':
        #         if self.values['ISPIN'] == 1:
        #             self.qe += '  {} = {},\n'.format(tag, 1)
        #         elif self.values['ISPIN'] == 2:
        #             self.qe += '  {} = {},\n'.format(tag, 2)
        #         else:
        #             raise ValueError('ISPIN is not 1 or 2')
        #     if tag == 'input_dft':
        #         if self.values['METAGGA'].lower() == 'scan':
        #             self.qe += '  {} = {},\n'.format(tag, '"scan"')
        #         else:
        #             raise ValueError('METAGGA is not SCAN')
        #     if tag == 'lda_plus_u':
        #         if self.values['LDAU']:
        #             self.qe += '  {} = {},\n'.format(tag, '.TRUE.')
        #     if tag == 'lda_plus_u_kind':
        #         self.qe += '  {} = {},\n'.format(tag, '0')
        #     if tag == 'Hubbard_U':
        #         for i, x in enumerate(self.values['LDAUU']):
        #             self.qe += '  {}({}) = {},\n'.format(tag, i + 1, x)
        # self.qe += '/\n\n'
        # for index, tag in enumerate(electrons_tags):
        #     if index == 0:
        #         self.qe += '&ELECTRONS\n'
        #     if tag == 'electron_maxstep':
        #         self.qe += '  {} = {},\n'.format(tag, self.values['NELM'])
        #     if tag == 'conv_thr':
        #         self.qe += '  {} = {},\n/\n\n'.format(tag,
        #                                               self.values['EDIFF'])
        #
        # # IONS
        # if self.values['ISIF'] in [2, 3]:
        #     self.qe += '&IONS\n/\n\n'
        #
        # # CELL
        # if self.values['ISIF'] == 3:
        #     self.qe += '&CELL\n/\n\n'
        #
        # # ATOMIC_SPECIES
        # self.qe += 'ATOMIC_SPECIES\n'
        # for x in self.values['symbols']:
        #     self.qe += '  {} {} {}\n'.format(x, 1, x + '.UPF')
        # self.qe += '\n'
        #
        # # K_POINTS
        # self.qe += 'K_POINTS automatic\n'
        # grid = self.values['grid']
        # if self.values['center'].lower() == 'gamma':
        #     self.qe += '  {} {} {} 0 0 0\n'.format(grid[0], grid[1], grid[2])
        # else:
        #     self.qe += '  {} {} {} 1 1 1\n'.format(grid[0], grid[1], grid[2])
        #
        # # CELL_PARAMETERS
        # self.qe += '\nCELL_PARAMETERS\n'
        # for x in self.values['lattice_vectors']:
        #     self.qe += '  {} {} {}\n'.format(x[0], x[1], x[2])

        # ATOMIC_POSITIONS
        self.qe = '\nATOMIC_POSITIONS crystal\n'
        k = 0
        coordinates = self.values['coordinates']
        for i, symbol in enumerate(self.values['symbols']):
            for j in range(self.values['number'][i]):
                self.qe += '  {} {} {} {}\n'.format(symbol, coordinates[k][0],
                                                    coordinates[k][1],
                                                    coordinates[k][2])
                k += 1

        # WRITE QE INPUT
        with open('qe.in', 'w') as file:
            file.write(self.qe)


def main():
    """
Main program
    """
    # INCAR
    incar = Incar(INCAR, TAGS, TYPES, UNITS)
    incar.read_incar()
    incar.parse_magmom()
    incar.clean_ldauu_values()
    values = incar.values

    # POSCAR/CONTCAR
    poscar = Poscar(POSCAR)
    poscar.read_poscar()
    values.update(poscar.values)

    # KPOINTS
    kpoints = Kpoints(KPOINTS)
    kpoints.read_kpoints()
    values.update(kpoints.values)

    # QE INPUT
    qe = Qe(values, PAW_PATH)
    qe.write_qe()


if __name__ == '__main__':
    main()
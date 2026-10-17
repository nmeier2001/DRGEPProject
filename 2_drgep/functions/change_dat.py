'''
Changing the data in the dat file to the reduced mechanism. The reduced mechanism is stored in the skeletal_mechanism.pkl file.
'''

import re
from pathlib import Path
import shutil
import pickle as pkl

BASE_DIR = Path(__file__).parent

# OUTPUT_BASE_DIR = BASE_DIR / 'RUNS'
FILE_PATH = BASE_DIR / 'Model'

class ChangeDat:
    '''
    Class to change the data in the dat file to the reduced mechanism. The reduced mechanism is stored in the skeletal_mechanism.pkl file.
    '''

    def __init__(self, dat_file_path, skeletal_mechanism_path, epsilon=0.0, sigma=0.01):
        self.dat_file_path = dat_file_path
        self.skeletal_mechanism_path = skeletal_mechanism_path
        self.epsilon = epsilon
        self.sigma = sigma

        self.copy_dat_file_path = Path(self.dat_file_path.parent)  / "modified_dat"
        print(f"Copy dat file path: {self.copy_dat_file_path}")

    # delete all files in the copy_dat_file_path folder
        if self.copy_dat_file_path.exists():
            for file in self.copy_dat_file_path.iterdir():
                if file.is_file():
                    file.unlink()

    # Load the skeletal mechanism from the pickle file
        with open(self.skeletal_mechanism_path, 'rb') as f:
            self.skeletal_rxn_info = pkl.load(f)
    # creat copy of the original dat file to change
        self.output_dat_path = self.copy_dat_file_path / f"model_{self.epsilon}_{self.sigma}.dat"
        shutil.copy(self.dat_file_path, self.output_dat_path)
    #read the copy dat file
        with open(self.output_dat_path, 'r', encoding='utf-8') as f:
            self.lines = f.read()

        self.change_species()
        self.change_reactions()

    def change_species(self):
        """
        Find the species between line 4 and 77 in the dat file and comment out the species that are not in the skeletal mechanism.
        """

        # create list of species from the skeletal mechanism
        skeletal_species = set()
        for rxn in self.skeletal_rxn_info.values():
            skeletal_species.update(rxn['reactants'])
            skeletal_species.update(rxn['products'])

        # find species between line 4 and 77 in dat file and comment out the species that are not in the skeletal mechanism
        for i in range(4, 77):
            line = self.lines.splitlines()[i]
            species = line.split()[0]
            if species not in skeletal_species:
                if not line.lstrip().startswith('!'):
                    self.lines = self.lines.replace(line, f'! {line}')
            else:
                self.lines = self.lines.replace(line, f'{line}')

        # update the dat file with the new lines
        with open(self.output_dat_path, 'w', encoding='utf-8') as f:
            f.write(self.lines)


    import re

    def change_reactions(self):
        """
        Kommentiert Reaktionen aus die nicht im Skelettmechanismus sind.
        Identifikation über die normalisierte Reaktionsgleichung (Edukte > Produkte),
        nicht über den Label-Tag – robuster gegenüber fehlenden oder mehrdeutigen Labels.
        """

        skeletal_rxn_strings = {
            self._normalize_reaction(info['reaction'])
            for info in self.skeletal_rxn_info.values()
        }

        with open(self.output_dat_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        rxn_pattern = re.compile(r'^!?\s*(.+?)\s*:')

        for i in range(91, len(lines)):
            line     = lines[i]
            stripped = line.lstrip()

            # Nur echte Reaktionszeilen – müssen '>' enthalten
            if '>' not in line:
                continue

            is_commented = stripped.startswith('!')

            match = rxn_pattern.match(stripped.lstrip('!').lstrip())
            if not match:
                continue

            rxn_str      = match.group(1).strip()
            rxn_norm     = self._normalize_reaction(rxn_str)
            in_skeletal  = rxn_norm in skeletal_rxn_strings

            if in_skeletal and is_commented:
                lines[i] = self._remove_leading_comment(line)
            elif not in_skeletal and not is_commented:
                lines[i] = f'!{line}'

        with open(self.output_dat_path, 'w', encoding='utf-8') as f:
            f.writelines(lines)

    @staticmethod
    def _normalize_reaction(rxn_str: str) -> str:
        """
        Normalisiert eine Reaktionsgleichung für robusten Vergleich:
        - Whitespace entfernen
        - Großschreibung vereinheitlichen
        - Edukte und Produkte jeweils sortieren (A+B = B+A)

        'E + CO2 > E + CO + O'  →  'CO2+E>CO+E+O'
        'E+CO2>E+CO+O'          →  'CO2+E>CO+E+O'
        """
        rxn_str = rxn_str.strip()

        if '!' in rxn_str:
            rxn_str = rxn_str[:rxn_str.index('!')].strip()

        if '>' not in rxn_str:
            return rxn_str

        lhs, rhs = rxn_str.split('>', 1)

        reactants = sorted(s.strip() for s in lhs.split('+') if s.strip())
        products  = sorted(s.strip() for s in rhs.split('+') if s.strip())

        return '+'.join(reactants) + '>' + '+'.join(products)

    @staticmethod
    def _remove_leading_comment(line: str) -> str:
        """Entfernt nur das erste '!' am Zeilenanfang."""
        stripped   = line.lstrip(' \t')
        leading_ws = line[:len(line) - len(stripped)]
        return leading_ws + stripped[1:]


if __name__ == "__main__":
    dat_file_path = FILE_PATH / 'model.dat'
    skeletal_mechanism_path = BASE_DIR / 'skeletal_mechanism.pkl'

    changer = ChangeDat(dat_file_path, skeletal_mechanism_path, epsilon=0.0, sigma=0.01)
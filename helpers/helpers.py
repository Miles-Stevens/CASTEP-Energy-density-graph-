from Imports import *
from ase import Atoms as ASEAtoms
from ase.io import read as ase_read
from ase.io import write as ase_write

class Helpers:
    def __init__(self):
        self._structure_cache = {}
        self._neighbor_cache = {}
        
    def _get_atomic_mass(self, element: str) -> float:
        """Get atomic mass with LRU caching"""
        return atomic_masses[atomic_numbers[element]]

    def _is_experimental_file(self, file_path: Path) -> bool:
        """Filename-based check for the experimental reference structure"""
        return 'experimental' in Path(file_path).name.lower()
    
    def _get_cached_structure(self, file_path: Path) -> Optional[ASEAtoms]:
        """Get cached structure or read and cache it"""
        if file_path in self._structure_cache:
            return self._structure_cache[file_path]
        
        try:
            with open(file_path, 'r') as f:
                content = f.read()
            
            # Use StringIO for in-memory processing
            atoms = ase_read(StringIO(content), format='res', index= 0)
            if isinstance(atoms, list):
                if not atoms:
                    return None
                atoms = atoms[0]

            self._structure_cache[file_path] = atoms
            return atoms
        except Exception as e:
            print(f"Error reading structure from {file_path}: {e}")
            return None
    
    def _get_cached_neighbors(self, atoms: ASEAtoms, cutoff: float) -> Dict[str, np.ndarray]:
        """Get cached neighbor list or compute and cache it"""
        # Create a simple hash for the atoms object based on positions and cell
        atoms_hash = hash((
            tuple(atoms.get_positions().flatten()),
            tuple(atoms.get_cell().array.flatten()) if atoms.cell is not None else (),
            tuple(atoms.get_chemical_symbols()),
            cutoff
        ))
        
        if atoms_hash in self._neighbor_cache:
            return self._neighbor_cache[atoms_hash]
        
        try:
            i, j, d = neighbor_list('ijd', atoms, cutoff=cutoff, self_interaction=False)
            neighbors = {'i': i, 'j': j, 'd': d}
            self._neighbor_cache[atoms_hash] = neighbors
            return neighbors
        except Exception:
            return {'i': np.array([]), 'j': np.array([]), 'd': np.array([])}
    
    def create_cif_from_file(self, file_path: Path) -> bool:
        """Create CIF file from input file with improved error handling and validation"""
        temporary_path = None
        try:
            base_cif = file_path.with_suffix('.cif')
            
            with open(file_path, 'r') as f:
                content = f.read()

            with tempfile.NamedTemporaryFile(mode='w+', suffix='.res', delete=False) as tmp:
                tmp.write(content)
                temporary_path = tmp.name

            atoms = ase_read(temporary_path, format='res', index=0)
            if isinstance(atoms, list):
                if not atoms:
                    raise ValueError(f'No structures found in {file_path}')
                atoms = atoms[0]
            if not isinstance(atoms, ASEAtoms):
                raise TypeError(f'Expected an ASE Atoms object, got {type(atoms).__name__}')

            ase_write(base_cif, atoms, format="cif")
            return True
        except Exception as e:
            print(f"Error creating CIF from {file_path}: {e}")
            return False
        finally:
            if temporary_path is not None:
                try:
                    os.unlink(temporary_path)
                except FileNotFoundError:
                    pass
    
    def _calculate_stoichiometry(self, atom_counts: Counter) -> str:
        """Calculate reduced stoichiometry from atom counts with validation"""
        if not atom_counts:
            return "?"
        
        # Remove any zero or negative counts
        valid_counts = {el: count for el, count in atom_counts.items() if count > 0}
        if not valid_counts:
            return "?"
        
        # More efficient calculation
        counts = list(valid_counts.values())
        if not counts:
            return "?"
        if len(counts) == 1:
            gcd_value = counts[0]
        else:
            gcd_value = counts[0]
            for count in counts[1:]:
                gcd_value = gcd(gcd_value, count)
        
        # Ensure gcd is at least 1
        if gcd_value <= 0:
            gcd_value = 1
        
        return ''.join(
            f"{el}{c // gcd_value if c // gcd_value > 1 else ''}" 
            for el, c in sorted(valid_counts.items())
        )
    
    def _get_gcd_from_atom_counts(self, atom_counts: Counter) -> int:
        """Get GCD from atom counts with validation"""
        if not atom_counts:
            return 1
        
        # Remove any zero or negative counts
        valid_counts = [count for count in atom_counts.values() if count > 0]
        if not valid_counts:
            return 1
        
        if len(valid_counts) == 1:
            return max(1, valid_counts[0])
        
        result = valid_counts[0]
        for count in valid_counts[1:]:
            result = gcd(result, count)
        
        return max(1, result)
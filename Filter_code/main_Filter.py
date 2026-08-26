from pathlib import Path
import os
import tempfile
from ase import Atoms
from ase.io import read
from helpers.helpers import Helpers
from Filter_code import Cyanurate_filter, Perovskite_filter


class Filter(Helpers):
    def __init__ (self): ## not sure if they will live here or go somewhere else 
        super().__init__()
        self.apply_cyanurate_filter = False
        self.apply_keep_octahedral = False

    def main(self, file_path: Path) -> bool:
            """OPTIMIZED: Apply structural filters with caching and early returns"""
            try:
                # OPTIMIZATION 1: Early return if no filters are active
                active_filters = [
                    self.apply_cyanurate_filter,
                    self.apply_keep_octahedral,
                ]
                
                if not any(active_filters):
                    return False

                # OPTIMIZATION 2: Use cached structure
                try:
                    with open(file_path, 'r') as f:
                        content = f.read()

                    with tempfile.NamedTemporaryFile(mode='w+', suffix='.res', delete=False) as tmp:
                        tmp.write(content)
                        tmp_path = tmp.name

                    try:
                        atoms = read(tmp_path, format='res', index=0)
                        if isinstance(atoms, list):
                            if not atoms:
                                return False
                            atoms = atoms[0]
                        if not isinstance(atoms, Atoms):
                            return False
                    finally:
                        os.unlink(tmp_path)
                except Exception:
                    return False
                    
                # OPTIMIZATION 3: Only compute neighbor lists for required cutoffs
                required_cutoffs = set()
                if self.apply_cyanurate_filter:
                    required_cutoffs.add(1.8)
                if self.apply_keep_octahedral:
                    required_cutoffs.add(3.0)


                # OPTIMIZATION 4: Compute only required neighbor lists with caching
                cutoff_map = {}
                for cutoff in required_cutoffs:
                    cutoff_map[cutoff] = self._get_cached_neighbors(atoms, cutoff)

                # OPTIMIZATION 5: Order filters by computational cost (cheapest first)
                # This allows for early exits
                
                # Cheapest: Check for presence of specific elements first
                if self.apply_cyanurate_filter:
                    symbols = atoms.get_chemical_symbols()
                    if 'C' not in symbols:
                        # No carbon atoms, can't have cyanurate structure
                        pass
                    elif Cyanurate_filter.Filter().has_cyanurate_structure(atoms, cutoff_map[1.8]):
                        return True
                
                if self.apply_keep_octahedral:
                    symbols = atoms.get_chemical_symbols()
                    if 'In' not in symbols:
                        # No indium atoms, filter passes
                        pass
                    elif not Perovskite_filter.Filter().has_Perovskite_indium(atoms, cutoff_map[3.0]):
                        return True
                

                return False

            except Exception as e:
                print(f"Error filtering {file_path}: {e}")
                return False


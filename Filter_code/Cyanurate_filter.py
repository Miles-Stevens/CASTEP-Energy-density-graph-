# this filter checks to see if the res files has a Cyanurate strucutre (rings) that we don't care about in this case
import numpy as np
from ase import Atoms

class Filter:
    def  has_cyanurate_structure(self,atoms, neighbors) -> bool: 
        """OPTIMIZED: Check for cyanurate structure"""
        try:
            symbols = Atoms.get_chemical_symbols(atoms)
            carbon_indices = [i for i, sym in enumerate(symbols) if sym == 'C']
            if not carbon_indices:
                return False
            
            # OPTIMIZATION: Use numpy for faster counting
            i_array = neighbors['i']
            if len(i_array) == 0:
                return False
            
            # Count coordination numbers only for carbon atoms
            carbon_coords = np.bincount(i_array[np.isin(i_array, carbon_indices)], 
                                        minlength=len(symbols))
            return bool (np.any(carbon_coords[carbon_indices] > 2, axis= None))
            
        except Exception:
            return False
# This filter looks to see if the struicture has a Perovskiter structure
#  with 1 indium surounded by 6 nitrogen or sulphur.

import numpy as np
from ase import Atoms
class Filter:
    def has_Perovskite_indium( self, atoms, neighbors) -> bool:
        """OPTIMIZED: Check for octahedral indium"""
        try:
            symbols = np.array(Atoms.get_chemical_symbols(atoms))
            indium_mask = symbols == 'In'
            if not np.any(indium_mask):
                return False
            
            i_array = neighbors['i']
            if len(i_array) == 0:
                return False
            
            # OPTIMIZATION: Use numpy operations
            indium_indices = np.where(indium_mask)[0]
            coordination = np.bincount(i_array, minlength=len(symbols))
            return bool(np.any(coordination[indium_indices] == 6))
            
        except Exception:
            return False

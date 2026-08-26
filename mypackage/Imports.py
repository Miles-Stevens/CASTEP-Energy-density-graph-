# this will hae all the needed imports for all the files 

import warnings
import os
import time
import subprocess
from pathlib import Path
from math import gcd
from itertools import cycle
from collections import Counter
from functools import partial, lru_cache
from typing import List, Tuple, Optional, Dict, Any
import concurrent.futures
import tempfile
import numpy as np
from scipy import constants, spatial
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm
import ipywidgets as widgets
from IPython.display import display

# Scientific computing imports
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from ase import io, Atoms
from io import StringIO
from ase.neighborlist import neighbor_list
from ase.io import read
from ase.data import atomic_masses
from ase.data import atomic_numbers
from ase.geometry import cellpar_to_cell
from ase.geometry import distance
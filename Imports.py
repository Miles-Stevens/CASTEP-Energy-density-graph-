
# these are the 3rd party libaries that are used in the program # I want to decrease these
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

from ase import io, Atoms
from io import StringIO
from ase.neighborlist import neighbor_list
from ase.io import read as ase_write
from ase.data import atomic_masses
from ase.data import atomic_numbers
from ase.geometry import cellpar_to_cell
from ase.geometry import distance



from helpers.helpers import Helpers
from itertools import cycle

class Graph(Helpers):
    
    def __init__(self):
        self.marker_cycle = cycle(['o', '*', '^', 'D', 'P', 'v', '<', '>', 's', 'X'])
        self.fu_to_color = {
            1: 'red', 2: 'green', 3: 'blue',
            4: 'purple', 6: 'magenta', 8: 'cyan',
            'other': 'black', 'Filtered out': 'lightgray'
        }
        self.vesta_path = r"C:\\VESTA-win64\\VESTA-win64\\VESTA.exe"


 
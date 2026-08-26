import numpy as np 
from itertools import cycle
import matplotlib.pyplot as plt
import ipywidgets as widgets
from IPython.display import display


from helpers.helpers import Helpers

%matplotlib widget


class Graph(Helpers):
    
    def __init__(self):
        self.marker_cycle = cycle(['o', '*', '^', 'D', 'P', 'v', '<', '>', 's', 'X'])
        self.fu_to_color = {
            1: 'red', 2: 'green', 3: 'blue',
            4: 'purple', 6: 'magenta', 8: 'cyan',
            'other': 'black', 'Filtered out': 'lightgray'
        }
        self.vesta_path = r"C:\\VESTA-win64\\VESTA-win64\\VESTA.exe"

    def create_plot(self, file_names: np.ndarray, energies: np.ndarray, 
                            densities: np.ndarray, formula_units: np.ndarray, 
                            stoichiometries: np.ndarray) -> None:
        """Create interactive plot with corrected hover labels"""
        # Create marker mapping
        unique_stoichiometries = sorted(set(stoichiometries))
        stoich_to_marker = {stoich: next(self.marker_cycle) for stoich in unique_stoichiometries}
        
        # Group data
        fu_labels = np.array([
            fu if fu in self.fu_to_color else 'other' 
            for fu in formula_units
        ])
        
        group_keys = np.array([f"{s}_{fu}" for s, fu in zip(stoichiometries, fu_labels)])
        unique_groups = np.unique(group_keys)
        
        grouped = {}
        for group in unique_groups:
            mask = group_keys == group
            indices = np.where(mask)[0]
            stoich, fu_label = group.rsplit('_', 1)
            grouped[(stoich, fu_label)] = indices
        
        hull_energies, hull_densities = self.Create_convex_hull(densities=densities, energies=energies)
        possible_energies= hull_energies+30


        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        plt.subplots_adjust(right=0.75)
        
        scatter_points = []
        ax.plot(hull_densities, hull_energies, label= 'leading edge', c='k', linestyle= '--')
        ax.plot(hull_densities, possible_energies, label= '30 meV/atom above leading edge', c='r', linestyle= '--')
        for (stoich, fu_label), idxs in grouped.items():
            marker = stoich_to_marker[stoich]
            if fu_label == 'Filtered out':
                color = self.fu_to_color['Filtered out']
            else:
                try:
                    fu_int = int(fu_label)
                    color = self.fu_to_color.get(fu_int, self.fu_to_color['other'])
                except ValueError:
                    color = self.fu_to_color['other']
            alpha = 0.8 if fu_label != 'Filtered out' else 0.5
            size = 60 if fu_label != 'Filtered out' else 30
            sc = ax.scatter(
                densities[idxs], energies[idxs],
                marker=marker, color=color, alpha=alpha,
                label=f"{stoich} | {fu_label} FU", picker=True, s=size ,edgecolors='black', linewidth=0.5
            )
            scatter_points.append((sc, idxs))
        
        # Add legend
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(
            by_label.values(), by_label.keys(),
            loc='upper left', bbox_to_anchor=(1.05, 1),
            title="Stoichiometry | F.U."
        )
        
        # Set labels with correct energy type
        ax.set_xlabel('Density (g/cm³)')
        ylabel = ('Relative energy (meV/atom)' if self.config.relative_energies_on 
                else 'Energy per atom (meV/atom)')
        ax.set_ylabel(ylabel)
        
        # Add interactivity with corrected labels
        self._add_plot_interactivity(fig, ax, scatter_points, file_names, 
                                        densities, energies, stoichiometries, fu_labels)
        
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0,4)
        ax.set_ylim(min(energies)-10, max(energies)+10)
        plt.tight_layout()
        plt.show()
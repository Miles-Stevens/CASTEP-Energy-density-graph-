import numpy as np
from scipy import spatial


class graph():
    def densify_polyline_xy(self, xv, yv, N):
            xv = np.asarray(xv, dtype=float)
            yv = np.asarray(yv, dtype=float)

            # Sort by x
            order = np.argsort(xv)
            xv = xv[order]
            yv = yv[order]

            # Base grid
            x_dense = np.linspace(xv[0], xv[-1], N)

            # FORCE vertex x-values into grid
            x_all = np.unique(np.concatenate([x_dense, xv]))

            # Interpolate
            y_all = np.interp(x_all, xv, yv)

            return x_all, y_all
            
    def Create_convex_hull(self, densities: np.ndarray, energies: np.ndarray):
        """Create convex hull from densities and energies"""
        #make coppies of the arrys
        densities_plus= densities.copy()
        energies_plus= energies.copy()
        # add points to connect them to the axies 
        energies_plus=np.append(energies_plus,(1050, 1050))
        densities_plus=np.append(densities_plus,(0, 4))
        points = np.column_stack((densities_plus, energies_plus))

        #create convex hull 
        hull = spatial.ConvexHull(points)
        vertecies= hull.vertices
        hull_points=points[vertecies]
        hull_points= hull_points[np.argsort(hull_points[:,0])] # this is a sort based on the density values
        hull_desnity, hull_energies= self.densify_polyline_xy(hull_points[:,0], hull_points[:,1], len(densities_plus))

        return hull_energies, hull_desnity
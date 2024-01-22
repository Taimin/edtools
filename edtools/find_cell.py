import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from collections import defaultdict
import yaml
from .utils import volume
from typing import Tuple

def unify(cell):
    if cell[3] > 90:
        cell[3] = 180 - cell[3]
    if cell[4] > 90:
        cell[4] = 180 - cell[4]
    if cell[5] > 90:
        cell[5] = 180 - cell[5]

def to_niggli_cell(cell):
    L = from_uc_to_lattice_matrix(cell)
    L = niggli_reduce_3d(L)
    cell = from_lattice_matrix_to_uc(L)
    return cell

def from_uc_to_lattice_matrix(cell):
    lattice = np.zeros((3, 3))
    cell[3] *= np.pi/180
    cell[4] *= np.pi/180
    cell[5] *= np.pi/180
    lattice[0, 0] = cell[0]
    lattice[1, 0] = cell[1] * np.cos(cell[5])
    lattice[2, 0] = cell[2] * np.cos(cell[4])
    lattice[1, 1] = cell[1] * np.sin(cell[5])
    lattice[1, 2] = -cell[2] * (np.cos(cell[4]) * np.cos(cell[5]) - np.cos(cell[3])) / np.sin(cell[5])
    lattice[2, 2] = cell[2]
    cell[3] /= np.pi/180
    cell[4] /= np.pi/180
    cell[5] /= np.pi/180
    return lattice

def from_lattice_matrix_to_uc(lattice):
    cell = [0, 0, 0, 0, 0, 0]
    cell[0] = lattice[0, 0]
    cell[2] = lattice[2, 2]
    try:
        gamma = np.arctan(lattice[1, 1] / lattice[1, 0])
    except ZeroDivisionError:
        gamma = np.pi
    if gamma < 0:
        gamma = gamma + np.pi
    cell[5] = np.rad2deg(gamma)
    cell[1] = lattice[1, 1] / np.sin(gamma)
    beta = np.arccos(lattice[2, 0] / lattice[2, 2])
    cell[4] = np.rad2deg(beta)
    alpha = lattice[2, 1] / lattice[2, 2] * np.sin(gamma) + np.cos(beta) * np.cos(gamma)
    cell[3] = np.rad2deg(alpha)
    return cell

def niggli_reduce_3d(lattice, eps: float=1e-5, loop_max=100):
    """
    niggli reduction

    :param lattice: Origin lattice
    :type lattice: list or np.ndarray with 9 elements
    :param eps: tolerance
    :type eps: float default 1e-5
    :return: a reduced lattice
    :rtype: 3x3 np.ndarray
    """
    try:
        reduced_lattice = np.array(lattice).reshape((3, 3))
    except:
        raise ValueError("Must 9 elements for 3x3 lattice")

    L = reduced_lattice
    G = _get_metric(L)

    # This sets an upper limit on the number of iterations.
    for _ in range(loop_max):
        reduced = True
        # step 0: get parameters for A1-A8
        # X, E, Z for xi, eta, zeta respectively
        G = _get_metric(L)
        A, B, C, X, E, Z = _get_G_param(G)

        # step 1
        if A > B + eps or (abs(A - B) < eps and abs(X) > abs(E) + eps):
            M = np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]])
            L = np.matmul(M, L)
            reduced = False
        # step 2
        if (B > C + eps) or (abs(B - C) < eps and abs(E) > abs(Z) + eps):
            M = np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]])
            L = np.matmul(M, L)
            reduced = False
            continue

        l, m, n = _get_angle_param(X, E, Z, eps)
        # step 3
        if l * m * n == 1:
            i, j, k = l, m, n
            M = np.array([[i, 0, 0], [0, j, 0], [0, 0, k]])
            L = np.matmul(M, L)
            reduced = False
        # step 4
        elif l * m * n == 0 or l * m * n == -1:
            i = -1 if l == 1 else 1
            j = -1 if m == 1 else 1
            k = -1 if n == 1 else 1

            if i * j * k == -1:
                if n == 0:
                    k = -1
                elif m == 0:
                    j = -1
                elif l == 0:
                    i = -1
            M = np.array([[i, 0, 0], [0, j, 0], [0, 0, k]])
            L = np.matmul(M, L)
            reduced = False

        G = _get_metric(L)
        A, B, C, X, E, Z = _get_G_param(G)
        # step 5
        if (
            abs(X) > B + eps
            or (abs(X - B) < eps and 2 * E < Z - eps)
            or (abs(X + B) < eps and Z < -eps)
        ):
            M = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [0, -np.sign(X), 1]])
            L = np.matmul(M, L)
            reduced = False
            continue

        # step 6
        if (
            abs(E) > A + eps
            or (abs(A - E) < eps and 2 * X < Z - eps)
            or (abs(A + E) < eps and Z < -eps)
        ):
            M = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [-np.sign(E), 0, 1]])
            L = np.matmul(M, L)
            reduced = False
            continue

        # step 7
        if (
            abs(Z) > A + eps
            or (abs(A - Z) < eps and 2 * X < E - eps)
            or (abs(A + Z) < eps and E < -eps)
        ):
            M = np.array([[1, 0, 0],
                          [-np.sign(Z), 1, 0],
                          [0, 0, 1]])
            L = np.matmul(M, L)
            reduced = False
            continue

        # step 8
        if X + E + Z + A + B < -eps or (abs(X + E + Z + A + B) < eps < Z + (A + E) * 2):
            M = np.array([[1, 0, 0],
                          [0, 1, 0],
                          [1, 1, 1]])
            L = np.matmul(M, L)
            reduced = False
            continue

        if reduced:
            break

    reduced_lattice = L
    return reduced_lattice

def niggli_check(L, eps=1e-5) -> bool:
    G = _get_metric(L)
    A, B, C, X, E, Z = _get_G_param(G)
    return _niggli_check(A, B, C, X, E, Z, eps)

def _niggli_check(A, B, C, X, E, Z, eps=1e-5):
    """Checks that the niggli reduced cell satisfies the niggli conditions.
    Conditions listed at: https://arxiv.org/pdf/1203.5146.pdf.
    Args:
        A (float): a.a
        B (float): b.b
        C (float): c.c
        xi (float): 2*b.c
        eta (float): 2*c.a
        zeta (float): 2*a.b

    Returns:
        False if niggli conditons aren't met.
    """
    if not (A-eps > 0 and (A < B-eps or np.allclose(A,B,atol=eps)) and
            (B < C-eps or np.allclose(B,C,atol=eps))):
        return False

    if np.allclose(A,B,atol=eps) and not (abs(X) < abs(E)-eps or
                                          np.allclose(abs(X),abs(E),atol=eps)):
        return False

    if np.allclose(B,C,atol=eps) and not (abs(E) < abs(Z)-eps
                                          or np.allclose(abs(E),abs(Z),atol=eps)):
        return False

    if not ((X-eps > 0 and E-eps > 0 and Z-eps > 0) or
            ((X < 0-eps or np.allclose(X,0,atol=eps))
             and (E < 0-eps or np.allclose(E,0,atol=eps))
             and (Z < 0-eps or np.allclose(Z,0,atol=eps)))):
        return False

    if not (abs(X) < B-eps or np.allclose(abs(X),B,atol=eps)):
        return False

    if not ((abs(E) < A-eps or np.allclose(abs(E),A,atol=eps)) and (abs(Z) < A-eps or
                                                           np.allclose(abs(Z),A,atol=eps))):
        return False

    if not (C < A+B+C+X+E+Z-eps or np.allclose(C, A+B+C+X+E+Z,atol=eps)):
        return False

    if np.allclose(X,B,atol=eps) and not (Z < 2.*E-eps or
                                           np.allclose(Z,2.*E,atol=eps)):
        return False

    if np.allclose(E,A,atol=eps) and not (Z < 2.*X-eps or
                                            np.allclose(Z,2.*X,atol=eps)):
        return False

    if np.allclose(Z,A,atol=eps) and not (E < 2.*X-eps or
                                             np.allclose(E,2.*X,atol=eps)):
        return False

    if np.allclose(X,-B,atol=eps) and not np.allclose(Z,0,atol=eps):
        return False

    if np.allclose(E,-A,atol=eps) and not np.allclose(Z,0,atol=eps):
        return False

    if np.allclose(Z,-A,atol=eps) and not np.allclose(E,0,atol=eps):
        return False

    if np.allclose(C,A+B+C+X+E+Z,rtol=0.0) and not ((2.*A+2.*E+Z) < 0-eps or
                                                 np.allclose(2.*A+2.*E+Z,0,atol=eps)):
        return False

    return True

def _get_G_param(G, eps=1e-5) -> Tuple[float, float, float, float, float, float]:
    """
    A = a.a = G[0, 0]
    B = b.b = G[1, 1]
    C = c.c = G[2, 2]
    xi = 2.b.c = 2 * G[1, 2]
    eta = 2.c.a = 2 * G[0, 2]
    zeta = 2.a.b = 2 * G[0, 1]
    """
    A = G[0, 0]
    B = G[1, 1]
    C = G[2, 2]
    X = 2 * G[1, 2]
    E = 2 * G[0, 2]
    Z = 2 * G[0, 1]

    return A, B, C, X, E, Z

def _get_angle_param(X, E, Z, eps=1e-5) -> Tuple[int, int, int]:
    # angle types
    l = _get_angle_type(X, eps)
    m = _get_angle_type(E, eps)
    n = _get_angle_type(Z, eps)

    return l, m, n


def _get_angle_type(angle, eps=1e-5) -> int:
    """
    -------------
    Angle   value
    -----   -----
    Acute   1
    Obtuse  -1
    Right   0
    -------------
    """
    if angle < -eps:
        return -1
    elif angle > eps:
        return 1
    else: # -eps < angle < eps
        return 0

def _get_metric(lattice) -> np.ndarray:
    M = lattice.reshape((3, 3))
    return np.matmul(M, M.T)

def weighted_average(values, weights=None):
    """Returns weighted mean and standard deviation"""
    mean = np.average(values, weights=weights)
    variance = np.average((values - mean)**2, weights=weights)
    std = np.sqrt(variance)
    return mean, std


def parse_cellparm(fn):
    cells = []
    weights = []
    
    with open(fn, "r") as f:
        for line in f:
            line = line.strip().rsplit("!")[0].upper()
            if not line:
                continue
            line = line.replace("=", "")
            line = line.split()
            i = line.index("UNIT_CELL_CONSTANTS")
            cell = [float(val) for val in line[i+1:i+7]]
            try:
                j = line.index("WEIGHT")
                weight = int(line[j+1])
            except ValueError:
                weight = 1
            
            cells.append(cell)
            weights.append(weight)

    print(f"Loaded {len(cells)} unit cells from file {fn}")
    cells = np.array(cells)
    weights = np.array(weights)

    return cells, weights


def find_cell(cells, weights, binsize=0.5):
    """Opens a plot with 6 subplots in which the cell parameter histogram is displayed.
    It will calculate the weighted mean of the unit cell parameters. The ranges can be
    adjusted by dragging on the plots.
    """
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    ang_par = cells[:,3:6]
    ang_xlim = int(np.percentile(ang_par, 5)) - 2, int(np.percentile(ang_par, 95)) + 2
    
    latt_parr = cells[:,0:3]
    latt_xlim = int(np.percentile(latt_parr, 5)) - 2, int(np.percentile(latt_parr, 95)) + 2
    
    spans = {}
    lines = {}
    variables  = {}
    names = "a b c \\alpha \\beta \\gamma".split()
    params = {}
    
    def get_spanfunc(i, ax):
        def onselect(xmin, xmax):
            # print(i, xmin, xmax)
            update(i, ax, xmin, xmax)
            fig.canvas.draw()
        return onselect
    
    def update(i, ax, xmin, xmax):
        par, bins = variables[i]
        idx = (par > xmin) & (par < xmax)
        sel_par = par[idx]
        sel_w = weights[idx]
    
        if len(sel_par) == 0:
            mu, sigma = 0.0, 0.0
        else:
            mu, sigma = weighted_average(sel_par, sel_w)
    
        if i in lines:
            for item in lines[i]:
                try:
                    item.remove()
                except ValueError:
                    pass

        if sigma > 0:
            x = np.arange(xmin-10, xmax+10, binsize/2)
            y = stats.norm.pdf(x, mu, sigma)
            l = ax.plot(x, y, 'r--', linewidth=1.5)
            lines[i] = l
    
        name = names[i]
        ax.set_title(f"${name}$: $\mu={mu:.2f}$, $\sigma={sigma:.2f}$")
        params[i] = mu, sigma
        return mu, sigma
    
    k = binsize/2  # displace by half a binsize to center bins on whole values

    for i in range(6):
        ax = axes[i]
        
        par = cells[:,i]
    
        median = np.median(par)
        bins = np.arange(min(par)-1.0-k, max(par)+1.0-k, binsize)  # pad 1 in case par values are all equal

        n, bins, patches = ax.hist(par, bins, rwidth=0.8, density=True)

        variables[i] = par, bins
    
        mu, sigma = update(i, ax, median-2, median+2)
    
        ax.set_ylabel("Frequency")
        if i < 3:
            xlim = latt_xlim
            ax.set_xlabel("Length ($\mathrm{\AA}$)")
        if i >=3:
            xlim = ang_xlim
            ax.set_xlabel("Angle ($\mathrm{^\circ}$)")
        
        ax.set_xlim(*xlim)
        onselect = get_spanfunc(i, ax)
        
        span = SpanSelector(ax, onselect, 'horizontal', useblit=True,
                    rectprops=dict(alpha=0.2, facecolor='red'), span_stays=True, minspan=1.0)
        
        spans[i] = span  # keep a reference in memory
        params[i] = mu, sigma
    
    plt.show()

    constants, esds = list(zip(*params.values()))

    return constants, esds


def d_calculator(cell: list) -> tuple:
    """Helper function for `unit_cell_lcv_distance`"""
    a, b, c, alpha, beta, gamma = cell
    d_ab = np.sqrt(a**2 + b**2 - 2*a*b*np.cos(np.radians(180 - gamma)))
    d_ac = np.sqrt(a**2 + c**2 - 2*a*c*np.cos(np.radians(180 - beta)))
    d_bc = np.sqrt(b**2 + c**2 - 2*b*c*np.cos(np.radians(180 - alpha)))
    return d_ab, d_ac, d_bc


def unit_cell_lcv_distance(cell1: list, cell2: list) -> float:
    """Implements Linear Cell Volume from Acta Cryst. (2013). D69, 1617-1632"""
    d_ab1, d_ac1, d_bc1 = d_calculator(cell1)
    d_ab2, d_ac2, d_bc2 = d_calculator(cell2)
    M_ab = abs(d_ab1 - d_ab2)/min(d_ab1, d_ab2)
    M_ac = abs(d_ac1 - d_ac2)/min(d_ac1, d_ac2)
    M_bc = abs(d_bc1 - d_bc2)/min(d_bc1, d_bc2)
    return max(M_ab, M_ac, M_bc)


def get_clusters(z, cells, distance=0.5):
    clusters = fcluster(z, distance, criterion='distance')
    grouped = defaultdict(list)
    for i, c in enumerate(clusters):
        grouped[c].append(i)
    
    print("-"*40)
    np.set_printoptions(formatter={'float': '{:7.2f}'.format})
    for i in sorted(grouped.keys()):
        cluster = grouped[i]
        clustsize = len(cluster)
        if clustsize == 1:
            del grouped[i]
            continue
        print(f"\nCluster #{i} ({clustsize} items)")
        vols = []
        for j in cluster:
            cell = cells[j]
            vol = volume(cell)
            vols.append(vol)
            print(f"{j+1:5d} {cell}  Vol.: {vol:6.1f}")
        print(" ---")
        print("Mean: {}  Vol.: {:6.1f}".format(np.mean(cells[cluster], axis=0), np.mean(vols)))
        print(" Min: {}  Vol.: {:6.1f}".format(np.min(cells[cluster], axis=0), np.min(vols)))
        print(" Max: {}  Vol.: {:6.1f}".format(np.max(cells[cluster], axis=0), np.max(vols)))
    
    print("")

    return grouped


def distance_from_dendrogram(z, ylabel: str="", initial_distance: float=None) -> float:
    """Takes a linkage object `z` from scipy.cluster.hierarchy.linkage and displays a
    dendrogram. The cutoff distance can be picked interactively, and is returned
    ylabel: sets the label for the y-axis
    initial_distance: initial cutoff distsance to display
    """
    if initial_distance == None:
        # corresponding with MATLAB behavior
        distance = round(0.7*max(z[:,2]), 4)
    else:
        distance = initial_distance

    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    tree = dendrogram(z, color_threshold=distance, ax=ax)

    # use 1-based indexing for display by incrementing label
    _, labels = plt.xticks()
    for l in labels:
        l.set_text(str(int(l.get_text())+1))

    ax.set_xlabel("Index")
    ax.set_ylabel(f"Distance ({ylabel})")
    ax.set_title(f"Dendrogram (cutoff={distance:.2f})")
    hline = ax.axhline(y=distance)
    
    def get_cutoff(event):
        nonlocal hline
        nonlocal tree
        nonlocal distance

        if event:
            distance = round(event.ydata, 4)
            ax.set_title(f"Dendrogram (cutoff={distance:.2f})")
            hline.remove()
            hline = ax.axhline(y=distance)
        
            for c in ax.collections:
                c.remove()

            tree = dendrogram(z, color_threshold=distance, ax=ax)

            fig.canvas.draw()
    
    fig.canvas.mpl_connect('button_press_event', get_cutoff)
    plt.show()

    return distance


def volume_difference(cell1: list, cell2: list):
    """Return the absolute difference in volumes between two unit cells"""
    v1 = volume(cell1)
    v2 = volume(cell2)
    return abs(v1-v2)


def cluster_cell(cells: list, 
                 distance: float=None, 
                 method: str="average", 
                 metric: str="euclidean", 
                 use_radian: bool=False,
                 use_sine: bool=False):
    """Perform hierarchical cluster analysis on a list of cells. 
    
    method: complete, average, weighted, centroid, median, ward, single
    metric: lcv, volume, euclidean
    distance: cutoff distance, if it is not given, pop up a dendrogram to
        interactively choose a cutoff distance
    use_radian: Use radian instead of degrees to downweight difference
    use_sine: Use sine for unit cell clustering (to disambiguousize the difference in angles)
    """

    from scipy.spatial.distance import pdist

    if use_sine:
        _cells = to_sin(cells)
    elif use_radian:
        _cells = to_radian(cells)
    else:
        _cells = cells

    if metric == "lcv":
        dist = pdist(_cells, metric=unit_cell_lcv_distance)
        z = linkage(dist,  method=method)
        initial_distance = None
    elif metric == "volume":
        dist = pdist(_cells, metric=volume_difference)
        z = linkage(dist,  method=method)
        initial_distance = 250.0
    else:
        z = linkage(_cells,  metric=metric, method=method)
        initial_distance = 2.0

    if not distance:
        distance = distance_from_dendrogram(z, ylabel=metric, initial_distance=initial_distance)

    print(f"Linkage method = {method}")
    print(f"Cutoff distance = {distance}")
    print(f"Distance metric = {metric}")
    print("")

    return get_clusters(z, cells, distance=distance)

def to_radian(cells):
    """convert all angles in unit cell parameter list to radians
    cells: the cell parameters that are parsed from cells.yaml as np array"""
    cells_radian = cells.copy()
    cells_radian[:, 3:6] = np.radians(cells[:, 3:6])

    return cells_radian

def to_sin(cells):
    """convert all angles in unit cell parameter list to sine
    cells: the cell parameters that are parsed from cells.yaml as np array"""
    angles = cells[:, 3:6]
    angles_sine = np.sin(np.radians(angles))

    cells_sine = cells.copy()
    cells_sine[:, 3:6] = angles_sine
    # convert also the cell angles using arcsin in order to avoid the <> 90 degree ambiguity thingy
    cells[:, 3:6] = np.degrees(np.arcsin(cells_sine[:, 3:6]))

    return cells_sine


def put_in_order(cells):
    """order cell parameters in order to eliminate difference in cell distance because of parameter order"""
    ordered_cells = []
    for i in range(0, len(cells)):
        cell = cells[i, :]
        cell_lengths = np.array(cell[:3])
        cell_angles = np.array(cell[3:])
        cell_array = np.vstack((cell_lengths, cell_angles))
        sortedArr = cell_array[:, np.argsort(cell_array[0, :])]
        sortedArr = sortedArr.ravel()
        ordered_cells.append(sortedArr)
    return np.array(ordered_cells)


def main():
    import argparse

    description = "Program for finding the unit cell from a serial crystallography experiment."
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
        
    parser.add_argument("args",
                        type=str, nargs="*", metavar="FILE",
                        help="Path to cells.yaml file")

    parser.add_argument("-b","--binsize",
                        action="store", type=float, dest="binsize",
                        help="Binsize for the histogram, default=0.5")

    parser.add_argument("-c","--cluster",
                        action="store_true", dest="cluster",
                        help="Apply cluster analysis")

    parser.add_argument("-d","--distance",
                        action="store", type=float, dest="distance",
                        help="Cutoff distance to use for clustering, bypass dendrogram")

    parser.add_argument("-m","--method",
                        action="store", type=str, dest="method",
                        choices="single average complete median weighted centroid ward".split(),
                        help="Linkage algorithm to use (see `scipy.cluster.hierarchy.linkage`)")

    parser.add_argument("-t","--metric",
                        action="store", type=str, dest="metric",
                        choices="euclidean lcv volume".split(),
                        help="Metric for calculating the distance between items (see `scipy.cluster.hierarchy.linkage`)")

    parser.add_argument("-l", "--use_raw_cell",
                        action="store_true", dest="use_raw_cell",
                        help="Use the reduced cell or not")

    parser.add_argument("-r", "--use_radian_for_angles",
                        action="store_true", dest="use_radian_for_clustering",
                        help="Use radians for unit cell clustering (to downweight the difference in angles)")

    parser.add_argument("-s", "--use_sine_for_angles",
                        action="store_true", dest="use_sine_for_clustering",
                        help="Use sine for unit cell clustering (to disambiguousize the difference in angles)")
 
    parser.set_defaults(binsize=0.5,
                        cluster=False,
                        distance=None,
                        method="average",
                        metric="euclidean",
                        use_raw_cell=False,
                        use_radian_for_clustering=False,
                        use_sine_for_clustering=False)
    
    options = parser.parse_args()

    distance = options.distance
    binsize = options.binsize
    cluster = options.cluster
    method = options.method
    metric = options.metric
    use_raw_cell = options.use_raw_cell
    use_radian = options.use_radian_for_clustering
    use_sine = options.use_sine_for_clustering
    args = options.args

    if args:
        fn = args[0]
    else:
        fn = "cells.yaml"

    ds = yaml.load(open(fn, "r"), Loader=yaml.Loader)

    cells = np.array([d['unit_cell'] for d in ds])
    cells = put_in_order(cells)
    if use_raw_cell:
        # transform cell into reduced cell
        raw_cells = []
        for i in range(0, len(cells)):
            raw_cell = unify(cells[i, :])
            raw_cells.append(raw_cell)
    try:
        weights = np.array([d["indexed"] for d in ds])
    except:
        weights = np.array([d["weight"] for d in ds])

    if cluster:
        clusters = cluster_cell(cells, distance=distance, method=method, metric=metric, use_radian=use_radian, use_sine=use_sine)
        for i, idx in clusters.items():
            clustered_ds = [ds[i] for i in idx]
            fout = f"cells_cluster_{i}_{len(idx)}-items.yaml"
            yaml.dump(clustered_ds, open(fout, "w"))
            print(f"Wrote cluster {i} to file `{fout}`")
    
    else:
        constants, esds = find_cell(cells, weights, binsize=binsize)
        
        print()
        print("Weighted mean of histogram analysis")
        print("---")
        print("Unit cell parameters: ", end="")
        for c in constants:
            print(f"{c:8.3f}", end="")
        print()
        print("Unit cell esds:       ", end="")
        for e in esds:
            print(f"{e:8.3f}", end="")
        print()

        try:
            import uncertainties as u
        except ImportError:
            pass
        else:
            print()
            names = (("a"," Å"), ("b"," Å"), ("c"," Å"),
                     ("α", "°"), ("β", "°"), ("γ", "°"))
            for i, (c, e) in enumerate(zip(constants, esds)):
                name, unit = names[i]
                val = u.ufloat(c, e)
                end = ", " if i < 5 else "\n"
                print(f"{name}={val:.2uS}{unit}", end=end)

        print()
        print("UNIT_CELL_CONSTANTS= " + " ".join(f"{val:.3f}" for val in constants))


if __name__ == '__main__':
    main()

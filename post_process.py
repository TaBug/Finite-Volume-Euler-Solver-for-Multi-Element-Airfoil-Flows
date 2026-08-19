import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from readgri import readgri


def figPath(datFile, meshFile, kind):
    """fig/<dat stem>__<mesh stem>_<kind>.png, creating the folder if needed.

    Naming the figure after both inputs means a plot always says which solution
    and which mesh produced it, so runs cannot be mixed up later.
    """
    name = f"{Path(datFile).stem}__{Path(meshFile).stem}_{kind}.png"
    out = Path('fig') / name
    out.parent.mkdir(parents=True, exist_ok=True)
    return out


def savePlot(f, savePath):
    """Write the figure out. Must be called before plt.show(), which hands the
    figure to the GUI and can leave nothing left to save."""
    if savePath is None:
        return
    f.savefig(savePath, dpi=150, bbox_inches='tight')
    print(f"saved {savePath}")


def plotMach(U, mesh, savePath=None):
    grid = readgri(mesh)  # read once; readgri rebuilds a sparse edge hash per call
    elem = grid['E']
    nodes = grid['V']
    X, Y = nodes.T
    # read the state matrix
    gamma = 1.4
    r = U[:, 0]
    u = U[:, 1] / r
    v = U[:, 2] / r
    q = np.sqrt(u ** 2 + v ** 2)
    p = (gamma - 1) * (U[:, 3] - 0.5 * r * q ** 2)
    c = np.sqrt(gamma * p / r)
    M = q / c

    f = plt.figure(figsize=(12, 6))
    plt.tripcolor(X, Y, triangles=elem, facecolors=M, shading='flat', edgecolor='black')
    cbar = plt.colorbar(orientation='vertical')
    cbar.ax.tick_params()
    plt.set_cmap('jet')
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel('y (m)', fontsize=16)
    f.tight_layout()
    plt.xlim([-0.3, 1.3])
    plt.ylim([-0.4, 0.4])
    savePlot(f, savePath)
    plt.show()
    plt.close(f)


# plot the cp and calculates the lift and drag coefficient
# INPUTS: U = state matrix (4xNe)
#         Uinf = free-stream state (1x4)
#         B2E = boundary face to element mapping matrix (Nbx3)
#         Bn = boundary face normal vector (Nbx2)
#         elem = nodes to elements mapping matrix (Nex3)
#         nodes = nodes coordinates matrix (number of nodes x 2)
#         alpha = AoA
# OUTPUTS: cd = drag coefficient
#          cl = lift coefficient
#          x = x-coordinates of the boundary faces (Nbx1)
#          cp = pressure coefficient of the boundary faces (Nbx1)
WALL_NAMES = ['main', 'slat', 'flap']  # the airfoil bodies; the rest is farfield


def getCoefficients(U, Uinf, mesh, alphaDeg):
    grid = readgri(mesh)  # read once; readgri rebuilds a sparse edge hash on every call
    V = grid['V']
    BE = grid['BE']        # [node1, node2, element, boundary group] per boundary face
    Bname = grid['Bname']

    # read the state matrix
    gamma = 1.4
    r = U[:, 0]
    u = U[:, 1] / r
    v = U[:, 2] / r
    q = np.sqrt(u ** 2 + v ** 2)
    p = (gamma - 1) * (U[:, 3] - 0.5 * r * q ** 2)

    # free-stream state
    c = 1
    rinf = Uinf[0]
    uinf = Uinf[1] / rinf
    vinf = Uinf[2] / rinf
    qinf = np.sqrt(uinf ** 2 + vinf ** 2)
    pinf = (gamma - 1) * (Uinf[3] - 0.5 * rinf * qinf ** 2)

    # Groups have different face counts, so the per-body results stay in dicts
    # rather than being stacked into one ragged array.
    pB = {name: [] for name in WALL_NAMES}
    xB = {name: [] for name in WALL_NAMES}

    Fx = 0.0
    Fy = 0.0
    # BE already pairs each boundary face with its element, and orients the node
    # pair to follow that element's traversal - so the normal below is outward
    # without having to trust the node order in the .gri file.
    for node1, node2, iElem, group in BE:
        name = Bname[group]
        if name not in WALL_NAMES:
            continue  # the farfield box exerts no force on the body

        node1Coord = V[node1]
        node2Coord = V[node2]
        l = np.sqrt((node2Coord[0] - node1Coord[0]) ** 2 + (node2Coord[1] - node1Coord[1]) ** 2)

        pB[name].append(p[iElem])
        xB[name].append((node1Coord[0] + node2Coord[0]) / 2)

        nx = (node2Coord[1] - node1Coord[1]) / l
        ny = -(node2Coord[0] - node1Coord[0]) / l
        Fx += l * p[iElem] * nx
        Fy += l * p[iElem] * ny

    alpha = np.radians(alphaDeg)  # the caller works in degrees
    D = np.cos(alpha) * Fx + np.sin(alpha) * Fy
    L = -np.sin(alpha) * Fx + np.cos(alpha) * Fy

    qdyn = 0.5 * rinf * qinf ** 2
    cl = L / (qdyn * c)
    cd = D / (qdyn * c)
    cp = {name: (np.array(pB[name]) - pinf) / qdyn for name in WALL_NAMES}
    x = {name: np.array(xB[name]) for name in WALL_NAMES}
    return cd, cl, x, cp


def computeFreestreamState(Minf, alpha):
    alpha = np.radians(alpha)
    gamma = 1.4
    uInf = np.array([1, Minf * np.cos(alpha), Minf * np.sin(alpha), (1 / (gamma ** 2 - gamma)) + 0.5 * (Minf ** 2)])
    return uInf


def plotCp(x, cp, savePath=None):
    # cp plot
    f = plt.figure(figsize=(8, 8))
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel(r'$c_p$', fontsize=16)
    for name in WALL_NAMES:
        plt.scatter(x[name], cp[name], label=name)
    f.tight_layout()
    plt.legend()
    savePlot(f, savePath)
    plt.show()
    plt.close(f)


def main():
    datFile = 'dat/rusanov_CFL0.9_secondOrder.dat' # converged solution from first-order FVM
    mesh = 'gri/smoothed_local_all.gri' # mesh file
    u = np.loadtxt(datFile)
    Minf = 0.25 # freestream Mach number - must match main.cpp
    alpha = 8 # angle of attack in degrees
    uInf = computeFreestreamState(Minf, alpha)

    cd, cl, x, cp = getCoefficients(u, uInf, mesh, alpha)
    print(f"cd = {cd}, cl = {cl}")
    plotCp(x, cp, figPath(datFile, mesh, 'cp'))
    plotMach(u, mesh, figPath(datFile, mesh, 'mach'))


if __name__ == "__main__":
    main()

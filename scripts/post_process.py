import argparse
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from readgri import readgri


LIMITER_NAMES = ('NONE', 'BJ', 'LCD')  # must match the limiterType strings main.cpp writes


def limiterTag(datFile):
    """Limiter suffix to put in a figure name, or None for a first-order result.

    main.cpp names second-order solutions
    dat/secondOrder/<flux>_CFL<cfl>_secondOrder_<limiter>.dat, so the limiter is
    normally already the last field of the stem and is returned
    as-is. Files written before that suffix existed are still second-order but
    record nothing about the limiter, and their plots would otherwise be
    indistinguishable from a tagged run - those are marked "limUnknown".
    """
    stem = Path(datFile).stem
    if 'secondOrder' not in stem and '2ndOrder' not in stem:
        return None
    tail = stem.rsplit('_', 1)[-1]
    return tail if tail in LIMITER_NAMES else 'limUnknown'


def figPath(datFile, meshFile, kind):
    """fig/<dat stem>__<mesh stem>_<kind>.png, creating the folder if needed.

    Naming the figure after both inputs means a plot always says which solution
    and which mesh produced it, so runs cannot be mixed up later. For a
    second-order solution that also means the limiter: it is appended only when
    the stem does not already end with it, so a properly named file does not come
    out as "..._BJ_BJ".
    """
    stem = Path(datFile).stem
    tag = limiterTag(datFile)
    if tag is not None and not stem.endswith(f"_{tag}"):
        stem = f"{stem}_{tag}"
    # main.cpp puts the mesh in the solution name too, so appending it again
    # would read "..._smoothed_local_all_LCD__smoothed_local_all_cp.png". Older
    # files written before that still need it, hence the check rather than a
    # plain removal.
    meshStem = Path(meshFile).stem
    name = f"{stem}_{kind}.png" if meshStem in stem else f"{stem}__{meshStem}_{kind}.png"
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


def plotMach(U, mesh, savePath=None, show=True):
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
    if show:
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


def plotCp(x, cp, savePath=None, show=True):
    # cp plot
    f = plt.figure(figsize=(8, 8))
    plt.xlabel('x (m)', fontsize=16)
    plt.ylabel(r'$c_p$', fontsize=16)
    for name in WALL_NAMES:
        plt.scatter(x[name], cp[name], label=name)
    f.tight_layout()
    plt.legend()
    savePlot(f, savePath)
    if show:
        plt.show()
    plt.close(f)


NO_INPUT = ("ERROR: no input available. Pass --dat, --mesh, --minf and --alpha"
            " to run without prompts.")


def chooseFile(folder, extension, prompt, recursive=False):
    """Numbered menu over the files in a folder, returning the one picked.

    Mirrors chooseFileFromFolder in processMesh.h so this reads the same way as
    the solver: identical sorted listing, identical 1-based numbering. Paths are
    relative to the working directory, so this - like the solver - expects to be
    run from the repository root.
    """
    root = Path(folder)
    pattern = '*' + extension
    # dat/ splits into firstOrder/ and secondOrder/, so solutions need a
    # recursive walk while the flat gri/ folder does not.
    names = sorted(root.rglob(pattern) if recursive else root.glob(pattern))
    if not names:
        raise SystemExit(
            f"ERROR: no {extension} files in {root}\n"
            "       Nothing has been written there yet, or you are not running"
            " from the project root."
        )

    print(f"{prompt}:")
    for i, name in enumerate(names, 1):
        print(f"  {i} - {name.as_posix()}")

    while True:
        try:
            answer = input("Enter Option: ").strip()
        except EOFError:
            raise SystemExit(NO_INPUT)
        if answer.isdigit() and 1 <= int(answer) <= len(names):
            return names[int(answer) - 1].as_posix()
        print(f"  Please enter a number from 1 to {len(names)}.")


def meshFromDat(datFile, folder='gri'):
    """Work out which mesh a solution was computed on, from its filename.

    main.cpp names solutions <flux>_CFL<cfl>_<order>_<mesh>[_<limiter>], so the
    mesh is everything between the order field and the optional limiter. Mesh
    names contain underscores themselves ("smoothed_local_all"), so the known
    fields are peeled off the two ends rather than counted from the front.

    Returns None when the name predates that scheme, or names a mesh that is not
    in the folder - the caller then falls back to asking, which is the right
    answer rather than guessing wrong. Plotting a solution against the wrong mesh
    produces a figure that looks plausible and is meaningless.
    """
    root = Path(folder)
    stem = Path(datFile).stem
    fields = stem.split('_')

    if len(fields) > 3 and fields[1].lower().startswith('cfl'):
        middle = fields[3:]                      # everything past <flux>_CFL<cfl>_<order>
        if middle and middle[-1] in LIMITER_NAMES:
            middle = middle[:-1]
        if middle:
            guess = root / ('_'.join(middle) + '.gri')
            if guess.is_file():
                return guess.as_posix()

    # Older names put no mesh in the stem, and hand-renamed files may not follow
    # the field layout at all, so fall back to any mesh whose name appears as a
    # whole field. Longest first, so a short name cannot win over a longer one
    # that contains it.
    for path in sorted(root.glob('*.gri'), key=lambda p: len(p.stem), reverse=True):
        if f"_{path.stem}_" in f"_{stem}_":
            return path.as_posix()
    return None


def promptFloat(prompt, default):
    """Read a number, taking a blank line as the default.

    Minf and the angle of attack are not recorded anywhere in a solution file -
    only the flux and CFL are, through the filename - so they have to come from
    the user. The defaults are the subsonic condition from the project statement.
    """
    while True:
        try:
            answer = input(f"{prompt} [{default}]: ").strip()
        except EOFError:
            raise SystemExit(NO_INPUT)
        if not answer:
            return default
        try:
            return float(answer)
        except ValueError:
            print("  Not a number.")


def main():
    parser = argparse.ArgumentParser(
        description="Plot cp and Mach contours for a solution, and report cl and cd."
    )
    parser.add_argument('--dat', help='solution file to plot (default: pick from a list)')
    parser.add_argument('--mesh',
                        help='mesh the solution was computed on'
                             ' (default: read from the solution name, else pick from a list)')
    parser.add_argument('--minf', type=float, help='freestream Mach number (default: ask)')
    parser.add_argument('--alpha', type=float, help='angle of attack in degrees (default: ask)')
    parser.add_argument('--no-show', action='store_true',
                        help='write the figures without opening a window')
    args = parser.parse_args()

    # Anything not given on the command line is asked for, so one script serves
    # both a quick interactive look and a scripted batch of runs.
    datFile = args.dat or chooseFile('dat', '.dat', 'Choose Solution File (.dat)', recursive=True)

    # The solution name records its own mesh, so the usual case needs no answer
    # here. Where it cannot be read - a file predating that naming - fall back to
    # asking rather than picking something that would plot without complaint.
    mesh = args.mesh or meshFromDat(datFile)
    meshSource = 'from solution name' if (mesh and not args.mesh) else None
    if mesh is None:
        print("Could not tell which mesh " + Path(datFile).name + " was computed on.")
        mesh = chooseFile('gri', '.gri', 'Choose Mesh File (.gri)')

    Minf = args.minf if args.minf is not None else promptFloat('Freestream Mach number', 0.25)
    alpha = args.alpha if args.alpha is not None else promptFloat('Angle of attack (degrees)', 8.0)

    print(f"solution: {datFile}")
    print(f"mesh:     {mesh}" + (f"  ({meshSource})" if meshSource else ""))
    print(f"Minf = {Minf}, alpha = {alpha} deg")

    u = np.loadtxt(datFile)
    uInf = computeFreestreamState(Minf, alpha)

    cd, cl, x, cp = getCoefficients(u, uInf, mesh, alpha)
    print(f"cd = {cd}, cl = {cl}")
    plotCp(x, cp, figPath(datFile, mesh, 'cp'), show=not args.no_show)
    plotMach(u, mesh, figPath(datFile, mesh, 'mach'), show=not args.no_show)


if __name__ == "__main__":
    main()



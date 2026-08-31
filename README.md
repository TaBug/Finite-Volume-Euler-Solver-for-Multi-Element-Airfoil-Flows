# AE 623 Project 2 — 2D Euler Finite Volume Solver

Cell-centered finite volume solver for the 2D Euler equations on unstructured
triangular meshes. First order in space to converge a starting field, then
second order (Green-Gauss reconstruction, optional Barth-Jespersen limiter)
marched to steady state with SSP-RK2.

## Layout

```
main.cpp              driver: reads the mesh, runs 1st order, then 2nd order
solver.h              2nd-order residual, face geometry cache, boundary states
getRes.h  FE_FVM.h    1st-order residual and its forward-Euler driver
RK2_FVM.h             SSP-RK2 time marching
fluxes.h              Roe / Rusanov / HLLE, plus the analytic wall flux
limiters.h            Barth-Jespersen limiter
processMesh.h         .gri reader and mesh topology (I2E, B2E, E2F, normals)
data_conversion.h     state <-> text, filename <-> settings
tools.h               small shared helpers

testers/test_flux.cpp flux function checks

gri/                  meshes the solver reads
msh/                  the gmsh sources the .gri files came from
dat/firstOrder/       converged first-order solutions written by main
dat/secondOrder/      converged second-order solutions written by main
fig/                  plots written by scripts/post_process.py
scripts/              post-processing and mesh plotting (Python)
docs/                 report
legacy/               older results and scripts, kept for reference; nothing builds against them
```

## Building

Requires MSVC and header-only [Eigen 3.4.0](https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip)
extracted to `%USERPROFILE%\eigen-3.4.0`.

```
build_opt.bat main.cpp main.exe      # optimized (/O2, no debug info)
.vscode\build.bat main.cpp main.exe  # same flags plus /Zi, objects into build/
```

In VS Code, `build active file (MSVC)` builds whichever file is open.

## Running

Run from the repository root — the solver looks for meshes in `gri/` and writes
to `dat/`, both relative to the working directory.

```
main.exe
```

It prompts for the flow condition, mesh, flux function, CFL, and limiter. The two
flow conditions are the ones the project statement specifies, both at 8 degrees
angle of attack:

| | Minf | alpha |
|---|---|---|
| subsonic  | 0.25 | 8 deg |
| transonic | 0.50 | 8 deg |

Solutions are named
`dat/<order>/<flux>_CFL<cfl>_<order>_<flow>_<mesh>[_<limiter>].dat`, so a file
records the settings that produced it; `data_conversion.h` reads the flux and CFL
back out of that name when a solution is reloaded, and `scripts/post_process.py`
reads the flow condition and mesh. A state file holds only the four conserved
variables, so the name is the only record of which condition produced it —
reloading a first-order solution written at the other condition prints a warning.

The limiter stays last because `post_process.py` reads it off the final field;
mesh and flow sit before it.

## Post-processing

```
python scripts/post_process.py       # cp plot, Mach contours, cl and cd
```

It lists the solutions in `dat/` and asks which one to plot. The mesh and the flow
condition are both read out of the solution's own name, so nothing else normally
has to be answered; it falls back to asking only when the name predates that
scheme.

Every answer can be given up front instead, which is also how to run it
unattended:

```
python scripts/post_process.py --dat dat/secondOrder/rusanov_CFL0.9_secondOrder_c0_BJ.dat \
                               --minf 0.25 --alpha 8 --no-show
```

`--mesh` overrides the mesh read from the name. Figures land in `fig/`, named
after the solution and the mesh.

## Checks

```
.vscode\build.bat testers\test_flux.cpp test_flux.exe
test_flux.exe
```

Evaluates all three fluxes across a boundary at a 30 degree angle of attack and
prints them for comparison against the hand calculation in the report.

The old `verify_solver` and `profile_*` harnesses are gone. They existed to check
and time the optimized residual against `secondOrderFV`, a nested-vector
reference implementation of the same scheme; that reference has been removed, so
they had nothing left to compare against. `git log` has them if the comparison is
ever needed again.

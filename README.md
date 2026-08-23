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

verify_solver.cpp     optimized path vs. reference implementation, plus timings
test_flux.cpp         flux function checks
profile_*.cpp         profiling harnesses for the residual and allocations

gri/                  meshes the solver reads
msh/                  the gmsh sources the .gri files came from
dat/                  converged solutions written by main
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

It prompts for the flux function, CFL, and limiter. Solutions are named
`dat/<flux>_CFL<cfl>_<order>[_<limiter>].dat`, so a file records the settings
that produced it; `data_conversion.h` reads the flux and CFL back out of that
name when a solution is reloaded.

## Post-processing

```
python scripts/post_process.py       # cp plot, Mach contours, cl and cd
```

Figures land in `fig/`, named after both the solution and the mesh. Edit
`datFile` and `mesh` at the bottom of the script to choose a run.

## Checks

```
build_opt.bat verify_solver.cpp verify_solver.exe
verify_solver.exe c1.gri 20
```

Marches the same initial condition through the reference implementation and the
optimized path and reports the largest difference between them, which should be
zero, along with the speedup.

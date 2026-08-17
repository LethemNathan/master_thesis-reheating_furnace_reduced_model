# Reheating Furnace – Radiation & Conduction Model

This project solves the coupled radiation, wall conduction, slab conduction,
and gas energy-balance problem inside a reheating furnace.

The code is a **pipeline**: the scripts must be run in a fixed order because
each script uses variables created by the setup script.

## Pipeline overview

```
full_Model_Setup.m   (Part I – setup)
        |
        |--->  DOM_Coupled_Solution.m                (Option A: full DOM)
        |
        '--->  direct_Exchange_Area_Computation.m
                        |
                        '--->  full_Reduced_Model_Coupling.m   (Option B: reduced model)
```

## How to run

You must always run the setup script first. Then choose **one** of the two options.

### Step 1 – Setup (always required)

Run `full_Model_Setup.m`.

This script loads the mesh, builds the geometry and topology, sets the material
properties and boundary conditions, and creates all the data needed by the next
scripts (geometry, mesh, radiation surfaces, contexts, initial temperatures,
gas energy model). Do not clear the workspace after this step.

### Step 2 – Choose one solver

**Option A – Full DOM solution**

Run `DOM_Coupled_Solution.m`.

This solves the coupled problem using the full Discrete Ordinates Method (DOM).
It iterates between radiation, wall conduction, and slab conduction until the
temperatures converge.

**Option B – Reduced model**

1. Run `direct_Exchange_Area_Computation.m`.
   This builds the direct exchange-area matrix `S` from the DOM.
2. Run `full_Reduced_Model_Coupling.m`.
   This solves the coupled problem using the reduced radiation model, plus wall
   conduction, slab conduction, and the gas-zone energy balance.

## Important notes

- Run the scripts in the **same MATLAB session**. Later scripts reuse variables
  from the setup script and do not create them again.
- Option B needs `S` and `exchOpts`, which come from
  `direct_Exchange_Area_Computation.m`. Run that script before
  `full_Reduced_Model_Coupling.m`.
- Wall conduction is turned off by default (`solveWalls = false` in the setup).
  Set it to `true` if you want to solve the walls.

## Input files

All input files are expected in the `input/` folder.

### Furnace mesh

The main furnace mesh `input/gmsh_reheating_furnace.msh` must first be generated
using **Gmsh** from the geometry and meshing script provided in the repository:

`input/gmsh_oneBlock_reheating_furnace.geo`

Open this `.geo` file in Gmsh, generate the mesh, and save/export the resulting
mesh as:

`input/gmsh_reheating_furnace.msh`

The MATLAB setup script expects the generated mesh to be available at this
location.

The second mesh file is provided separately:

- `input/nonParticipatingMedia.msh` – markers for non-participating media

### Definition scripts (run inside the setup)

- `input/subZonesDefinition.m`
- `input/bramestoreSlab.m`
- `input/physicalConstant.m`
- `input/wallDataFurnace_Generation.m`

### Data files (.mat)

- `input/postprocZones.mat` – zone temperatures
- `input/extWallReport.mat` – external wall report
- `input/skid.mat` – skid data
- `input/intSurfReport.mat` – inter-zone flow data

## Main parameters (set in `full_Model_Setup.m`)

- Angular quadrature: `Nth = 4`, `Nph = 8`
- Initial gas temperature: `T_init = 1300` K
- Solve walls: `solveWalls = false`
- Gas absorption coefficient: `kappa = 0.21`
- Imposed heat flux: `q_imposed = 30500` W/m²
- Coupled loop: `maxOuterIter = 30`, `outerTol = 1e-3`

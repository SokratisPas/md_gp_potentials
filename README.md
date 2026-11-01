# Molecular Dynamics with a Lennard-Jones Potential

This repository contains a C++ molecular-dynamics simulation. It's a cleaner and improved version of the https://github.com/SokratisPas/Melt-Simulation" simulation without the rendering part. 
It reads an extended XYZ structure, initializes particle velocities, and advances the system using velocity Verlet integration, periodic boundary conditions, a neighbor list, and a Berendsen thermostat.

## Build and run

You need a C++17-compatible compiler. From the repository root, build and run:

```sh
g++ -std=c++17 src/main.cpp -o md-gp-pt
./md-gp-pt
```

The program reads `input/YOUR_DATA_FILE.txt` and writes its results to
`output/lj.dat`.

## Input format

The input is an extended XYZ file. See `input/dataExample.txt`.


The current implementation uses the three diagonal lattice values as the
orthogonal box dimensions.

## Output

`output/lj.dat` is comma-separated and contains a header followed by sampled
simulation data: step, time, kinetic energy, potential energy, total energy,
temperature, and pressure.

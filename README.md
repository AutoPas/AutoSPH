# AutoSPH

An SPH simulator built on top of the [AutoPas](https://github.com/AutoPas/AutoPas) library.

This project was developed for the Master's Thesis in Computational Science and Engineering: **Implementation of an Auto-Tuned SPH Simulator with AutoPas** by Iñigo Ruiz, advised by Manish Mishra at TUM. The thesis document is available [here](https://mediatum.ub.tum.de/1864913).

It is intended for 3D simulaton of incompressible flows.
The `input` folder contain some examples of initial configuration files.
Below are shown animations and screenshots of the output visualized with [ParaView](https://www.paraview.org/).

## Benchmark cases

### Dam Break

![dambreak_benchmark](https://github.com/AutoPas/AutoSPH/blob/main/docs/animations/dam_break.gif?raw=true)

### Sloshing

![sloshing_benchmark](https://github.com/AutoPas/AutoSPH/blob/main/docs/animations/sloshing.gif?raw=true)


## Paraview Analysis

![dambreak_paraview](https://github.com/AutoPas/AutoSPH/blob/main/docs/paraview_screenshots/paraview_dambreak.png?raw=true)

![sloshing_paraview](https://github.com/AutoPas/AutoSPH/blob/main/docs/paraview_screenshots/paraview_sloshing.png?raw=true)


## Variable Number of Particles

### 1728 Particles

![dambreak_1728](https://github.com/AutoPas/AutoSPH/blob/main/docs/animations/db_h12n06.gif?raw=true)

### 32768 Particles

![dambreak_32768](https://github.com/AutoPas/AutoSPH/blob/main/docs/animations/db_h12n16.gif?raw=true)


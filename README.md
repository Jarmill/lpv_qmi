# lpv_qmi
Data-Driven Control of Linear Parameter Varying (LPV) systems under Process Noise with Quadratic Matrix Inequalities (QMI).

Continuous-time systems supply noisy state-derivative-parameter observations. Discrete-time process require noisy state-parameter-next state observations.

A gain scheduling controller is synthesized by enforcing a QMI at each vertex of the parameter polytope.

The code is currently constructed for elementwise-L2 noise, but other types of noise models (energy bounds) may be enforced by modification of the QMIs.


## Instructions
Generate a trajectory (with elementwise-L2 noise of bound epsilon) using `lpvsim.sim`. Define the vertices of the polytope (such as a manual definition or using lcon2vert from https://www.mathworks.com/matlabcentral/fileexchange/30892-analyze-n-dimensional-convex-polyhedra). Pass the trajectory into the `lpvstab` object for discrete-time stabilization, and then attempt generation of a gain scheduled controller on these vertices with the method `lpvstab.stab`. 

For continuous-time systems use `lpvstab_cont` rather than `lpvstab`.

For H2 control in discrete-time, use `lpvh2` rather than `lpvstab`.

## Dependencies

- YALMIP: https://yalmip.github.io/
- Mosek: https://www.mosek.com/ (or any solver compatible with YALMIP)

All code is written and tested on Matlab R2021a.

## Reference

https://ieeexplore.ieee.org/abstract/document/9971732



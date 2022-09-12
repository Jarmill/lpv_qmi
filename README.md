# lpv_qmi
Data-Driven Control of Linear Parameter Varying (LPV) systems under Process Noise with Quadratic Matrix Inequalities (QMI).

Continuous-time systems supply noisy state-derivative-parameter observations. Discrete-time process require noisy state-parameter-next state observations.

A gain scheduling controller is synthesized by enforcing a QMI at each vertex of the parameter polytope.

The code is currently constructed for elementwise-L2 noise, but other types of noise models (energy bounds) may be enforced by modification of the QMIs.

## Dependencies

- YALMIP: https://yalmip.github.io/
- Mosek: https://www.mosek.com/ (or any solver compatible with YALMIP)

All code is written and tested on Matlab R2021a.



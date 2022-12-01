# Reservoir_Simulator

To execute file, please run Main.m file
All the reservoir data can be found under the same folder


Description:
Development of a three-phase two-dimensional reservoir model

Input Data Routine
• Developed a routine to read input data files. The routine is coded as general as possible to irregular grid inputting with heterogeneous and anisotropic property distribution. For uniformly distributed properties, the subroutine has the option for reading only one entry and assigning it over the entire domain. Also a checking mechanism is provided for improper input distribution.

Reservoir Model
• It can run on 1-D and 2-D (with elevation change) geometries for single-, two- and three-phase simulation studies.
• It includes a subroutine for handling the PVT, relative permeability and capillary pressure data. Uses a table look-up technique with proper interpolation. Implement “the out-of-the range of the table” checks.
• It incorporates the variable bubble point formulation.
• It constructs the Jacobian matrix and the right-hand-side vector for each Newton-Raphson iteration. Uses numerical differentiation to calculate entries of the Jacobian. For the relative permeability calculations, it allows for single- or two-point upstream weighting.
• It solves the resulting Jacobian matrix together with the right-hand-side vector using any [A][X]=[B] type solver.

For more details please see the Final_Report.pdf in the repository.

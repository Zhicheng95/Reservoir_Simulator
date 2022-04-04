# Reservoir_Simulator

Development of a three-phase two-dimensional reservoir model
General Guidelines
Input Data Routine
• Develop a routine to read input data files. The routine should be coded as general as possible to irregular grid inputting with heterogeneous and anisotropic property distribution. For uniformly distributed properties, the subroutine should have the option for reading only one entry and assigning it over the entire domain. Provide a checking mechanism for improper input distribution (required).
• Provide options to echo print or not to print data read (required)

Reservoir Model
• Design your simulator in such a way that it can run on 1-D and 2-D geometries for single-, two- and three-phase simulation studies (optional).
• All the constant parts of all coefficients should be calculated at the beginning of the simulation (required).
• Develop the subroutine for handling the PVT, relative permeability and capillary pressure data. Use a table look-up technique with proper interpolation. Implement “the out-of-the range of the table” checks (required).
• Incorporate the variable bubble point formulation (required).
• Construct the Jacobian matrix and the right-hand-side vector for each Newton-Raphson iteration. Use numerical differentiation to calculate entries of the Jacobian. For the relative permeability calculations, allow for single- or two-point upstream weighting (required).
• Solve the resulting Jacobian matrix together with the right-hand-side vector using any [A][X]=[B] type solver (required).

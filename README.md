# Gauss-Seidel

Function to solve a linear system of equations, defined by a square matrix `A` and a vector `B`, using the Gauss-Seidel method

##### Outputs
`x` - n x 1 array - solutions of the linear system, where n is the number of solutions

##### Inputs
`A`   - n x n array - matrix containing coefficients of the variables for each equation of the linear system

`B`   - n x 1 array - vector containing the constant of each equation of the system

`x_0` - 1 x n array - array containing initial guesses for each of the unkowns, used for the Gauss-Seidel method (usually filled with zeros)

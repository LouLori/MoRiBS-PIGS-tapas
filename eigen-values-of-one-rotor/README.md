<<<<<<< eb5a87762ed599c52ca6527acaec59ec6259cb15
# eigenvalues-of-two-linear-rotors
=======
# DVR-for-spherical-harmonics
Here I have several codes:
1. rs.f 
It gives us eigenvalues and eigenfunctions of a real hermitian matrix;

2. normalization-test-for-spherical-harmonics.f
Here we can check the normalization condition for spherical harmonics by creating <j1,m1|j2,m2> matrix, where j1,j2 are the angular momentum quantum numbers and m1,m2 are the azimuthal quantum numbers.

3. gauss-legendre-quadrature-by-matrix-diagonalization.f
Here we created gauss-legendre quadrature points by diagonalizing <j1|x|j2> matrix and tested \sum_{j} P_{l1}(x_{j})P_{l2}(x_{j}) =delta_{l1,l2}

4. gauss-quadrature-for-legendre-polynomials-and-phi.f
Here I have used gauleg.for subroutines (taken from Numerical recipes book) and Eqs. (A10a - A10b in J. Chem. Phys. vol. 96, page - 1982, year 1992) and Eq. 16 (J. chem. Phys. vol 102, page 3622, year 1995) to generate DVR points and weights for angular and azimuthal quantum numbers, respectively.

5. vh2h2.f
This is the potential routines for two linea rotors interactions. all-coefs, short-range and long-range files are needed for initialization of the vh2h2.f.

6. associated-legendre.cc call-in-c++-legendre.f
is c++ code for associated legendre polynomial where call-in-c++-legendre.f is called. By the way, it will not give correct results. 

7. script.sh
it is for compiler of 6 items
>>>>>>> eigen final

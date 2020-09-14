# Matrix-Class
Matrix Class for C++ including stable and smart inverse.

Efficient implementation of a C++ matrix (and vector) class.
The implementation provides the following auxiliary operations:

\+ adds matrices or vectors,

\+ adds a scalars to a vector,

\+ adds a scalar on the matrix diagonal,

\- subtracts matrices or vectors

\- subtracts a scalar from a vector

\- subtracts the identity of a matrix

\* multiply matrices

\\ solve a linear system: This major implementation employs Householder reflections to find the solution of the linear system A\*x= b.
The implementation works for rank deficient and rectangular matrices. For regular matrices, the solution is x= A<sup>-1</sup> b; for singular matrices, x= A<sup>\+</sup>b, where A<sup>\+</sup> is the (possibly rank deficient) pseudoinverse. The decomposition can be reused (recycled).

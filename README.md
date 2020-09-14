# Matrix-Class
Matrix Class for C++ including a stable and smart inverse.

Efficient implementation of a C++ matrix (and vector) class.

The implementation provides the following auxiliary operators:

\+ adds matrices or vectors,

\+ adds a scalars to a vector,

\+ adds the identity matrix,

\- subtracts matrices or vectors,

\- subtracts a scalar from a vector,

\- subtracts the identity of a matrix,

\* multiply matrices,

\\ solves a linear system: this major implementation employs Householder reflections to efficiently find the solution of the linear system A\*x= b.
The implementation works for rank deficient and rectangular matrices.
For regular matrices, the solution is x= A<sup>-1</sup> b; for singular matrices, x= A<sup>\+</sup>b, where A<sup>\+</sup> is the (possibly rank deficient) pseudoinverse. The decomposition can be reused (recycled).
The complete orthogonal decomposition pivots the rows, cf. [notes](https://www.tu-chemnitz.de/mathematik/fima/public/mathematischeStatistik.pdf#page=113).

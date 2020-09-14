# Matrix-Class
Matrix Class for C++ including a stable and smart inverse.

Efficient implementation of a C++ matrix (and vector) class. 
The implementation provides the following operators on matrices and vectors:

\+ adds matrices or vectors,

\+ adds a scalar to a vector,

\+ adds the identity matrix,

\- subtracts matrices or vectors,

\- subtracts a scalar from a vector,

\- subtracts the identity from the matrix,

\* multiply matrices,

\/ solves a linear system: this main implementation of the class employs Householder reflections (one of the top 10 algorithms of the 20th century) to efficiently find the solution of the linear system A\*x= b.
The rank-revealing implementation works for rank deficient and rectangular matrices.
For regular matrices, the solution is x= A<sup>-1</sup> b; for singular matrices, x= A<sup>\+</sup>b, where A<sup>\+</sup> is the generalized inverse of the possibly rank deficient matrix A. The second decomposition is efficient by exploiting the structure of the resulting matrix. The decomposition can be reused (recycled).
The complete orthogonal decomposition pivots the rows, cf. [notes](https://www.tu-chemnitz.de/mathematik/fima/public/mathematischeStatistik.pdf#page=113).

Invoke the library by calling
```cpp
#include "matrix.h"
```

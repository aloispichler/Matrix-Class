#ifndef matrix_H		//	© 2020, Alois Pichler
#define matrix_H


#include <cstring>		// enable memcpy
#include <cassert>		// enable assert


class Matrix;			// forward declaration (see below)


template<typename T= double>				// the default type is double
class Vector			// declaration of the class Vector
{	public:
	Vector<T>(unsigned _dim= 0)				// default constructor
	{	vLen= _dim; vData= new T[vLen];}
	Vector<T>(const Vector<T>&);			// copy constructor
	~Vector<T>(){	delete [] vData;}		// destructor

	template<unsigned dim>					// construct from array
	Vector<T>(const T (&vec)[dim]) : Vector<T>(dim)	// v inherits from Vector
	{	if (vData != vec)	// Vector<double> vec({ -1, 7, 1});
			std::memcpy(vData, vec, dim*sizeof(T));}

	Vector<T>(const Matrix&){};				// construct from Matrix

	Vector<T>& operator= (const Vector<T>&);// assignment operator=
  	T& operator[] (unsigned) const;			// element access operator
	T& operator() (unsigned) const;			// the indexes are one-based, not zero based.

	unsigned Length() const {return vLen;}	// count elements
	Vector<T> Zeros();						// fill Vector with 0
	Vector<T> operator+ (T);				// add a scalar
	Vector<T> operator+ (const Vector<T>&);	// Vector addition
	Vector<T> operator- (const Vector<T>&);	// Vector subtraction
	Vector<T> operator* (T);				// multiplication by scalar
	unsigned maxPosition();					// maximum [] position

	private:
		T* vData; unsigned vLen;
};


class Matrix
{	public:

	Matrix(): mRows(0), mCols(0), mData(nullptr) {};	// default constructor
	Matrix(const unsigned);					// constructor of a square matrix
	Matrix(const unsigned, const unsigned);	// constructor of a rectangular matrix
	Matrix Let(std::initializer_list<double>);
	Matrix(const Matrix&);					// copy constructor
	~Matrix(){	delete [] mData;}			// destructor
	Matrix(const Vector<double>&);			// construct from Vector


	template<unsigned mRows, unsigned mCols>	// construct from array
	Matrix(const double (&arr)[mRows][mCols]) : Matrix(mRows, mCols)
	{	std::memcpy(mData, arr, mRows* mCols* sizeof(double));}
	//	Matrix A((double[][3]){{ -1.6, 7, 13}, {5, 9.3, 2}}); enable array initialization

	Matrix& operator= (const Matrix&);		// assignment operator=

	double* operator[] (unsigned) const;	// element access operator
	double& operator()(unsigned, unsigned) const;//	the indexes are ONE-based, NOT zero based.
	unsigned rows() const {return mRows;}	// the number of rows
	unsigned cols() const {return mCols;}	// the number of columns

	Matrix operator+ (const Matrix&);		// matrix addition
	Matrix operator+ (double);				// add a multiple of the unit matrix
	Matrix operator- (const Matrix&);		// matrix subtraction
	Matrix operator* (double);				// multiplication by scalar
	Matrix operator* (const Matrix&);		// matrix multiplication
	Matrix operator/ (Matrix);				// Solve A/b
	Vector<double> operator/ (Vector<double>);	// Solve A/b

	Matrix Fill(const double);
	
	private:
	double *mData; unsigned mRows, mCols;
};


static Matrix Eye(const unsigned);			// identity matrix
Matrix Transpose(Matrix);					// transpose matrix
Matrix Inverse(Matrix, double);				// pseudoinverse
std::ostream& operator << (std::ostream&, const Matrix);


struct QResult	// container to hold result of QR decomposition
{	Matrix Mat;		// the matrix Mat holds Q and R
	unsigned Rank;	// the rank of the matrix
	Vector<double>    V1, V2;	// auxiliary vectors for the decomposition
	Vector<unsigned>  Permutation;		// permutation to pivot rows
};


template<typename T>	//	C++ lacks Square
const T Square(const T& a)
{	return a* a;}


#include "vector.cpp"	// detailed definitions
#include "inverse.cpp"
#include "matrix.cpp"


#endif		// matrix_H
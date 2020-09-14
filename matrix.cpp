#ifndef matrix_CPP			//	© 2020, Alois Pichler
#define matrix_CPP


#include <iomanip>			// std::cout, std::endl
#include "matrix.h"			// header file


// class MatrixException
// {	public:
// 		std::string msg;
// 		MatrixException(std::string arg) : msg(arg) {}
// };


Matrix::Matrix(unsigned dim)					//	Construct a square matrix
{	mRows= dim; mCols= dim; mData= new double [dim* dim];}


Matrix::Matrix(unsigned dimR, unsigned dimC)	// construct a rectangular matrix
{	mRows= dimR; mCols= dimC; mData= new double [dimR* dimC];}
	
	
Matrix::Matrix(const Matrix& src)	// copy constructor
{	mRows= src.mRows; mCols= src.mCols;
	mData= new double[mRows* mCols];		// deep copy
	std::memcpy(mData, src.mData, mRows* mCols* sizeof(double));}


Matrix::Matrix(const Vector<double>& src)	// construct from Vector
{	mRows= src.Length(); mCols= 1;
	mData= new double[mRows];	// deep copy from Vector source
	std::memcpy(mData, &src[0], mRows* sizeof(double));}


Matrix& Matrix::operator=(const Matrix& src)	// assignment operator=
{	if (this != &src)		// protect against invalid self-assignment
	{	if (mRows* mCols != src.mRows* src.mCols)
		{	delete[] mData;			// destroy old stack
			mData= new double[src.mRows* src.mCols];}
		mRows= src.mRows; mCols= src.mCols;	// deep copy
		std::memcpy(mData, src.mData, mRows* mCols* sizeof(double));}
	return *this;}


double* Matrix::operator[](const unsigned _row) const 	// row access operator
{	assert((_row < mRows) && "Matrix: Row out of range.");
	return & mData[mCols* _row];}		//	row-major


double &Matrix::operator()(const unsigned row, const unsigned col) const 	//	element access operator
{	assert((row > 0 && row <= mRows) && "Matrix: row out of range.");
	assert((col > 0 && col <= mCols) && "Matrix: col out of range.");
	return mData[(row-1)*mCols+ col-1];}	//	the indexes are one-based, not zero based.


Matrix Matrix::Fill(const double val= 0)	 //	fill this matrix with entries val
{	for (unsigned i= 0; i< mRows; i++)
		for (unsigned j= 0; j< mCols; j++)
			(*this)[i][j]= val;
	return *this;}


Matrix Eye(const unsigned mRows)	// provides the identity matrix
{	Matrix c(mRows, mRows);
	for(unsigned i= 0; i< mRows; i++)
	{	for(unsigned j= 0; j< mRows; j++) c[i][j]= 0;
		c[i][i]= 1;}
	return c;}


Matrix Matrix::operator+ (double lambda)	// add lambda* identity
{	Matrix c(mRows, mCols);
	unsigned k= 0;
	for	(unsigned i= 0; i< mRows; i++)
		for (unsigned j= 0; j< mCols; j++)
		{	c.mData[k]= this->mData[k]; if (i==j) c.mData[k]+= lambda;
			k++;}
	return c;}


Matrix Matrix::operator+ (const Matrix& other)	// add matrices
{	assert((other.rows()== mRows && other.cols()== mCols) && "Matrices can't be added.");
	Matrix c(mRows, mCols);
	for	(unsigned i= 0; i< mRows* mCols; i++)
		c.mData[i]= (*this).mData[i] + other.mData[i];
	return c;}


Matrix Matrix::operator- (const Matrix& other)	// subtract matrices
{	assert((other.rows()== mRows && other.cols()== mCols) && "Matrices cannot be substracted.");
	Matrix diff(mRows, mCols);
	for (unsigned i= 0; i< mRows* mCols; i++)
		diff.mData[i]= (*this).mData[i] - other.mData[i];
	return diff;}


Matrix Matrix::operator* (double lambda)	// multiplication by scalar
{	Matrix c(mRows, mCols);
	for (unsigned i= 0; i< mRows; i++)
		for (unsigned j= 0; j< mCols; j++)
			c.mData[i* mCols+ j]= lambda* this->mData[i* mCols+ j];
	return c;}


Matrix Transpose(Matrix A)					// transposition
{	Matrix ATranspose(A.cols(), A.rows());
	for (unsigned i= 0; i< ATranspose.rows(); i++)
		for (unsigned j= 0; j< ATranspose.cols(); j++)
			ATranspose[i][j]= A[j][i];
	return ATranspose;}


Matrix Matrix::operator* (const Matrix& other)	// multiply this matrix and another
{	assert((mCols== other.rows()) && "Matrix multiplication: dimension mismatch.");
	Matrix c(mRows, other.cols());
	for (unsigned i= 0; i< mRows; i++)
		for (unsigned j= 0; j< other.cols(); j++)
		{	double *s= &c[i][j]; *s= 0;
			for (unsigned k= 0; k< mCols; k++)
				*s+= (*this).mData[i* mCols+ k] * other.mData[k* other.mCols +j];}
	return c;}


Matrix Matrix::operator/ (Matrix b)	// Solve A* x= b
{	return SolveQRQ(decomposeQRQ(*this), b);}


Vector<double> Matrix::operator/ (Vector<double> b)	// find A/ b
{	return Vector<double>(SolveQRQ(decomposeQRQ(*this), Matrix(b)));}


double Frobenius(const Matrix& mat)	// weighted Frobenius norm
{	double *tmp= mat[0]+ mat.rows()* mat.cols();
	double tmpNorm= 0;				// Hilbert-Schmidt
	for (double* s= mat[0]; s< tmp; s++)
		tmpNorm+= Square(*s);
	return pow(tmpNorm/ mat.rows()/ mat.cols(), 1/ 2);}


std::ostream& operator << (std::ostream& os, const Matrix m)	// output
{	os << "(" << m.rows() << "x" << m.cols() << " matrix):";
	for (unsigned i= 0; i< m.rows(); i++)
	{	os << std::endl;		// new row
		for (unsigned j= 0; j< m.cols(); j++)
			os << std::setw(11) << std::right << m[i][j] << " " ;}
	return os;}


// void Matrix::input()			//	manual input of matrix via console
// {	std::cout << "Enter the elements of the " << mRows << " x " << mCols << " Matrix: " << std::endl;
// 	for(unsigned i=0; i< mRows ; i++)
// 	{	std::cout << "row " << i+1 << ": ";
// 		for(unsigned j= 0; j< mCols; j++)
// 			std::cin >> mData[i* mCols+ j];}}


// Matrix input()					//	manual input of matrix via console
// {	unsigned	mRows, mCols;
// 	std::cout << std::endl << "Matrix input. Enter the dimension of Matrix: "; std::cin >> mRows >> mCols;
// 	Matrix c(mRows, mCols); c.input(); return c;}


#endif	// matrix_CPP
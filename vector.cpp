#ifndef vector_CPP			//	© 2020, Alois Pichler
#define vector_CPP


//#include "matrix.h"		// header file
#include <iostream>


template<typename T>		// copy constructor
Vector<T>::Vector(const Vector<T>& src)
{	vLen= src.vLen;
	vData= new T[vLen];	// deep copy
	std::memcpy(vData, src.vData, vLen*sizeof(T));}


template<>				// construct a Vector from Matrix
Vector<double>::Vector<double>(const Matrix& src)
{	vLen= src.rows()* src.cols();
	vData= new double[vLen];	// deep copy from Vector source
	std::memcpy(vData, src[0], vLen* sizeof(double));}


template<typename T>		// assignment operator=
Vector<T>& Vector<T>::operator=(const Vector<T>& src)
{	if (this != &src)		// protect against invalid self-assignment
	{	if (vLen != src.vLen)
		{	delete[] vData;				// destroy old stack?
			vLen= src.vLen;
			vData = new T[src.vLen];}	// deep copy
		std::memcpy(vData, src.vData, vLen*sizeof(T));}
	return *this;}


template<typename T> 		// element access operator
T& Vector<T>::operator[](const unsigned _index) const
{	return vData[_index];}


template<typename T>		// element access operator
T &Vector<T>::operator()(const unsigned index) const
{	return vData[index- 1];}// the indexes are one-based, not zero based.


template<typename T>
Vector<T> Vector<T>::Zeros()//	set all entries to 0
{	for(unsigned i= 0; i< vLen; i++)
		(*this)[i]= (T) 0;
	return *this;}


template<typename T>		// find the max [] position in the vector
unsigned Vector<T>::maxPosition()
{	unsigned j= 0;
	for (unsigned i= 1; i< vLen; i++)
		if (vData[i]> vData[j]) j= i;
	return j;}


template<typename T>
Vector<T> Vector<T>::operator+ (T lambda)	// add a scalar to this vector
{	Vector<T> c(vLen);
	for	(unsigned i= 0; i< vLen; i++)
		c.vData[i]= this->vData[i]+ lambda;
 	return c;}


template<typename T>
Vector<T> Vector<T>::operator+ (const Vector<T>& other)	// add Vectors
{	assert((vLen== other.Length()) && "Vector length mismatch. ");
	Vector c(vLen);
	for	(unsigned i= 0; i< vLen; i++)
		c.vData[i]= this->vData[i] + other.vData[i];
	return c;}


template<typename T>
Vector<T> Vector<T>::operator- (const Vector<T>& other)	// subtract Vectors
{	assert((vLen== other.Length()) && "Vector length mismatch. ");
	Vector c(vLen);
	for(unsigned i= 0; i< vLen; i++)
		c.vData[i]= (*this).vData[i] - other.vData[i];
	return c;}


template<typename T>
Vector<T> Vector<T>::operator* (T lambda)	// multiply this vector by scalar
{	Vector c(vLen);
	for(unsigned i= 0; i< vLen; i++)
		c.vData[i]= lambda* this->vData[i];
	return c;}


double Norm(const Vector<double> vec, const double p= 2)	// p-norm of a Vector
{	if (p < 1) std::cout << "Norm: 1 > p= " << p << std::endl;
	double tmpNorm= 0;
	if (p == 1)			//	1-norm, weighted
	{	for (unsigned i= 0; i< vec.Length(); i++)
			tmpNorm+= fabs(vec[i]);
		return tmpNorm/ vec.Length();}
	else if (p == 2)	//	2-norm, weighted
	{	for (unsigned i= 0; i< vec.Length(); i++)
			tmpNorm+= Square(vec[i]);
		return sqrt(tmpNorm/ vec.Length());}
	else if (isinf(p))	//	∞-norm
	{	for (unsigned i= 0; i< vec.Length(); i++)
			tmpNorm= std::max(tmpNorm, fabs(vec[i]));
		return tmpNorm;}
	else				//	all other weighted p-norms
	{	for (unsigned i= 0; i< vec.Length(); i++)
			if (vec[i] != 0.) tmpNorm+= pow(fabs(vec[i]), p);
		return pow(tmpNorm/ vec.Length(), 1/ p);}}


template<typename T>		// print
std::ostream& operator<< (std::ostream& os, const Vector<T> &src)
{	os << "(" << src.Length() << "-vector):";
	for(unsigned i= 0; i< src.Length(); i++)
		os << " " << src[i];
	return os;}

#endif	// vector_CPP
/*	This algorithm provides a complete orthogonal decomposition by pivoting the rows,
	see https://www.tu-chemnitz.de/mathematik/fima/public/mathematischeStatistik.pdf#page=113
	For column pivoting and complete orthogonal decomposition
	cf. https://www.math.usm.edu/lambers/mat610/sum10/lecture11.pdf
*/


//	Decompose the Matrix A so that A= P* Q2* L* Q1.
//	The decomposition can be recycled.

QResult decomposeQRQ(Matrix &A, double aTol= 1e-6)	//	QRQ 
{	QResult QRQ; QRQ.Mat= A;			// in-place decomposition. check: Is this A?
	unsigned maxRank= std::min(QRQ.Mat.cols(), QRQ.Mat.rows());	// maximal rank;
	QRQ.Rank= 0;
	QRQ.Permutation= Vector<unsigned>(maxRank);
	for (unsigned k= 0; k< QRQ.Permutation.Length(); k++)
		QRQ.Permutation[k]= k;				// initialize P
	QRQ.V1= Vector<double>(QRQ.Mat.rows());	// auxiliary vector

	Vector<double> normRow2(QRQ.Mat.rows());	// squares of the row norms
	for (unsigned i= 0; i< QRQ.Mat.rows(); i++)	// compute the row norms
	{	normRow2[i]= 0.0;
		for (unsigned j= 0; j< QRQ.Mat.cols(); j++)// compute the norms
			normRow2[i]+= Square(QRQ.Mat[i][j]);}

	for (unsigned k= 0; k< maxRank; k++)	// Householder iteration step
	{	unsigned rMax= k;					// Pivot: find largest row
		for (unsigned i= k+1; i< QRQ.Mat.rows(); i++)
			if (normRow2[i]> normRow2[rMax]) rMax= i;

		if (k != rMax)		//	row pivoting: swap rows
		{	QRQ.Permutation[k]= rMax;
			std::swap(normRow2[k], normRow2[rMax]);
			for (unsigned j= 0; j< QRQ.Mat.cols(); j++)
	 			std::swap(QRQ.Mat[k][j], QRQ.Mat[rMax][j]);}
		if (normRow2[k] < aTol* aTol)	// numerical check
		{	std::cout << "QR: Rank deficient " << A.rows() << "x" << QRQ.Mat.cols() << " matrix. Rank: " << QRQ.Rank << std::endl;
			break;}
		++QRQ.Rank;	// increase the matrix rank

		double normTmp= sqrt(normRow2[k]);		// start Householder Q1
		if (QRQ.Mat[k][k] > 0)
		{	QRQ.V1[k]= QRQ.Mat[k][k]+ normTmp; QRQ.Mat[k][k]= -normTmp;
			normTmp= sqrt(2* normTmp* QRQ.V1[k]);}
		else	// norm of Householder vector
		{	QRQ.V1[k]= QRQ.Mat[k][k]- normTmp; QRQ.Mat[k][k]= normTmp;
			normTmp= -sqrt(-2* normTmp* QRQ.V1[k]);}
		QRQ.V1[k]/= normTmp;		// normalize with V1 > 0
		for (unsigned j= k+1; j< QRQ.Mat.cols(); j++)
			QRQ.Mat[k][j]/= normTmp;// normalize the reflector

		for (unsigned i= k+1; i< QRQ.Mat.rows(); i++)	// Householder reflection-----
		{	double sum= QRQ.V1[k]* QRQ.Mat[i][k];		// compute the inner product
			for (unsigned j= k+1; j< QRQ.Mat.cols(); j++)
				sum+= QRQ.Mat[k][j]* QRQ.Mat[i][j];

			sum*= 2; QRQ.Mat[i][k] -= sum* QRQ.V1[k];	// apply reflector
			for (unsigned j= k+1; j< QRQ.Mat.cols(); j++)
				QRQ.Mat[i][j] -= sum* QRQ.Mat[k][j];
			normRow2[i]-= Square(QRQ.Mat[i][k]);}	// update all norms
	}

	unsigned Delta= QRQ.Mat.rows()- QRQ.Rank;	// rank deficiency
	if (Delta)		// apply smart Householder
	{	QRQ.V2= Vector<double>(QRQ.Rank);			// auxiliary vector
		for (unsigned k= QRQ.Rank- 1; k< -1 ; --k)	// start Householder Q2
		{	double normTmp= 0;
			for (unsigned i= 0; i<= Delta; ++i)
				normTmp+= Square(QRQ.Mat[k+i][k]);
			normTmp= sqrt(normTmp); double *s= &QRQ.Mat[k+Delta][k];
			if (*s > 0)			// sign of Householder transformation
			{	QRQ.V2[k]= *s+ normTmp; *s= -normTmp;
				normTmp= sqrt(2* normTmp* QRQ.V2[k]);}
			else
			{	QRQ.V2[k]= *s- normTmp; *s= normTmp;
				normTmp= -sqrt(-2* normTmp* QRQ.V2[k]);}
			QRQ.V2[k]/= normTmp;				// now, V2 > 0
			for (unsigned i= 0; i< Delta; ++i)
				QRQ.Mat[k+i][k]/= normTmp;		// normalize the reflector

			for (unsigned j= k- 1; j< -1; --j)	// all Householder reflections-------
			{	double sum= QRQ.V2[k]* QRQ.Mat[k+Delta][j];	// compute the inner product
				for (unsigned i= 0; i< Delta; ++i)
					sum+= QRQ.Mat[k+i][k]* QRQ.Mat[k+i][j];

				sum*= 2; QRQ.Mat[k+Delta][j] -= QRQ.V2[k]* sum;	// apply reflector
				for (unsigned i= 0; i< Delta; ++i)
					QRQ.Mat[k+i][j] -= QRQ.Mat[k+i][k]* sum;}}}
	return QRQ;}		/* return the decomposition	*/





Matrix SolveQRQ(QResult QRQ, Matrix b)
{	assert((QRQ.Mat.rows() == b.rows()) && "SolveQRQ. Dimensions do not match.");
	Matrix xRes= Matrix(QRQ.Mat.cols(), b.cols());// the sizes of x and b differ!
	for (unsigned c= 0; c< b.cols(); c++)		// cover all columns of b
	{	for (unsigned i= 0; i< QRQ.Permutation.Length(); i++)
			if (i != QRQ.Permutation[i])			// Step 1: Permutation P
				std::swap(b[i][c], b[QRQ.Permutation[i]][c]);

		unsigned Delta= QRQ.Mat.rows()- QRQ.Rank;
		if (Delta)	// Householder Q2, if applicable
			for (unsigned k= QRQ.Rank- 1; k< -1; --k)	// apply all Householders
			{	double sum= QRQ.V2[k]* b[k+Delta][c];	// inner product Q2
				for (unsigned i= 0; i< Delta; ++i)
					sum+= QRQ.Mat[k+i][k]* b[k+i][c];
				sum*= 2;
				b[k+Delta][c]-= QRQ.V2[k]* sum;	// apply Q2
				for (unsigned i= 0; i< Delta; ++i)
					b[k+i][c]-= QRQ.Mat[k+i][k]* sum;}

		for (unsigned i= 0; i< xRes.rows(); i++)	// Step 4: L^{-1}
			if (i < QRQ.Rank)				// compute xRes[i]
			{	double *s= &xRes[i][c]; *s= b[i+Delta][c];
				for (unsigned j= 0; j< i; j++)	// substitute backwards
					*s -= QRQ.Mat[i+Delta][j]* xRes[j][c];
				*s/= QRQ.Mat[i+Delta][i];}
			else
				xRes[i][c]= 0;

		for (unsigned k= QRQ.Rank- 1; k< -1 ; k--)	// Step 5: all Householders, Q1
		{	double sum= QRQ.V1[k]* xRes[k][c];	// the inner product
			for (unsigned j= k+1; j< QRQ.Mat.cols(); j++)
				sum+= QRQ.Mat[k][j]* xRes[j][c];

			sum*= 2; xRes[k][c] -= sum* QRQ.V1[k];	// apply Q1
			for (unsigned j= k+1; j< QRQ.Mat.cols(); j++)
				xRes[j][c] -= sum* QRQ.Mat[k][j];}}
	return xRes;}


Matrix Inverse(Matrix A, double aTol= 1e-6)	//	pseudoinverse
{	return SolveQRQ(decomposeQRQ(A, aTol), Eye(A.rows()));}
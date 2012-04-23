#include "stdafx.h"
#include "matrix_reshape.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
void testMatrixReshape() {
	using namespace boost::numeric::ublas;
	typedef matrix<double,column_major> matrix_type;

	matrix_type test (4,5);
	matrix_type matprod (2,2);

	// Reshape a 4 x 5 matrix and acess it as 10 x 2 matrix
	matrix_reshape < matrix_type > add1(test, 10, 2);

	// It is OK to initialize here since no data was copied 
	// in previous step
	for (std::size_t i = 0; i < test.size1 (); ++ i)
		for(std::size_t j =0; j < test.size2 (); ++j) test (i,j) = (i+1);

	// Multiply 2 x 10 matrix by 10 x 2 matrix
	matprod.assign(prod <matrix_type> ( matrix_reshape<matrix_type> (test,2,10),
		matrix_reshape<matrix_type> (test,10,2)));

	std::cout << "matrix" << test << std::endl;
	std::cout << "matprod" << matprod << std::endl;

	// Acess each element of the reshaped matrix with iterators
	typedef matrix_reshape <matrix_type > :: iterator1 i1_t;
	typedef matrix_reshape <matrix_type > :: iterator2 i2_t;
	for (i1_t i1 = add1.begin1(); i1 != add1.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			std::cout << "(" << i2.index1() << "," << i2.index2()
			<< ":" << *i2 << ")  ";
		std::cout << std::endl;
	}
	std::cout << "Test for matrix reshape complete" << std::endl;

}

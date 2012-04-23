#include "stdafx.h"
#include "vector_reshape.hpp"
	#include <boost/numeric/ublas/matrix.hpp>
	#include <boost/numeric/ublas/vector.hpp>
	#include <boost/numeric/ublas/io.hpp>
	
void testVectorReshape () {
	  using namespace boost::numeric::ublas;
	  typedef matrix<double,column_major> matrix_type;

	  vector<double> v (10);
	  matrix_type matprod (2,2);

	  // Adapt a 1 x 10 vector as 2 x 5 matrix
	  vector_reshape < vector <double> > add1(v, 2, 5);

	  // It is OK to initialize here since the adaptor
	  // doees not copy data
	  for (std::size_t i = 0; i < v.size (); ++ i)   v (i) = (i+1)+0.00;

	  // Multiply 2 x 5 matrix with 5 x 2 matrix
	  matprod.assign(prod <matrix_type> ( vec_resh (v,2,5), vec_resh (v,5,2)));

	  std::cout << "vector" << v << std::endl;
	  std::cout << "matprod" << matprod << std::endl;

	  // Acess each element using iterators
	  typedef vector_reshape <vector <double> > :: iterator1 i1_t;
	  typedef vector_reshape <vector <double> > :: iterator2 i2_t;
	  for (i1_t i1 = add1.begin1(); i1 != add1.end1(); ++i1) {
		   for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
				   std::cout << "(" << i2.index1() << "," << i2.index2()
								<< ":" << *i2 << ")  ";
			std::cout << std::endl;
	  }
			  std::cout << "Test for vector reshape complete" << std::endl;

}

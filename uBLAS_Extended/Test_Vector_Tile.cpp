#include "stdafx.h"
#include "tile_vector.hpp"
	#include <boost/numeric/ublas/matrix.hpp>
	#include <boost/numeric/ublas/vector.hpp>
	#include <boost/numeric/ublas/io.hpp>

	void testVectorTile () {
	  using namespace boost::numeric::ublas;
	  matrix<double,column_major> matprod (1,10);
	  vector<double> v (10);

	  // Create 10 x 100 matrix tiling a 1x10 vector
	  vector_tile < vector <double> > add1(v, 10,10);

	  for (std::size_t i = 0; i < v.size (); ++ i)   v (i) = (i+1);
	  // Multiplies a ( 10 x 100 ) matrix by (100 x 10) matrix
	  matprod.assign(prod <matrix<double,column_major> > (matvec(v,1,10),matvec(v,100,1)));

	  std::cout << "vector" << v << std::endl;
	  std::cout << "matprod" << matprod << std::endl;

	  // Access every element of a tile using iterators
	  typedef vector_tile <vector <double> > :: iterator1 i1_t;
	  typedef vector_tile <vector <double> > :: iterator2 i2_t;
	  for (i1_t i1 = add1.begin1(); i1 != add1.end1(); ++i1) {
		for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			  std::cout << "(" << i2.index1() << "," << i2.index2()
						<< ":" << *i2 << ")  ";
			std::cout << std::endl;
	  }
	  std::cout << "Test for Vector tile complete" << std::endl;
}

	#include "stdafx.h"	
	#include "tile_matrix.hpp"
	#include <boost/numeric/ublas/matrix.hpp>
	#include <boost/numeric/ublas/io.hpp>

	void testMatrixTile() {
	  using namespace boost::numeric::ublas;
	  typedef matrix<double,column_major> matrix_type;

	  matrix_type m(10,10);
	  matrix_type matprod (10,10);
	  //Tile a 10 x 10 matrix to form a 40 x 40 matrix
	  matrix_tile < matrix_type > add1(m, 4, 4);

	  // No data was copied during tiling so it is
	  // OK to initialize the matrix being tiled here!
	  for (std::size_t i = 0; i < m. size1 () ; ++i)
		  for (std::size_t j =0; j < m.size2 () ; ++j)
			  m (i,j) = i + j;

	  // Multiply 10 x 20 matrix by 20 x 10 matrix
	  matprod.assign(prod <matrix_type> (matmat(m,1,2),matmat(m,2,1)));

	  std::cout << "matrix" << m << std::endl;
	  std::cout << "matprod" << matprod << std::endl;

	  //Access each element of the tile by iterators
	  typedef matrix_tile < matrix_type > :: iterator1 i1_t;
	  typedef matrix_tile < matrix_type > :: iterator2 i2_t;
		for (i1_t i1 = add1.begin1(); i1 != add1.end1(); ++i1) {
		   for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
			  std::cout << "(" << i2.index1() << "," << i2.index2()
							<< ":" << *i2 << ")";
		std::cout << std::endl;
	  }
			  std::cout << "Test for matrix tile complete" << std::endl;

	}

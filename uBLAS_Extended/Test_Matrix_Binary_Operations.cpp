	#include "stdafx.h"
	#include "matrix_expression_ext.hpp"
	#include <boost/numeric/ublas/matrix.hpp>
	#include <boost/numeric/ublas/io.hpp>
	void testMatrixBinaryOperations () {
		using namespace boost::numeric::ublas;
		matrix<double> m (3, 3), mm(3,3);
		matrix<bool> mb1(3,3), mb2(3,3);
		for (unsigned i = 0; i < m.size1 (); ++ i){
			for (unsigned j = 0; j < m.size2 (); ++ j){
				m (i, j) = (pow(-1.0,(int)j))*(i+j);
				mm(i, j) = (i+j);
				mb1 (i,j)= true;
				mb2 (i,j)= false;
			}
		}
		std::cout << 1 + m << std::endl;
		std::cout << m + 1 << std::endl;
		std::cout << 1 - m << std::endl;
		std::cout << m - 1 << std::endl;
		std::cout << (m ^ 2) << std::endl;
		std::cout << (true && mb1) << std::endl;
		std::cout << (mb1 && true) << std::endl;
		std::cout << (mb1 && mb2) << std::endl;
		std::cout << (true || mb2) << std::endl;
		std::cout << (mb2 || true) << std::endl;
		std::cout << (mb1 || mb2) << std::endl;
		std::cout << (!mb1) << std::endl;
		std::cout << (m > 3) << std::endl;
		std::cout << (m >= 3) << std::endl;
		std::cout << (m > mm) << std::endl;
		std::cout << (5 > m) << std::endl;
		std::cout << (5 >= m) << std::endl;
		std::cout << (m >= mm) << std::endl;
		std::cout << (m < 3) << std::endl;
		std::cout << (3 < m) << std::endl;
		std::cout << (m < mm) << std::endl;
		std::cout << (3 <= m) << std::endl;
		std::cout << (m <= 3) << std::endl;
		std::cout << (m <= mm) << std::endl;
		std::cout << (3 != m) << std::endl;
		std::cout << (m != 3) << std::endl;
		std::cout << (m != mm) << std::endl;
		std::cout << ((m <= 3) && (m > mm)) << std::endl;
		std::cout << " Testing of Matrix binary operations complete " << std::endl;
}

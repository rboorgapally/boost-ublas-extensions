#include "stdafx.h"
#include "matrix_expression_ext.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
void testMatrixUnaryOperations () {
	using namespace boost::numeric::ublas;
	matrix<double> m (3, 3), mm(3,3);
	for (unsigned i = 0; i < m.size1 (); ++ i){
		for (unsigned j = 0; j < m.size2 (); ++ j){
			m (i, j) = (pow(-1.0f,(int)j))*(i+j);
			mm(i, j) = (i+j + 1);
		}
	}
	std::cout << abs(m) << std::endl;
	std::cout << log10(mm) << std::endl;
	std::cout << exp(m) << std::endl;
	std::cout << log(mm) << std::endl;
	std::cout << sin(m) << std::endl;
	std::cout << cos(m) << std::endl;
	std::cout << sqrt(mm) << std::endl;
	std::cout << sum(m) << std::endl;
	std::cout << norm_L1(m) << std::endl;
	std::cout << norm_L2(m) << std::endl;
	std::cout << col_norm_inf(m) << std::endl;
	std::cout << max(m) << std::endl;
	std::cout << min(m) << std::endl;
	std::cout << abs_max(m) << std::endl;
	std::cout << abs_min(m) << std::endl;
	std::cout << mean(m) << std::endl;
	std::cout << abs_mean(m) << std::endl;
	std::cout << rms (m) << std::endl;
	std::cout << "Test for Matrix unary operators complete" << std::endl;

}

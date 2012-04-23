#include "stdafx.h"
#include "vector_expression_ext.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

void testVectorUnaryReturningScalar () {
	using namespace boost::numeric::ublas;
	vector<double>  v(10);
	for (unsigned i = 0; i < v.size (); ++ i)
		v (i) = (i+1);

	std::cout << log_min (v) << std::endl;
	std::cout << log_max (v) << std::endl;
	std::cout << min(v) << std::endl;
	std::cout << max(v) << std::endl;
	std::cout << sum(v) << std::endl;
	std::cout << abs_min(v) << std::endl;
	std::cout << abs_max(v) << std::endl;
	std::cout << mean(v) << std::endl;
	std::cout << abs_mean(v) << std::endl;
	std::cout << rms(v) << std::endl;

	std::cout << " Unary vector returning scalar test complete" << std::endl;
}

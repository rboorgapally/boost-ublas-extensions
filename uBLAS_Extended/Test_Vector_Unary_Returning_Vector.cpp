#include "stdafx.h"
#include "vector_expression_ext.hpp"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

void testVectorUnaryReturningVector () {
	using namespace boost::numeric::ublas;
	vector<double>  v(10);
	for (unsigned i = 0; i < v.size (); ++ i)
		v (i) = (i+1);

	std::cout << abs(v) << std::endl;
	std::cout << exp(v) << std::endl;
	std::cout << log(v) << std::endl;
	std::cout << log10(v) << std::endl;
	std::cout << sqrt(v) << std::endl;
	std::cout << sin(v)  << std::endl;
	std::cout << cos(v)  << std::endl;
	//std::cout << log_min (v) << std::endl;
	//std::cout << log_max (v) << std::endl;
	//std::cout << min(v) << std::endl;
	//std::cout << max(v) << std::endl;
	//std::cout << sum(v) << std::endl;
	//std::cout << abs_min(v) << std::endl;
	//std::cout << abs_max(v) << std::endl;
	//std::cout << mean(v) << std::endl;
	//std::cout << abs_mean(v) << std::endl;
	//std::cout << rms(v) << std::endl;

	std::cout << " Unary vector returning vector test complete" << std::endl;
}
